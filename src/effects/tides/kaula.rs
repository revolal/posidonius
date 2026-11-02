// Implemented by Alexandre Revol alexandre.revol@unige.ch
use super::{TidalModel, TidesEffect};
use crate::constants::{G, MAX_PARTICLES};
use crate::{Axes, Particle};
use serde::{Deserialize, Serialize};

mod love_number;
pub use love_number::LoveNumber;

// Planetary and stellar tides (docs/tides/kaula_tides.pdf, Chapters 8 to 10)
pub mod kaula_tide;

// Structure declare the intput parameters needed for the Kaula model calculation.
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq, Default)]
#[serde(deny_unknown_fields)]
pub struct KaulaParameters {
    // Externally provided spectrum for excitation frequency, real and imaginary part.
    #[serde(default)]
    pub love_numbers: LoveNumber,
    // Secular force of the former Kaula code: no longer used (always zero), kept so that the .json files
    // written by earlier versions remain readable. Still written: the .bin snapshots (bincode) need the
    // same fields when written and read.
    #[serde(default)]
    pub tidal_force: Axes,

    // Tide of this body (kaula_tide, Section 10.4): prepared spectrum, mean-rate history updated at the
    // accepted steps, and last force/torque with the state that produced them. The aliases read the
    // .json files written by development versions of this model.
    #[serde(default, alias = "planet_spectrum")]
    pub spectrum: kaula_tide::TideSpectrum,
    #[serde(default)]
    pub rate_history: kaula_tide::RateHistory,
    #[serde(default, alias = "planet_tide")]
    pub tide_state: TideState,
}

// Last planetary tide computed, keyed by the exact state (r, v, s) that produced it: the torque is
// reused only for the same state, otherwise recomputed (LeapFrog asks for the torque without the force).
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq, Default)]
#[serde(default)]
pub struct TideState {
    pub valid: bool,
    pub position: Axes,
    pub velocity: Axes,
    pub spin: Axes,
    pub tide: kaula_tide::Tide,
}

// Planetary tide of the particle from its tidal coordinates (relative to the tidal host), its spin and
// its mean-rate history: eqs. (8.17), (8.18) and (8.21).
fn calculate_planetary_tide(
    tidal_host_particle: &Particle,
    particle: &Particle,
    spectrum: &kaula_tide::TideSpectrum,
) -> kaula_tide::Tide {
    let kaula = particle.tides.get_kaula();
    let position = particle.tides.coordinates.position;
    let geometry = kaula_tide::nodal_geometry(
        position,
        particle.tides.coordinates.velocity,
        particle.spin,
        G * (tidal_host_particle.mass + particle.mass),
    );
    let rates = kaula_tide::mean_rates(&kaula.rate_history, geometry.n_osc);
    kaula_tide::tide_tabulated(
        position,
        &geometry,
        &rates,
        &kaula.love_numbers,
        spectrum,
        G,
        tidal_host_particle.mass,
        particle.radius,
        kaula_tide::q_max(geometry.e),
    )
}

fn tide_state_matches(particle: &Particle) -> bool {
    let state = &particle.tides.get_kaula().tide_state;
    state.valid
        && state.position == particle.tides.coordinates.position
        && state.velocity == particle.tides.coordinates.velocity
        && state.spin == particle.spin
}

// Stellar tide of the pair (tidal host, orbiting body), stored in the orbiting body (one per body,
// Section 10.4): mean-rate history of the orbit seen in the equator of the host, and
// last force/torque with the state that produced them (body tidal coordinates, host spin).
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq, Default)]
#[serde(default)]
pub struct HostTide {
    pub rate_history: kaula_tide::RateHistory,
    pub tide_state: TideState,
}

// r = r_D - r_P and v for the stellar tide (D = host, P = orbiting body): opposite of the tidal
// coordinates of the body
fn stellar_tide_coordinates(particle: &Particle) -> (Axes, Axes) {
    let mut position = particle.tides.coordinates.position;
    position.negate();
    let mut velocity = particle.tides.coordinates.velocity;
    velocity.negate();
    (position, velocity)
}

// Tide raised in the tidal host by the orbiting body: the equations of the planetary tide with D = host
// and P = body (Section 9.1), signed spectrum of the host,
// mean rates of the pair.
fn calculate_stellar_tide(
    tidal_host_particle: &Particle,
    particle: &Particle,
    spectrum: &kaula_tide::TideSpectrum,
) -> kaula_tide::Tide {
    let kaula = tidal_host_particle.tides.get_kaula();
    let (position, velocity) = stellar_tide_coordinates(particle);
    let geometry = kaula_tide::nodal_geometry(
        position,
        velocity,
        tidal_host_particle.spin,
        G * (tidal_host_particle.mass + particle.mass),
    );
    let rates = kaula_tide::mean_rates(&particle.tides.host_kaula_tide.rate_history, geometry.n_osc);
    kaula_tide::tide_tabulated(
        position,
        &geometry,
        &rates,
        &kaula.love_numbers,
        spectrum,
        G,
        particle.mass,
        tidal_host_particle.radius,
        kaula_tide::q_max(geometry.e),
    )
}

fn host_tide_state_matches(tidal_host_particle: &Particle, particle: &Particle) -> bool {
    let state = &particle.tides.host_kaula_tide.tide_state;
    state.valid
        && state.position == particle.tides.coordinates.position
        && state.velocity == particle.tides.coordinates.velocity
        && state.spin == tidal_host_particle.spin
}

fn host_spectrum(tidal_host_particle: &Particle) -> kaula_tide::TideSpectrum {
    let kaula = tidal_host_particle.tides.get_kaula();
    if kaula.spectrum.prepared {
        kaula.spectrum
    } else {
        kaula_tide::prepare_signed_spectrum(&kaula.love_numbers)
    }
}

/// Force of the orbiting body on the tidal host due to the tide raised in the host (stellar tide,
/// eq. (8.17) with D = host, P = body); the body receives the opposite force.
pub fn calculate_stellar_tidal_force(
    tidal_host_particle: &mut Particle,
    particle: &mut Particle,
) -> Axes {
    let spectrum = host_spectrum(tidal_host_particle);
    tidal_host_particle.tides.get_kaula_mut().spectrum = spectrum;
    let tide = calculate_stellar_tide(tidal_host_particle, particle, &spectrum);
    particle.tides.host_kaula_tide.tide_state = TideState {
        valid: true,
        position: particle.tides.coordinates.position,
        velocity: particle.tides.coordinates.velocity,
        spin: tidal_host_particle.spin,
        tide,
    };
    tide.force
}

/// Dissipated power of the stellar tide of one orbiting body, eq. (8.21) with D = host:
/// P = -(F_Im . (v_host - v_body) + Gamma_Im . s_host), at the state of the last force evaluation.
pub fn stellar_tide_power(particle: &Particle) -> f64 {
    let state = particle.tides.host_kaula_tide.tide_state;
    if !state.valid {
        return 0.;
    }
    // state.velocity = v_body - v_host
    state.tide.dissipative_force.dot(&state.velocity)
        - state.tide.dissipative_torque.dot(&state.spin)
}

/// Update the mean-rate histories once per accepted step (Sections 8.4.1 and 10.4): of every planet
/// with a Kaula tide (orbit seen in its equator) and, if the tidal host has a Kaula tide, of every pair
/// (host, orbiting body) (orbit seen in the equator of the host). The averaging window (8.12) uses the
/// osculating mean motions of all the bodies orbiting the tidal host, with or without tides, computed
/// from the inertial coordinates relative to the host (tides.coordinates is set only for the bodies
/// with tides). The spins are L / I of the accepted step.
pub fn update_rate_histories(
    tidal_host_index: usize,
    particles: &mut [Particle],
    current_time: f64,
) {
    let host = particles[tidal_host_index];
    let relative = |particle: &Particle| {
        let mut r = particle.inertial_position;
        r.sub(&host.inertial_position);
        let mut v = particle.inertial_velocity;
        v.sub(&host.inertial_velocity);
        (r, v)
    };
    let mut mean_motions = [0.; MAX_PARTICLES];
    for (j, particle) in particles.iter().enumerate() {
        if j != tidal_host_index {
            let (r, v) = relative(particle);
            mean_motions[j] = kaula_tide::mean_motion_about_host(
                r,
                v,
                G * (host.mass + particle.mass),
            );
        }
    }
    let n_particles = particles.len();
    let host_kaula = matches!(
        host.tides.effect,
        TidesEffect::CentralBody(TidalModel::Kaula(_))
    );
    let mut host_spin = host.angular_momentum;
    host_spin.div(host.moment_of_inertia);
    for (i, particle) in particles.iter_mut().enumerate() {
        if i == tidal_host_index {
            continue;
        }
        let planet_kaula = matches!(
            particle.tides.effect,
            TidesEffect::OrbitingBody(TidalModel::Kaula(_))
        );
        // Bodies that feel the stellar tide: every OrbitingBody (Section 10.2)
        let stellar = host_kaula && matches!(particle.tides.effect, TidesEffect::OrbitingBody(_));
        if !planet_kaula && !stellar {
            continue;
        }
        let (r, v) = relative(particle);
        let mu = G * (host.mass + particle.mass);
        let mut others = [0.; MAX_PARTICLES];
        for j in 0..n_particles {
            if j != i && j != tidal_host_index {
                others[j] = mean_motions[j];
            }
        }
        let window = kaula_tide::averaging_window(mean_motions[i], &others[..n_particles]);
        if planet_kaula {
            let mut spin = particle.angular_momentum;
            spin.div(particle.moment_of_inertia);
            let geometry = kaula_tide::nodal_geometry(r, v, spin, mu);
            kaula_tide::update_rate_history(
                &mut particle.tides.get_kaula_mut().rate_history,
                current_time,
                &geometry,
                window,
            );
        }
        if stellar {
            let (mut r_host, mut v_host) = (r, v);
            r_host.negate();
            v_host.negate();
            let geometry = kaula_tide::nodal_geometry(r_host, v_host, host_spin, mu);
            kaula_tide::update_rate_history(
                &mut particle.tides.host_kaula_tide.rate_history,
                current_time,
                &geometry,
                window,
            );
        }
    }
}

/// Force of the tidal host on the orbiting body due to the tide raised in the body (planetary tide,
/// eq. (8.17)); the host receives the opposite force.
pub fn calculate_planetary_tidal_force(
    tidal_host_particle: &Particle,
    particle: &mut Particle,
) -> Axes {
    if !particle.tides.get_kaula().spectrum.prepared {
        let spectrum = kaula_tide::prepare_spectrum(&particle.tides.get_kaula().love_numbers);
        particle.tides.get_kaula_mut().spectrum = spectrum;
    }
    let spectrum = particle.tides.get_kaula().spectrum;
    let tide = calculate_planetary_tide(tidal_host_particle, particle, &spectrum);
    let (position, velocity, spin) = (
        particle.tides.coordinates.position,
        particle.tides.coordinates.velocity,
        particle.spin,
    );
    particle.tides.get_kaula_mut().tide_state = TideState {
        valid: true,
        position,
        velocity,
        spin,
        tide,
    };
    // Spherical components are not computed; denergy_dt has its own Kaula branch, eq. (8.21)
    particle
        .tides
        .parameters
        .internal
        .orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 0.;
    tide.force
}

/// Torque of the tide, of the same state as the force (Section 10.2): planetary tide (central_body =
/// false, on the spin of the orbiting body) or stellar tide (central_body = true, on the spin of the
/// host, due to this body). Returns r x F = -Gamma, which tides.rs negates.
pub fn calculate_torque_due_to_tides(
    tidal_host_particle: &Particle,
    particle: &Particle,
    central_body: bool,
) -> Axes {
    let tide = if !central_body {
        if tide_state_matches(particle) {
            particle.tides.get_kaula().tide_state.tide
        } else {
            let kaula = particle.tides.get_kaula();
            let spectrum = if kaula.spectrum.prepared {
                kaula.spectrum
            } else {
                kaula_tide::prepare_spectrum(&kaula.love_numbers)
            };
            calculate_planetary_tide(tidal_host_particle, particle, &spectrum)
        }
    } else {
        // No stellar tide for the bodies without tidal coordinates (TidesEffect::Disabled)
        if !matches!(particle.tides.effect, TidesEffect::OrbitingBody(_)) {
            return Axes::new();
        }
        if host_tide_state_matches(tidal_host_particle, particle) {
            particle.tides.host_kaula_tide.tide_state.tide
        } else {
            calculate_stellar_tide(tidal_host_particle, particle, &host_spectrum(tidal_host_particle))
        }
    };
    let mut torque = tide.torque;
    torque.negate();
    torque
}
