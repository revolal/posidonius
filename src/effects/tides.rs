use std::collections::HashMap;
use super::super::tools;
use super::super::constants::{K2};
use super::super::{Particle};
use super::super::{Axes};
use super::{EvolutionType};

use crate::constants::{G, TWO_PI, DAY, DEG2RAD};//, AU, M_SUN, M2AU,R_SUN,R_EARTH};//, G_SI, AU, M_SUN};

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
//#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TidesParticleInputParameters {
    pub dissipation_factor: f64,
    pub dissipation_factor_scale: f64, // to scale the dissipation factor (multiply)
    pub love_number: f64,   // Love number of degree 2 (i.e., k2). Dimensionless parameters that measure the rigidity of a planetary body and the 
                            // susceptibility of its shape to change in response to a tidal potential.
    pub kaula_coplanar_tides_input_parameters: KaulaCoplanarTidesInputParticleParameters,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
//#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct KaulaCoplanarTidesInputParticleParameters {
    pub love_number_excitation_frequency: [[f64;32];32],
    pub real_part_love_number: [[f64;32];32],
    pub imaginary_part_love_number: [[f64;32];32],
    pub num_datapoints: f64,
}
///////////////////////////////////////////////////////////////////////////////////////////////END MODIFICATION

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
//#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TidesParticleInternalParameters {
    pub distance: f64,
    pub radial_velocity: f64,
    pub scaled_dissipation_factor: f64, // sigma (dissipation_factor_scale*dissipation_factor)
    pub scalar_product_of_vector_position_with_stellar_spin: f64,
    pub scalar_product_of_vector_position_with_planetary_spin: f64,
    pub orthogonal_component_of_the_tidal_force_due_to_stellar_tide: f64,
    pub orthogonal_component_of_the_tidal_force_due_to_planetary_tide: f64,
    pub radial_component_of_the_tidal_force: f64,
    pub radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass: f64, // Needed to compute denergy_dt
    pub denergy_dt: f64, // Only for history output
    pub lag_angle: f64, // Used by EvolutionType::BolmontMathis2016, EvolutionType::GalletBolmont2017 and EvolutionType::LeconteChabrier2013(true)

    pub a: f64,
    pub spin: f64,
    pub orbital_frequency: f64,

    pub tidal_torque_z: f64,

    pub sigma220_2_excitative_frequency: f64,
    pub sigma220_1_excitative_frequency: f64,
    pub sigma2200_excitative_frequency: f64,
    pub sigma2201_excitative_frequency: f64,
    pub sigma2202_excitative_frequency: f64,
    pub im_love_number_sigma220_2: f64,
    pub im_love_number_sigma220_1: f64,
    pub im_love_number_sigma2200: f64,
    pub im_love_number_sigma2201: f64,
    pub im_love_number_sigma2202: f64,
    pub re_love_number_sigma220_2: f64,
    pub re_love_number_sigma220_1: f64,
    pub re_love_number_sigma2200: f64,
    pub re_love_number_sigma2201: f64,
    pub re_love_number_sigma2202: f64,

    pub sigma201_2_excitative_frequency: f64,
    pub sigma201_1_excitative_frequency: f64,
    pub sigma2010_excitative_frequency: f64,
    pub sigma2011_excitative_frequency: f64,
    pub sigma2012_excitative_frequency: f64,
    pub im_love_number_sigma201_2: f64,
    pub im_love_number_sigma201_1: f64,
    pub im_love_number_sigma2010: f64,
    pub im_love_number_sigma2011: f64,
    pub im_love_number_sigma2012: f64,
    pub re_love_number_sigma201_2: f64,
    pub re_love_number_sigma201_1: f64,
    pub re_love_number_sigma2010: f64,
    pub re_love_number_sigma2011: f64,
    pub re_love_number_sigma2012: f64,

    // pub im_love_number:f64,
    // pub re_love_number:f64,

    pub check_excitative_frequ: f64,

    pub excitative_frequ_sigma220q: [f64;5],
    pub re_love_number_sigma220q: [f64;5],
    pub im_love_number_sigma220q: [f64;5],

    pub excitative_frequ_sigma201q: [f64;5],
    pub re_love_number_sigma201q: [f64;5],
    pub im_love_number_sigma201q: [f64;5],
    
    pub force_by_tides: Axes, 
    pub dangular_momentum_dt_by_tides: Axes, 
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct TidesParticleOutputParameters {
    pub acceleration: Axes,
    pub dangular_momentum_dt: Axes, // Force
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
//#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TidesParticleParameters {
    pub input: TidesParticleInputParameters,
    //pub input: TidesTypes,
    pub internal: TidesParticleInternalParameters,
    pub output: TidesParticleOutputParameters,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct TidesParticleCoordinates {
    // Positions/velocities in a heliocentric frame 
    // (i.e., the host is at rest with respect to the origin of the coordinate system)
    pub position: Axes,
    pub velocity: Axes,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum TidesEffect {
    CentralBody,
    OrbitingBody,
    KaulaCoplanarOrbitingBody,
    //ConstTimeLagOrbitingBody,
    Disabled,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
//#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Tides {
    pub effect: TidesEffect,
    pub parameters: TidesParticleParameters,
    pub coordinates: TidesParticleCoordinates,
}

//impl TidesEffect{.....? prendre en compte l'enum dans 
impl Tides {
    //pub fn new(effect: TidesEffect, dissipation_factor: f64, dissipation_factor_scale: f64, love_number: f64) -> Tides {
    pub fn new(effect: TidesEffect, dissipation_factor: f64, dissipation_factor_scale: f64, love_number: f64, love_number_excitation_frequency: [[f64;32];32], real_part_love_number: [[f64;32];32], imaginary_part_love_number: [[f64;32];32], num_datapoints: f64) -> Tides {
    //pub fn new(effect: TidesEffect, dissipation_factor: f64, dissipation_factor_scale: f64, love_number: f64, love_number_excitation_frequency: LinkedList<f64>, real_part_love_number: LinkedList<f64>, imaginary_part_love_number: LinkedList<f64>, num_datapoints: f64) -> Tides {
    //pub fn new(effect: TidesEffect, dissipation_factor: f64, dissipation_factor_scale: f64, love_number: f64, love_number_excitation_frequency: f64, real_part_love_number: f64, imaginary_part_love_number: f64, num_datapoints: f64) -> Tides {
        Tides {
            effect: effect,
            parameters: TidesParticleParameters {
                input: TidesParticleInputParameters {
                    dissipation_factor: dissipation_factor,
                    dissipation_factor_scale: dissipation_factor_scale,
                    love_number: love_number,
                    // ////////////////////////////////////////////////////////////////////////////////
                    kaula_coplanar_tides_input_parameters: KaulaCoplanarTidesInputParticleParameters {
                        love_number_excitation_frequency: love_number_excitation_frequency,
                        real_part_love_number: real_part_love_number,
                        imaginary_part_love_number: imaginary_part_love_number,
                        num_datapoints: num_datapoints,
                    },
                    // ////////////////////////////////////////////////////////////////////////////////
                },
                internal: TidesParticleInternalParameters {
                    distance: 0.,
                    radial_velocity: 0.,
                    scaled_dissipation_factor: dissipation_factor_scale*dissipation_factor,
                    scalar_product_of_vector_position_with_stellar_spin: 0.,
                    scalar_product_of_vector_position_with_planetary_spin: 0.,
                    orthogonal_component_of_the_tidal_force_due_to_stellar_tide: 0.,
                    orthogonal_component_of_the_tidal_force_due_to_planetary_tide: 0.,
                    radial_component_of_the_tidal_force: 0.,
                    radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass: 0.,
                    denergy_dt: 0., // Only for history output
                    lag_angle: 0., // It will be initialized the first time the evolver is called
                
                    a: 0.,
                    spin: 0.,
                    orbital_frequency: 0.,

                    tidal_torque_z: 0.,

                    sigma220_2_excitative_frequency:0.,
                    sigma220_1_excitative_frequency:0.,
                    sigma2200_excitative_frequency:0.,
                    sigma2201_excitative_frequency:0.,
                    sigma2202_excitative_frequency:0.,
                    im_love_number_sigma220_2: 0.,
                    im_love_number_sigma220_1: 0.,
                    im_love_number_sigma2200: 0.,
                    im_love_number_sigma2201: 0.,
                    im_love_number_sigma2202: 0.,
                    re_love_number_sigma220_2: 0.,
                    re_love_number_sigma220_1: 0.,
                    re_love_number_sigma2200: 0.,
                    re_love_number_sigma2201: 0.,
                    re_love_number_sigma2202: 0.,

                    sigma201_2_excitative_frequency:0.,
                    sigma201_1_excitative_frequency:0.,
                    sigma2010_excitative_frequency:0.,
                    sigma2011_excitative_frequency:0.,
                    sigma2012_excitative_frequency:0.,
                    im_love_number_sigma201_2: 0.,
                    im_love_number_sigma201_1: 0.,
                    im_love_number_sigma2010: 0.,
                    im_love_number_sigma2011: 0.,
                    im_love_number_sigma2012: 0.,
                    re_love_number_sigma201_2: 0.,
                    re_love_number_sigma201_1: 0.,
                    re_love_number_sigma2010: 0.,
                    re_love_number_sigma2011: 0.,
                    re_love_number_sigma2012: 0.,

                    // im_love_number:0.,
                    // re_love_number:0.,

                    check_excitative_frequ:0., // if check/sigma.abs() > 10% -> recalcul   Calcul check with derive K2

                    excitative_frequ_sigma220q: [0.;5],
                    re_love_number_sigma220q: [0.;5],
                    im_love_number_sigma220q: [0.;5],

                    excitative_frequ_sigma201q: [0.;5],
                    re_love_number_sigma201q: [0.;5],
                    im_love_number_sigma201q: [0.;5],

                    force_by_tides: Axes{x: 0., y: 0., z: 0.},
                    dangular_momentum_dt_by_tides: Axes{x: 0., y: 0., z: 0.},
                },
                output: TidesParticleOutputParameters {
                    acceleration: Axes{x: 0., y: 0., z: 0.},
                    dangular_momentum_dt: Axes{x: 0., y: 0., z: 0.},
                },
            },
            coordinates: TidesParticleCoordinates {
                position: Axes{x: 0., y: 0., z: 0.},
                velocity: Axes{x: 0., y: 0., z: 0.},
            },
        }
    }
}
/////////////////////////:

pub fn initialize(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    if let TidesEffect::CentralBody = host_particle.tides.effect {
        host_particle.tides.parameters.internal.scalar_product_of_vector_position_with_stellar_spin = 0.;
        host_particle.tides.parameters.internal.scalar_product_of_vector_position_with_planetary_spin = 0.;
        host_particle.tides.parameters.output.acceleration.x = 0.;
        host_particle.tides.parameters.output.acceleration.y = 0.;
        host_particle.tides.parameters.output.acceleration.z = 0.;
        host_particle.tides.parameters.output.dangular_momentum_dt.x = 0.;
        host_particle.tides.parameters.output.dangular_momentum_dt.y = 0.;
        host_particle.tides.parameters.output.dangular_momentum_dt.z = 0.;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let TidesEffect::OrbitingBody = particle.tides.effect {
            //if let TidesEffect::KaulaCoplanarOrbitingBody = particle.tides.effect {
            //if let TidesEffect::CTOrbitingBody = particle.tides.effect {
                particle.tides.parameters.internal.scalar_product_of_vector_position_with_stellar_spin = particle.tides.coordinates.position.x * host_particle.spin.x 
                                + particle.tides.coordinates.position.y * host_particle.spin.y
                                + particle.tides.coordinates.position.z * host_particle.spin.z;
                particle.tides.parameters.internal.scalar_product_of_vector_position_with_planetary_spin = particle.tides.coordinates.position.x * particle.spin.x 
                                + particle.tides.coordinates.position.y * particle.spin.y
                                + particle.tides.coordinates.position.z * particle.spin.z;
                particle.tides.parameters.output.acceleration.x = 0.;
                particle.tides.parameters.output.acceleration.y = 0.;
                particle.tides.parameters.output.acceleration.z = 0.;
                particle.tides.parameters.output.dangular_momentum_dt.x = 0.;
                particle.tides.parameters.output.dangular_momentum_dt.y = 0.;
                particle.tides.parameters.output.dangular_momentum_dt.z = 0.;
            }
            if let TidesEffect::KaulaCoplanarOrbitingBody = particle.tides.effect {
                    particle.tides.parameters.internal.scalar_product_of_vector_position_with_stellar_spin = particle.tides.coordinates.position.x * host_particle.spin.x 
                                    + particle.tides.coordinates.position.y * host_particle.spin.y
                                    + particle.tides.coordinates.position.z * host_particle.spin.z;
                    particle.tides.parameters.internal.scalar_product_of_vector_position_with_planetary_spin = particle.tides.coordinates.position.x * particle.spin.x 
                                    + particle.tides.coordinates.position.y * particle.spin.y
                                    + particle.tides.coordinates.position.z * particle.spin.z;
                    particle.tides.parameters.output.acceleration.x = 0.;
                    particle.tides.parameters.output.acceleration.y = 0.;
                    particle.tides.parameters.output.acceleration.z = 0.;
                    particle.tides.parameters.output.dangular_momentum_dt.x = 0.;
                    particle.tides.parameters.output.dangular_momentum_dt.y = 0.;
                    particle.tides.parameters.output.dangular_momentum_dt.z = 0.;
                }
        }
    }
}

pub fn inertial_to_heliocentric_coordinates(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    if let TidesEffect::CentralBody = host_particle.tides.effect {
        // Inertial to Heliocentric positions/velocities
        host_particle.tides.coordinates.position.x = 0.;
        host_particle.tides.coordinates.position.y = 0.;
        host_particle.tides.coordinates.position.z = 0.;
        host_particle.tides.coordinates.velocity.x = 0.;
        host_particle.tides.coordinates.velocity.y = 0.;
        host_particle.tides.coordinates.velocity.z = 0.;
        host_particle.tides.parameters.internal.distance = 0.;
        host_particle.tides.parameters.internal.radial_velocity = 0.;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let TidesEffect::OrbitingBody = particle.tides.effect {
            //if let Tides::CTLOrbitingBody = particle.tides.effect {
                //if let TidesEffect::KaulaCoplanarOrbitingBody = particle.tides.effect {
                particle.tides.coordinates.position.x = particle.inertial_position.x - host_particle.inertial_position.x;
                particle.tides.coordinates.position.y = particle.inertial_position.y - host_particle.inertial_position.y;
                particle.tides.coordinates.position.z = particle.inertial_position.z - host_particle.inertial_position.z;
                particle.tides.coordinates.velocity.x = particle.inertial_velocity.x - host_particle.inertial_velocity.x;
                particle.tides.coordinates.velocity.y = particle.inertial_velocity.y - host_particle.inertial_velocity.y;
                particle.tides.coordinates.velocity.z = particle.inertial_velocity.z - host_particle.inertial_velocity.z;
                particle.tides.parameters.internal.distance = (particle.tides.coordinates.position.x.powi(2) 
                                                                + particle.tides.coordinates.position.y.powi(2)
                                                                + particle.tides.coordinates.position.z.powi(2)).sqrt();
                particle.tides.parameters.internal.radial_velocity = (particle.tides.coordinates.position.x*particle.tides.coordinates.velocity.x +
                                                                        particle.tides.coordinates.position.y*particle.tides.coordinates.velocity.y +
                                                                        particle.tides.coordinates.position.z*particle.tides.coordinates.velocity.z) / particle.tides.parameters.internal.distance;
            }
            if let TidesEffect::KaulaCoplanarOrbitingBody = particle.tides.effect {
                    particle.tides.coordinates.position.x = particle.inertial_position.x - host_particle.inertial_position.x;
                    particle.tides.coordinates.position.y = particle.inertial_position.y - host_particle.inertial_position.y;
                    particle.tides.coordinates.position.z = particle.inertial_position.z - host_particle.inertial_position.z;
                    particle.tides.coordinates.velocity.x = particle.inertial_velocity.x - host_particle.inertial_velocity.x;
                    particle.tides.coordinates.velocity.y = particle.inertial_velocity.y - host_particle.inertial_velocity.y;
                    particle.tides.coordinates.velocity.z = particle.inertial_velocity.z - host_particle.inertial_velocity.z;
                    particle.tides.parameters.internal.distance = (particle.tides.coordinates.position.x.powi(2) 
                                                                    + particle.tides.coordinates.position.y.powi(2)
                                                                    + particle.tides.coordinates.position.z.powi(2)).sqrt();
                    particle.tides.parameters.internal.radial_velocity = (particle.tides.coordinates.position.x*particle.tides.coordinates.velocity.x +
                                                                            particle.tides.coordinates.position.y*particle.tides.coordinates.velocity.y +
                                                                            particle.tides.coordinates.position.z*particle.tides.coordinates.velocity.z) / particle.tides.parameters.internal.distance;
                }
        }
    }
}

pub fn copy_heliocentric_coordinates(host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    if let TidesEffect::CentralBody = host_particle.tides.effect {
        host_particle.tides.coordinates.position = host_particle.heliocentric_position;
        host_particle.tides.coordinates.velocity = host_particle.heliocentric_velocity;
        host_particle.tides.parameters.internal.distance = host_particle.heliocentric_distance;
        host_particle.tides.parameters.internal.radial_velocity = host_particle.heliocentric_radial_velocity;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let TidesEffect::OrbitingBody = particle.tides.effect {
            //if let Tides::CTLOrbitingBody = particle.tides.effect {
            //if let TidesEffect::KaulaCoplanarOrbitingBody = particle.tides.effect {
                particle.tides.coordinates.position = particle.heliocentric_position;
                particle.tides.coordinates.velocity = particle.heliocentric_velocity;
                particle.tides.parameters.internal.distance = particle.heliocentric_distance;
                particle.tides.parameters.internal.radial_velocity = particle.heliocentric_radial_velocity;
            }
            if let TidesEffect::KaulaCoplanarOrbitingBody = particle.tides.effect {
                    particle.tides.coordinates.position = particle.heliocentric_position;
                    particle.tides.coordinates.velocity = particle.heliocentric_velocity;
                    particle.tides.parameters.internal.distance = particle.heliocentric_distance;
                    particle.tides.parameters.internal.radial_velocity = particle.heliocentric_radial_velocity;
                }
        }
    }
}


pub fn calculate_planet_dependent_dissipation_factors(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>) {
    match tidal_host_particle.evolution {
        EvolutionType::BolmontMathis2016(_) | EvolutionType::GalletBolmont2017(_) | EvolutionType::LeconteChabrier2013(true) => {
            let star_norm_spin_vector = tidal_host_particle.norm_spin_vector_2.sqrt();
            for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
                //
                //// Excitation frequency needed by the model based on the
                // instantaneous frequency (using positions, velocities and spins)
                //let frequency = (particle.tides.coordinates.velocity.x - tidal_host_particle.spin.y*particle.tides.coordinates.position.z + tidal_host_particle.spin.z*particle.tides.coordinates.position.y).powi(2)
                            //+ (particle.tides.coordinates.velocity.y - tidal_host_particle.spin.z*particle.tides.coordinates.position.x + tidal_host_particle.spin.x*particle.tides.coordinates.position.z).powi(2)
                            //+ (particle.tides.coordinates.velocity.z - tidal_host_particle.spin.x*particle.tides.coordinates.position.y + tidal_host_particle.spin.y*particle.tides.coordinates.position.x).powi(2);
                //let inverse_of_half_the_excitation_frequency = particle.tides.parameters.internal.distance / frequency;
                // NOTE:  two_times_the_inverse_of_the_excitation_frequency: 2/w
                //        inverse_of_half_the_excitation_frequency : 1/(w/2)
                //
                //// Excitation frequency needed by the model based on the
                // mean frequency (using mean motion and spin). 
                //
                // NOTE: The model is already here being used outside the 
                // validity domain, it seems not justified to use an 
                // instantaneous frequency.
                let gm = tidal_host_particle.mass_g+particle.mass_g;
                let (perihelion_distance, eccentricity) = tools::calculate_perihelion_distance_and_eccentricity(gm, particle.tides.coordinates.position, particle.tides.coordinates.velocity);
                let mean_motion = gm.sqrt() * (perihelion_distance/(1.0 - eccentricity)).powf(-1.5);
                let half_the_excitation_frequency = (star_norm_spin_vector - mean_motion).abs();
                let inverse_of_half_the_excitation_frequency = 1./half_the_excitation_frequency;

                let planet_dependent_dissipation_factor = tidal_host_particle.tides.parameters.input.dissipation_factor_scale * 2.0 * K2
                    * tidal_host_particle.tides.parameters.internal.lag_angle * inverse_of_half_the_excitation_frequency / (3.0*tidal_host_particle.radius.powi(5));

                star_planet_dependent_dissipation_factors.insert(particle.id.clone(), planet_dependent_dissipation_factor);
                ////////println!("Insert {} in {}", planet_dependent_dissipation_factor, particle.id);
            }
            //panic!("Please, contact Posidonius authors before using BolmontMathis2016/GalletBolmont2017/LeconteChabrier2013(true) evolutionary models. They may not be ready yet for scientific explotation.")
        },
        _ => {},
    }

}

pub fn planet_dependent_dissipation_factor(star_planet_dependent_dissipation_factors: &HashMap<usize, f64>,  id: &usize, evolution: EvolutionType, scaled_dissipation_factor: f64) -> f64 {
    match evolution {
        EvolutionType::BolmontMathis2016(_) | EvolutionType::GalletBolmont2017(_) | EvolutionType::LeconteChabrier2013(true) => {
            match star_planet_dependent_dissipation_factors.get(id) {
                Some(&value) => value,
                _ => scaled_dissipation_factor // This should not happen
            }
        },
        _ => scaled_dissipation_factor,
    }
}


//////////////////////////////////////////////////////////////////////////////
//// TIDES
pub fn calculate_torque_due_to_tides(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], central_body:bool, current_time:f64) {
    let mut dangular_momentum_dt = Axes{x: 0., y: 0., z:0.};
    let mut reference_spin = tidal_host_particle.spin.clone();
    let mut orthogonal_component_of_the_tidal_force: f64;
    let mut reference_rscalspin: f64;
    //panic!("Does it reach this point?-calculate_torque_due_to_tides");
    //posidoinus --silent 

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody = particle.tides.effect {
            // //////println!("Troque due to tides, OrbitingBody");
            //if let Tides::CTLOrbitingBody = particle.tides.effect {
            //if let TidesEffect::KaulaCoplanarOrbitingBody = particle.tides.effect {
            if !central_body {
                reference_spin = particle.spin.clone();
                reference_rscalspin = particle.tides.parameters.internal.scalar_product_of_vector_position_with_planetary_spin;
                orthogonal_component_of_the_tidal_force = particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide;
            } else {
                reference_rscalspin = particle.tides.parameters.internal.scalar_product_of_vector_position_with_stellar_spin;
                orthogonal_component_of_the_tidal_force = particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide;
            }

            // distance to star
            let distance = particle.tides.parameters.internal.distance;

            //// Torque calculation (star)
            // - Equation 8-9 from Bolmont et al. 2015
            let torque_due_to_tides_x: f64 = orthogonal_component_of_the_tidal_force 
                            * (distance * reference_spin.x - reference_rscalspin*particle.tides.coordinates.position.x/distance - 1.0/distance
                            * (particle.tides.coordinates.position.y*particle.tides.coordinates.velocity.z - particle.tides.coordinates.position.z*particle.tides.coordinates.velocity.y) );

            let torque_due_to_tides_y :f64 = orthogonal_component_of_the_tidal_force 
                            * (distance * reference_spin.y - reference_rscalspin*particle.tides.coordinates.position.y/distance - 1.0/distance
                            * (particle.tides.coordinates.position.z*particle.tides.coordinates.velocity.x - particle.tides.coordinates.position.x*particle.tides.coordinates.velocity.z) );

            let torque_due_to_tidez_z: f64 = orthogonal_component_of_the_tidal_force 
                            * (distance * reference_spin.z - reference_rscalspin*particle.tides.coordinates.position.z/distance - 1.0/distance
                            * (particle.tides.coordinates.position.x*particle.tides.coordinates.velocity.y - particle.tides.coordinates.position.y*particle.tides.coordinates.velocity.x) );

            let factor = -1.0;
            if central_body {
                // Integration of the spin (total torque tides):
                dangular_momentum_dt.x += factor * torque_due_to_tides_x;
                dangular_momentum_dt.y += factor * torque_due_to_tides_y;
                dangular_momentum_dt.z += factor * torque_due_to_tidez_z;
            } else {
                particle.tides.parameters.output.dangular_momentum_dt.x = factor * torque_due_to_tides_x;
                particle.tides.parameters.output.dangular_momentum_dt.y = factor * torque_due_to_tides_y;
                particle.tides.parameters.output.dangular_momentum_dt.z = factor * torque_due_to_tidez_z;
                //println!("dAngular_momentum dt in particle {:?}",particle.tides.parameters.output.dangular_momentum_dt );
            }
        }
        
        if let TidesEffect::KaulaCoplanarOrbitingBody = particle.tides.effect {

            // let orbital_elements = tools::calculate_keplerian_orbital_elements_2(G*(tidal_host_particle.mass+particle.mass), particle.heliocentric_position, particle.heliocentric_velocity);
            // //calculate keplerian elements return (a, q, eccentricity, i, p, n, l, orbital_period)
            // let semi_major_axis:f64                 = orbital_elements.0; //in AU
            // // println!("The SemiMajor Axis {} ", semi_major_axis);
            // // panic!("At the disco");
            // let eccentricity:f64;//                = orbital_elements.2;
            // if orbital_elements.2 < 1e-7 {
            //     eccentricity                        = 0.;
            // } else {
            //     eccentricity                        = orbital_elements.2;
            // }
            // let mut inclination:f64                 = orbital_elements.3;
            // if inclination<1e-7{
            //     inclination                         = 0.;}
            // let mean_anomaly:f64                    = orbital_elements.6;
            // let orbital_period:f64                  = orbital_elements.7;

            let gm = G*(tidal_host_particle.mass+particle.mass);
            let ( semi_major_axis, perihelion_distance, eccentricity, inclination, longitude_perihelion, longitude_of_ascending_node, mean_anomaly, orbital_period 
            ) = tools::calculate_keplerian_orbital_elements_2( gm, particle.heliocentric_position, particle.heliocentric_velocity);

            println!("\n\n from Tide.rs --- orbital elements \n\t a={} q={} e={}  \n\t i={}rad ={}deg \n\t p={}rad ={}deg  \n\t n={}rad ={}deg \n\t l={}rad ={}deg", 
            semi_major_axis, perihelion_distance, eccentricity, inclination, inclination/DEG2RAD, longitude_perihelion, longitude_perihelion/DEG2RAD, longitude_of_ascending_node, longitude_of_ascending_node/DEG2RAD, mean_anomaly, mean_anomaly/DEG2RAD );
            println!("\t \n");
            panic!(" At the Disco ");

            let _argument_perihelion:f64             = longitude_of_ascending_node -longitude_perihelion;
            let orbital_frequency:f64               = TWO_PI / (orbital_period *DAY); // in s^-1

            let eccentricity_function_g_20q  = eccentricty_function_g_20q(eccentricity);
            let _inclination_function        = inclination_function_f2lm(inclination);
            let _associated_legendre_p_2m    = associated_legendre_p_2m(0.); // legendre (cos(theta=\pi/2)=0)) . in coplanar approximation
            let im_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.imaginary_part_love_number;
            let re_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.real_part_love_number;
            let w_lmpq      = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.love_number_excitation_frequency;
            let nm_data:f64 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.num_datapoints;
            let spin:f64 = ( particle.spin.x.powi(2) + particle.spin.y.powi(2) + particle.spin.z.powi(2) ).sqrt() / DAY;
            let _t = current_time*24.;
            let n = orbital_frequency;
            let l = mean_anomaly;
            let mut sum_over_q:f64 = 0.;
            let mut sum_over_j:f64 = 0.;

            particle.tides.parameters.internal.a = semi_major_axis;
            particle.tides.parameters.internal.spin = spin;
            particle.tides.parameters.internal.orbital_frequency = orbital_frequency;

            //=============================================================================== 
            let mut tidal_torque_kaula_coplanar:f64=0.;
            let mut do_calculation:bool=false;
            if (1.- particle.tides.parameters.internal.check_excitative_frequ / sigma_2mpq(2., 0., 0., spin, orbital_frequency)).abs() > 0.001 {
                // //println!("Excitative frequency Variation {}", (1.- particle.tides.parameters.internal.check_excitative_frequ / sigma_2mpq(2., 0., 0., spin, orbital_frequency)).abs() );
                particle.tides.parameters.internal.check_excitative_frequ = sigma_2mpq(2., 0., 0., spin, orbital_frequency);
                do_calculation = true;
            }
            do_calculation = true;
            if do_calculation {
                // if (1. -spin/orbital_frequency).abs() < 0.00001 {
                //     spin = orbital_frequency;
                // }//===========================================================================================================================================================================
                
                //======================================================================================================== Q = -2
                // println!("\nsigma220_2");
                particle.tides.parameters.internal.excitative_frequ_sigma220q[0] = sigma_2mpq(2., 0., -2., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma220q[0], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma220q[0] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma220q[0] = kaula.0 ;
                //=============go to output
                particle.tides.parameters.internal.re_love_number_sigma220_2 = kaula.0;
                particle.tides.parameters.internal.im_love_number_sigma220_2 = kaula.1;
                particle.tides.parameters.internal.sigma220_2_excitative_frequency = particle.tides.parameters.internal.excitative_frequ_sigma220q[0];
                // println!("TIDES.RS");
                // println!("ImK2 {} ", particle.tides.parameters.internal.im_love_number_sigma220_2, );
                // println!("Freq {}", particle.tides.parameters.internal.sigma220_2_excitative_frequency);
                // println!("\nsigma201_2");
                particle.tides.parameters.internal.excitative_frequ_sigma201q[0] = sigma_2mpq(0., 1., -2., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma201q[0], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma201q[0] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma201q[0] = kaula.0 ;
                //=============go to output
                particle.tides.parameters.internal.re_love_number_sigma201_2 = kaula.0;
                particle.tides.parameters.internal.im_love_number_sigma201_2 = kaula.1;
                particle.tides.parameters.internal.sigma201_2_excitative_frequency = particle.tides.parameters.internal.excitative_frequ_sigma201q[0];


                // //println!("sigma220_2:{} \t Im:{} \t Re:{}", particle.tides.parameters.internal.sigma220_2_excitative_frequency, particle.tides.parameters.internal.im_love_number_sigma220_2, particle.tides.parameters.internal.re_love_number_sigma220_2);
                
                //======================================================================================================== Q = -1
                // println!("sigma220_1");
                particle.tides.parameters.internal.excitative_frequ_sigma220q[1] = sigma_2mpq(2., 0., -1., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma220q[1], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma220q[1] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma220q[1] = kaula.0 ;
                //=============go to output
                particle.tides.parameters.internal.re_love_number_sigma220_1 = kaula.0;
                particle.tides.parameters.internal.im_love_number_sigma220_1 = kaula.1;
                particle.tides.parameters.internal.sigma220_1_excitative_frequency = particle.tides.parameters.internal.excitative_frequ_sigma220q[1];
                // println!("\nsigma201_1");
                particle.tides.parameters.internal.excitative_frequ_sigma201q[1] = sigma_2mpq(0., 1., -1., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma201q[1], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma201q[1] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma201q[1] = kaula.0 ;
                //=============go to output
                particle.tides.parameters.internal.re_love_number_sigma201_1 = kaula.0;
                particle.tides.parameters.internal.im_love_number_sigma201_1 = kaula.1;
                particle.tides.parameters.internal.sigma201_1_excitative_frequency = particle.tides.parameters.internal.excitative_frequ_sigma201q[1];


                //======================================================================================================== Q = 0
                // println!("sigma2200");
                particle.tides.parameters.internal.excitative_frequ_sigma220q[2] = sigma_2mpq(2., 0., 0., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma220q[2], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma220q[2] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma220q[2] = kaula.0 ;
                //=============go to output
                particle.tides.parameters.internal.re_love_number_sigma2200 = kaula.0;
                particle.tides.parameters.internal.im_love_number_sigma2200 = kaula.1;
                particle.tides.parameters.internal.sigma2200_excitative_frequency = particle.tides.parameters.internal.excitative_frequ_sigma220q[2];

                // println!("\nsigma2010");
                particle.tides.parameters.internal.excitative_frequ_sigma201q[2] = sigma_2mpq(0., 1., 0., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma201q[2], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma201q[2] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma201q[2] = kaula.0 ;
                //=============go to output
                particle.tides.parameters.internal.re_love_number_sigma2010 = kaula.0;
                particle.tides.parameters.internal.im_love_number_sigma2010 = kaula.1;
                particle.tides.parameters.internal.sigma2010_excitative_frequency = particle.tides.parameters.internal.excitative_frequ_sigma201q[2];

                //======================================================================================================== Q = 1
                // println!("sigma2201");
                particle.tides.parameters.internal.excitative_frequ_sigma220q[3] = sigma_2mpq(2., 0., 1., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma220q[3], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma220q[3] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma220q[3] = kaula.0 ;
                //=============go to output
                particle.tides.parameters.internal.re_love_number_sigma2201 = kaula.0;
                particle.tides.parameters.internal.im_love_number_sigma2201 = kaula.1;
                particle.tides.parameters.internal.sigma2201_excitative_frequency = particle.tides.parameters.internal.excitative_frequ_sigma220q[3];
                // println!("\nsigma2011");
                particle.tides.parameters.internal.excitative_frequ_sigma201q[3] = sigma_2mpq(0., 1., 1., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma201q[3], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma201q[3] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma201q[3] = kaula.0 ;
                //=============go to output
                particle.tides.parameters.internal.re_love_number_sigma2011 = kaula.0;
                particle.tides.parameters.internal.im_love_number_sigma2011 = kaula.1;
                particle.tides.parameters.internal.sigma2011_excitative_frequency = particle.tides.parameters.internal.excitative_frequ_sigma201q[3];

                //======================================================================================================== Q = 2
                // println!("sigma2202");
                particle.tides.parameters.internal.excitative_frequ_sigma220q[4] = sigma_2mpq(2., 0., 2., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma220q[4], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma220q[4] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma220q[4] = kaula.0 ;
                //=============go to output
                particle.tides.parameters.internal.re_love_number_sigma2202 = kaula.0;
                particle.tides.parameters.internal.im_love_number_sigma2202 = kaula.1;
                particle.tides.parameters.internal.sigma2202_excitative_frequency = particle.tides.parameters.internal.excitative_frequ_sigma220q[4];
                // println!("\nsigma2012");
                particle.tides.parameters.internal.excitative_frequ_sigma201q[4] = sigma_2mpq(0., 1., 2., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma201q[4], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma201q[4] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma201q[4] = kaula.0 ;
                //=============go to output
                particle.tides.parameters.internal.re_love_number_sigma2012 = kaula.0;
                particle.tides.parameters.internal.im_love_number_sigma2012 = kaula.1;
                particle.tides.parameters.internal.sigma2012_excitative_frequency = particle.tides.parameters.internal.excitative_frequ_sigma201q[4];

            }
            //------------------------------------------------------------------------------- Eq Torque
            // let mut sum_g_im_k2_over_q:f64=0.;
            // for x in 0..5{
            //     sum_g_im_k2_over_q = sum_g_im_k2_over_q + ( eccentricity_function[x].powi(2) * particle.tides.parameters.internal.im_love_number_sigma220q[x] ) ;
            // }
            // tidal_torque_kaula_coplanar[1] = - (3./2.)*( ( G* tidal_host_particle.mass.powi(2)* particle.radius.powi(5) ) / semi_major_axis.powi(6) ) *sum_g_im_k2_over_q;
            //------------------------------------------------------------------------------- Torque No simplification Test implementation
            sum_over_q=0.;
            for q in -2..3{
                sum_over_j=0.;
                for j in -2..3{
                    // let phase_term:f64 = ((j as f64)-(q as f64))*n*t;
                    let phase_term:f64 = ((j as f64)-(q as f64))*l;
                    // let phase_term:f64 = 0.;
                    let frequence:f64 = sigma_2mpq(2., 0., (q as f64), spin, orbital_frequency);
                    // sum_over_j = sum_over_j +eccentricity_function_g_20q[(j+2) as usize]*( phase_term.sin()*particle.tides.parameters.internal.re_love_number_sigma220q[(q+2) as usize] +phase_term.cos()*particle.tides.parameters.internal.im_love_number_sigma220q[(q+2) as usize]);
                    sum_over_j = sum_over_j +eccentricity_function_g_20q[(j+2) as usize]*( phase_term.sin()*find_rek2(frequence, nm_data, re_k2, w_lmpq) +phase_term.cos()*find_imk2(frequence, nm_data, im_k2, w_lmpq) ) ;
                }
                sum_over_q = sum_over_q +eccentricity_function_g_20q[(q+2) as usize]*sum_over_j;
            }
            tidal_torque_kaula_coplanar = -(3./2.)*( (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5) ) / semi_major_axis.powi(6) ) *sum_over_q;

            particle.tides.parameters.internal.tidal_torque_z = tidal_torque_kaula_coplanar;
            //===============================================================================

            let torque_due_to_tides_x: f64 = 0.;

            let torque_due_to_tides_y: f64 = 0.;    

            // let torque_due_to_tidez_z: f64 = particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide *particle.heliocentric_distance;
            let torque_due_to_tidez_z: f64 = tidal_torque_kaula_coplanar;

            let factor = -1.0;
            if central_body {
                // Integration of the spin (total torque tides):
                dangular_momentum_dt.x += factor * torque_due_to_tides_x;
                dangular_momentum_dt.y += factor * torque_due_to_tides_y;
                dangular_momentum_dt.z += factor * torque_due_to_tidez_z;
            } else {
                particle.tides.parameters.output.dangular_momentum_dt.x = factor * torque_due_to_tides_x;
                particle.tides.parameters.output.dangular_momentum_dt.y = factor * torque_due_to_tides_y;
                particle.tides.parameters.output.dangular_momentum_dt.z = factor * torque_due_to_tidez_z;
                //println!("dAngular_momentum dt in particle {:?}",particle.tides.parameters.output.dangular_momentum_dt );
            }

        }
        particle.tides.parameters.internal.dangular_momentum_dt_by_tides = particle.tides.parameters.output.dangular_momentum_dt;
            
    }

    // if central_body {
    //     // - Equation 25 from Bolmont et al. 2015
    //     tidal_host_particle.tides.parameters.output.dangular_momentum_dt.x = dangular_momentum_dt.x;
    //     tidal_host_particle.tides.parameters.output.dangular_momentum_dt.y = dangular_momentum_dt.y;
    //     tidal_host_particle.tides.parameters.output.dangular_momentum_dt.z = dangular_momentum_dt.z;
    // }
    // println!("dRotationnal_Angular_momentum dt in central body {:?}",tidal_host_particle.tides.parameters.output.dangular_momentum_dt );
}

pub fn calculate_orthogonal_component_of_the_tidal_force(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>, current_time:f64) {
    let mut tidal_host_particle = tidal_host_particle;
    let mut particles = particles;
    let mut more_particles = more_particles;
    let mut star_planet_dependent_dissipation_factors = star_planet_dependent_dissipation_factors;
    let central_body = true;
    calculate_orthogonal_component_of_the_tidal_force_for(central_body, &mut tidal_host_particle, &mut particles, &mut more_particles, &mut star_planet_dependent_dissipation_factors, current_time);
    calculate_orthogonal_component_of_the_tidal_force_for(!central_body, &mut tidal_host_particle, &mut particles, &mut more_particles, &mut star_planet_dependent_dissipation_factors, current_time);
}

fn calculate_orthogonal_component_of_the_tidal_force_for(central_body:bool, tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>, current_time: f64) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody = particle.tides.effect {

            //// Only calculate tides if planet is not in disk
            //if particle.disk_interaction_time == 0.0 {

            // (distance to star)^7
            let distance_7 = particle.tides.parameters.internal.distance.powi(7);

            //// Tidal force calculation (star) :: Only orthogonal component is needed
            if central_body {
                // - Third line of Equation 5 from Bolmont et al. 2015
                //   This expression has R**10 (instead of R**5 in Eq. 5) 
                //   because it uses sigma (i.e., scaled_dissipation_factor) 
                //   and not k2$\Delta$t (between k2$\Delta$t and sigma 
                //   there is a R**5 factor as shown in Equation 28)
                //   - k2 is love number
                let star_scaled_dissipation_factor = planet_dependent_dissipation_factor(&star_planet_dependent_dissipation_factors, &tidal_host_particle.id, tidal_host_particle.evolution, tidal_host_particle.tides.parameters.internal.scaled_dissipation_factor);
                particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 4.5 * (particle.mass.powi(2))
                                                * (tidal_host_particle.radius.powi(10)) 
                                                * star_scaled_dissipation_factor / distance_7;
            } else {
                // - Second line of Equation 5 from Bolmont et al. 2015
                //   This expression has R**10 (instead of R**5 in Eq. 5) 
                //   because it uses sigma (i.e., scaled_dissipation_factor) 
                //   and not k2$\Delta$t (between k2$\Delta$t and sigma 
                //   there is a R**5 factor as shown in Equation 28)
                //   - k2 is love number
                particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 4.5 * (tidal_host_particle.mass.powi(2))
                                                * (particle.radius.powi(10))
                                                * particle.tides.parameters.internal.scaled_dissipation_factor / distance_7;
                // println!("Ortho force {}", particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide);
                // SBC
                ////////println!("> {:e} {:e} {:e} {:e}", tidal_host_particle.mass_g, particle.radius.powi(10), particle.scaled_dissipation_factor, distance_7);
                ////////println!("> {:e} {:e} {:e}", particle.tides.coordinates.position.x, particle.tides.coordinates.position.y, particle.tides.coordinates.position.z);
            }
            //} else {
                //if central_body {
                    //particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 0.0
                //} else{
                    //particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 0.0
                //}
            //}
        }
        if let TidesEffect::KaulaCoplanarOrbitingBody = particle.tides.effect {
            if central_body {
                particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 0.;
            } else {

                // let orbital_elements = tools::calculate_keplerian_orbital_elements_2(G*(tidal_host_particle.mass+particle.mass), particle.heliocentric_position, particle.heliocentric_velocity);
                // //calculate keplerian elements return (a, q, eccentricity, i, p, n, l, orbital_period)
                // let semi_major_axis:f64                 = orbital_elements.0; //in AU
                // let mut eccentricity:f64                = orbital_elements.2;
                // if eccentricity<1e-7{
                //     eccentricity                        = 0.;}
                // let mut inclination:f64                 = orbital_elements.3;
                // if inclination<1e-7{
                //     inclination                         = 0.;}
                // let longitude_perihelion:f64            = orbital_elements.4;
                // let longitude_of_ascending_node:f64     = orbital_elements.5;
                // let argument_perihelion:f64             = longitude_of_ascending_node -longitude_perihelion;
                // let mean_anomaly:f64                    = orbital_elements.6;
                // let orbital_period:f64                  = orbital_elements.7;

                let gm = G*(tidal_host_particle.mass+particle.mass);
                let ( semi_major_axis, _perihelion_distance, eccentricity, inclination, longitude_perihelion, longitude_of_ascending_node, mean_anomaly, orbital_period 
                ) = tools::calculate_keplerian_orbital_elements_2( gm, particle.heliocentric_position, particle.heliocentric_velocity);

                let argument_perihelion:f64             = longitude_of_ascending_node -longitude_perihelion;
                let orbital_frequency:f64               = TWO_PI / (orbital_period *DAY); // in s^-1

                // let phi:f64 =0.; // atan(particle.spin.y / particle.spin.x ) ; // in rad
                let spin:f64 = ( particle.spin.x.powi(2) + particle.spin.y.powi(2) + particle.spin.z.powi(2) ).sqrt() / DAY; //this is [rad.s^-1]
                let im_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.imaginary_part_love_number;
                let re_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.real_part_love_number;
                let w_lmpq      = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.love_number_excitation_frequency;
                let nm_data:f64 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.num_datapoints;
                
                let eccentricity_function_g_20q = eccentricty_function_g_20q(eccentricity);
                let inclination_function        = inclination_function_f2lm(inclination);
                let associated_legendre_p_2m    = associated_legendre_p_2m(0.); // legendre (cos(theta=\pi/2)=0)) . in coplanar approximation

                let mut sum_over:f64 = 0.;
                let n:f64 = orbital_frequency;
                let time:f64 = current_time*24.;
                let phi:f64 =spin*time;
                let l:f64 = mean_anomaly;
                let mut _integer_q:f64 = -2.;
                let mut _integer_j:f64 = -2.;

                // ---  Is 1_kaula or 2 kaula transformation ? -- //
                let kaula_1:bool=false;
                let kaula_2:bool=true;

                if kaula_1 {
                    // -------------------------------------------------------------------------------------------------------------------------------Eq48 Ortho-Force 1Kaula 2nd implementation
                    let mut sum_over_q:f64 = 0.;

                    sum_over_q=0.;
                    for q in -2..3{
                        let phase_term_220q:f64 = (2.+(q as f64))*l -2.*spin*time +2.*argument_perihelion +2.*longitude_of_ascending_node -2.*phi;
                        // let phase_term_220q:f64 = 0.;
                        let frequence:f64=sigma_2mpq(2., 0., (q as f64), spin, orbital_frequency);
                        // sum_over_q = sum_over_q +eccentricity_function_g_20q[(q+2) as usize] *( phase_term_220q.cos()*find_imk2(frequence, nm_data, im_k2, w_lmpq) -phase_term_220q.sin()*find_rek2(frequence, nm_data, re_k2, w_lmpq)  );
                        sum_over_q = sum_over_q +eccentricity_function_g_20q[(q+2) as usize] *( phase_term_220q.cos()*find_imk2(frequence, nm_data, im_k2, w_lmpq) +phase_term_220q.sin()*find_rek2(frequence, nm_data, re_k2, w_lmpq)  ); // Gab sign
                    }
                    let sum_term:f64 = -(1./2.)*associated_legendre_p_2m[2]*sum_over_q;
                    particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = ( (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(3)*particle.heliocentric_distance.powi(4) )) *sum_term;
                }
                
                if kaula_2 {
                    // -------------------------------------------------------------------------------------------------------------------------------Eq48 Ortho-Force 2Kaula 2nd implementation
                    let mut sum_over_q:f64 = 0.;
                    let mut sum_over_j:f64 = 0.;

                    sum_over_q = 0.;
                    for q in -2..3{
                        sum_over_j = 0.;
                        for j in -2..3{
                            let phase_term:f64 = ((j as f64)-(q as f64))*l;
                            // let phase_term:f64 = 0.;
                            let frequence:f64 = sigma_2mpq(2., 0., (q as f64), spin, orbital_frequency);
                            sum_over_j = sum_over_j +eccentricity_function_g_20q[(j+2) as usize] *( phase_term.sin()*find_rek2(frequence, nm_data, re_k2, w_lmpq) +phase_term.cos()*find_imk2(frequence, nm_data, im_k2, w_lmpq)  );
                        }
                        sum_over_q = sum_over_q +eccentricity_function_g_20q[(q+2) as usize] *sum_over_j;
                    }
                    particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = -(3./2.) *( (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(6)*particle.heliocentric_distance )) *sum_over_q;
                    // particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = particle.tides.parameters.internal.tidal_torque_z / particle.heliocentric_distance ;
                }
                
                // particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 0.0; //ORTHO FORCE NULL
            }
        }
    }
}

pub fn calculate_radial_component_of_the_tidal_force(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>, current_time:f64) {
    let star_mass_2 = tidal_host_particle.mass * tidal_host_particle.mass;

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody = particle.tides.effect {
            let planet_mass_2 = particle.mass * particle.mass;
            // Conservative part of the radial tidal force
            let radial_component_of_the_tidal_force_conservative_part = -3.0 * K2 / particle.tides.parameters.internal.distance.powi(7)
                        * (planet_mass_2 * tidal_host_particle.radius.powi(5) * tidal_host_particle.tides.parameters.input.love_number 
                        + star_mass_2 * particle.radius.powi(5) * particle.tides.parameters.input.love_number);

            // Dissipative part of the radial tidal force:
            let factor1 = -13.5 * particle.tides.parameters.internal.radial_velocity / particle.tides.parameters.internal.distance.powi(8);
            let star_scaled_dissipation_factor = planet_dependent_dissipation_factor(&star_planet_dependent_dissipation_factors, &particle.id, tidal_host_particle.evolution, tidal_host_particle.tides.parameters.internal.scaled_dissipation_factor);
            let term1 = planet_mass_2
                        * tidal_host_particle.radius.powi(10)
                        * star_scaled_dissipation_factor;
            let term2 = star_mass_2
                        * particle.radius.powi(10)
                        * particle.tides.parameters.internal.scaled_dissipation_factor;
            // If we consider the star as a point mass (used for denergy_dt calculation):
            particle.tides.parameters.internal.radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass = factor1 * term2;
            let radial_component_of_the_tidal_force_dissipative_part = particle.tides.parameters.internal.radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass + factor1 * term1;

            // Sum of the dissipative and conservative part of the radial force
            // - First line Equation 5 from Bolmont et al. 2015
            particle.tides.parameters.internal.radial_component_of_the_tidal_force = radial_component_of_the_tidal_force_conservative_part + radial_component_of_the_tidal_force_dissipative_part;
            // println!("Radial force {}", particle.tides.parameters.internal.radial_component_of_the_tidal_force);
        }

        if let TidesEffect::KaulaCoplanarOrbitingBody = particle.tides.effect {

            // let orbital_elements = tools::calculate_keplerian_orbital_elements_2(G*(tidal_host_particle.mass+particle.mass), particle.heliocentric_position, particle.heliocentric_velocity);
            // //calculate keplerian el return (a, q, eccentricity, i, p, n, l, orbital_period) 
            // let semi_major_axis:f64                 = orbital_elements.0; //in AU
            // let mut eccentricity:f64                = orbital_elements.2;
            // if eccentricity<1e-7
            // {
            //     eccentricity                        = 0.;   
            // }
            // let mut inclination:f64                 = orbital_elements.3;
            // if inclination<1e-7{
            //     inclination                         = 0.;
            // }
            // let longitude_perihelion:f64            = orbital_elements.4;
            // let longitude_of_ascending_node:f64     = orbital_elements.5;
            // let argument_perihelion:f64             = longitude_of_ascending_node -longitude_perihelion;
            // let mean_anomaly:f64                    = orbital_elements.6;
            // let orbital_period:f64                  = orbital_elements.7;

            let gm = G*(tidal_host_particle.mass+particle.mass);
            let ( semi_major_axis, _perihelion_distance, eccentricity, inclination, longitude_perihelion, longitude_of_ascending_node, mean_anomaly, orbital_period 
            ) = tools::calculate_keplerian_orbital_elements_2( gm, particle.heliocentric_position, particle.heliocentric_velocity);

            let argument_perihelion:f64             = longitude_of_ascending_node -longitude_perihelion;
            let orbital_frequency:f64               = TWO_PI / (orbital_period *DAY); // in s^-1
            let spin:f64                            = ( particle.spin.x.powi(2) + particle.spin.y.powi(2) + particle.spin.z.powi(2) ).sqrt() / DAY; //this is [rad.s^-1]

            let eccentricity_function_g_20q = eccentricty_function_g_20q(eccentricity);
            let eccentricity_function_g_21q = eccentricty_function_g_21q(eccentricity);
            let inclination_function        = inclination_function_f2lm(inclination);
            let associated_legendre_p_2m    = associated_legendre_p_2m(0.); //0. in coplanar approximation

            let im_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.imaginary_part_love_number;
            let re_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.real_part_love_number;
            let w_lmpq      = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.love_number_excitation_frequency;
            let nm_data:f64 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.num_datapoints;

            // let mut _excitative_frequency: f64;
            let mut sum_over:f64 = 0.;
            let l:f64 = mean_anomaly;
            let n:f64 = orbital_frequency;
            let time:f64 = current_time*24.; // in seconds
            let phi:f64 =spin*time;
            let mut _integer_q:f64 = -2.;
            let mut _integer_j:f64 = -2.;

            // ---  Is 1_kaula or 2 kaula transformation ? -- //
            let kaula_1:bool=false;
            let kaula_2:bool=true;

            if kaula_1 { 
                // ------------------------------------------------------------------------------------------------------------------------------------Radial-Force 1Kaula 2nd implementation
                let mut sum_over_q:f64=0.;
                sum_over_q=0.;
                for q in -2..3{
                    let phase_term_201q:f64 = (q as f64)*l;
                    // let phase_term_201q:f64 = 0.;
                    let frequence:f64=sigma_2mpq(0., 1., (q as f64), spin, orbital_frequency);
                    sum_over_q = sum_over_q +eccentricity_function_g_21q[(q+2) as usize]*( phase_term_201q.cos()*find_rek2(frequence, nm_data, re_k2, w_lmpq) +phase_term_201q.sin()*find_imk2(frequence, nm_data, im_k2, w_lmpq)  );
                }
                let term1:f64 = (3./2.)*associated_legendre_p_2m[0]*sum_over_q;
                
                sum_over_q=0.;
                for q in -2..3{
                    let phase_term_220q:f64 = (2.+(q as f64))*l -2.*spin*time +2.*argument_perihelion +2.*longitude_of_ascending_node -2.*phi;
                    // let phase_term_220q:f64 = 0.;
                    let frequence:f64=sigma_2mpq(2., 0., (q as f64), spin, orbital_frequency);
                    sum_over_q = sum_over_q +eccentricity_function_g_20q[(q+2) as usize]*( phase_term_220q.cos()*find_rek2(frequence, nm_data, re_k2, w_lmpq) +phase_term_220q.sin()*find_imk2(frequence, nm_data, im_k2, w_lmpq)  );
                }
                let term2:f64 = -(3./4.)*associated_legendre_p_2m[2]*sum_over_q;
                particle.tides.parameters.internal.radial_component_of_the_tidal_force = ( (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(3) *particle.heliocentric_distance.powi(4)) ) *(term1+term2);
            }

            if kaula_2 {
                // ------------------------------------------------------------------------------------------------------------------------------------Eq48 Radial-Force 2Kaula 2nd implementation
                let mut sum_over_q:f64 = 0.;
                let mut sum_over_j:f64 = 0.;
                
                sum_over_q = 0.;
                for q in -2..3{
                    sum_over_j = 0.;
                    for j in -2..3{
                        let phase_term:f64 = ((j as f64)-(q as f64))*l;
                        // let phase_term:f64 = 0.;
                        let frequence:f64 = sigma_2mpq(0., 1., (q as f64), spin, orbital_frequency);
                        sum_over_j = sum_over_j +eccentricity_function_g_21q[(j+2) as usize]*( phase_term.cos()*find_rek2(frequence, nm_data, re_k2, w_lmpq) -phase_term.sin()*find_imk2(frequence, nm_data, im_k2, w_lmpq)  );
                    }
                    sum_over_q = sum_over_q +(3./4.)*eccentricity_function_g_21q[(q+2) as usize]*sum_over_j;
                }
                let term1:f64 = sum_over_q;
                
                sum_over_q = 0.;
                for q in -2..3{
                    sum_over_j = 0.;
                    for j in -2..3{
                        let phase_term:f64 = ((j as f64)-(q as f64))*l;
                        // let phase_term:f64 = 0.;
                        let frequence:f64 = sigma_2mpq(2., 0., (q as f64), spin, orbital_frequency);
                        sum_over_j = sum_over_j +eccentricity_function_g_20q[(j+2) as usize]*( phase_term.cos()*find_rek2(frequence, nm_data, re_k2, w_lmpq) -phase_term.sin()*find_imk2(frequence, nm_data, im_k2, w_lmpq)  );
                    }
                    sum_over_q = sum_over_q +(9./4.)*eccentricity_function_g_20q[(q+2) as usize]*sum_over_j;
                }
                let term2:f64 = sum_over_q;
                particle.tides.parameters.internal.radial_component_of_the_tidal_force = -( (G *tidal_host_particle.mass.powi(2) *particle.radius.powi(5)) / (semi_major_axis.powi(6) *particle.heliocentric_distance) ) *(term1+term2);

            }
            

            // particle.tides.parameters.internal.radial_component_of_the_tidal_force = 0.0; //RADIAL FORCE NULL
        }
    }
}

pub fn calculate_denergy_dt(particles: &mut [Particle], more_particles: &mut [Particle]) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody = particle.tides.effect {
            // - Equation 32 from Bolmont et al. 2015
            //// Instantaneous energy loss dE/dt due to tides
            //// in Msun.AU^2.day^(-3)
            //radial_tidal_force_for_energy_loss_calculation = factor1 * term2; // Ftidr_diss
            let factor2 = particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.tides.parameters.internal.distance;
            particle.tides.parameters.internal.denergy_dt = -((1.0 / particle.tides.parameters.internal.distance * (particle.tides.parameters.internal.radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass + factor2 * particle.tides.parameters.internal.radial_velocity))
                        * (particle.tides.coordinates.position.x*particle.tides.coordinates.velocity.x + particle.tides.coordinates.position.y*particle.tides.coordinates.velocity.y + particle.tides.coordinates.position.z*particle.tides.coordinates.velocity.z)
                        + factor2 
                        * ((particle.spin.y*particle.tides.coordinates.position.z - particle.spin.z*particle.tides.coordinates.position.y - particle.tides.coordinates.velocity.x) * particle.tides.coordinates.velocity.x
                        + (particle.spin.z*particle.tides.coordinates.position.x - particle.spin.x*particle.tides.coordinates.position.z - particle.tides.coordinates.velocity.y) * particle.tides.coordinates.velocity.y
                        + (particle.spin.x*particle.tides.coordinates.position.y - particle.spin.y*particle.tides.coordinates.position.x - particle.tides.coordinates.velocity.z) * particle.tides.coordinates.velocity.z))
                        - (particle.tides.parameters.output.dangular_momentum_dt.x*particle.spin.x + particle.tides.parameters.output.dangular_momentum_dt.y*particle.spin.y + particle.tides.parameters.output.dangular_momentum_dt.z*particle.spin.z);
        }
        if let TidesEffect::KaulaCoplanarOrbitingBody = particle.tides.effect {
                // - Equation 32 from Bolmont et al. 2015
                //// Instantaneous energy loss dE/dt due to tides
                //// in Msun.AU^2.day^(-3)
                //radial_tidal_force_for_energy_loss_calculation = factor1 * term2; // Ftidr_diss
                let factor2 = particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.tides.parameters.internal.distance;
                particle.tides.parameters.internal.denergy_dt = -((1.0 / particle.tides.parameters.internal.distance * (particle.tides.parameters.internal.radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass + factor2 * particle.tides.parameters.internal.radial_velocity))
                            * (particle.tides.coordinates.position.x*particle.tides.coordinates.velocity.x + particle.tides.coordinates.position.y*particle.tides.coordinates.velocity.y + particle.tides.coordinates.position.z*particle.tides.coordinates.velocity.z)
                            + factor2 
                            * ((particle.spin.y*particle.tides.coordinates.position.z - particle.spin.z*particle.tides.coordinates.position.y - particle.tides.coordinates.velocity.x) * particle.tides.coordinates.velocity.x
                            + (particle.spin.z*particle.tides.coordinates.position.x - particle.spin.x*particle.tides.coordinates.position.z - particle.tides.coordinates.velocity.y) * particle.tides.coordinates.velocity.y
                            + (particle.spin.x*particle.tides.coordinates.position.y - particle.spin.y*particle.tides.coordinates.position.x - particle.tides.coordinates.velocity.z) * particle.tides.coordinates.velocity.z))
                            - (particle.tides.parameters.output.dangular_momentum_dt.x*particle.spin.x + particle.tides.parameters.output.dangular_momentum_dt.y*particle.spin.y + particle.tides.parameters.output.dangular_momentum_dt.z*particle.spin.z);
            }
    }
}


pub fn calculate_tidal_acceleration(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], current_time:f64) {
    let factor2 = 1. / tidal_host_particle.mass;
    let mut sum_total_tidal_force = Axes{x:0., y:0., z:0.};

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody = particle.tides.effect {
            let factor1 = 1. / particle.mass;

            // - Equation 6 from Bolmont et al. 2015
            let factor3 = particle.tides.parameters.internal.radial_component_of_the_tidal_force
                            + (particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide + particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide) * particle.tides.parameters.internal.radial_velocity / particle.tides.parameters.internal.distance;
            let total_tidal_force_x = factor3 * particle.tides.coordinates.position.x / particle.tides.parameters.internal.distance
                                    + particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide / particle.tides.parameters.internal.distance
                                        * (tidal_host_particle.spin.y * particle.tides.coordinates.position.z  - tidal_host_particle.spin.z * particle.tides.coordinates.position.y - particle.tides.coordinates.velocity.x)

                                    + particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.tides.parameters.internal.distance 
                                        * (particle.spin.y * particle.tides.coordinates.position.z  - particle.spin.z * particle.tides.coordinates.position.y - particle.tides.coordinates.velocity.x);

            let total_tidal_force_y = factor3 * particle.tides.coordinates.position.y / particle.tides.parameters.internal.distance
                                    + particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide / particle.tides.parameters.internal.distance 
                                        * (tidal_host_particle.spin.z * particle.tides.coordinates.position.x  - tidal_host_particle.spin.x * particle.tides.coordinates.position.z - particle.tides.coordinates.velocity.y)
                                    + particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.tides.parameters.internal.distance 
                                        * (particle.spin.z * particle.tides.coordinates.position.x  - particle.spin.x * particle.tides.coordinates.position.z - particle.tides.coordinates.velocity.y);
            let total_tidal_force_z = factor3 * particle.tides.coordinates.position.z / particle.tides.parameters.internal.distance
                                    + particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide / particle.tides.parameters.internal.distance 
                                        * (tidal_host_particle.spin.x * particle.tides.coordinates.position.y  - tidal_host_particle.spin.y * particle.tides.coordinates.position.x - particle.tides.coordinates.velocity.z)
                                    + particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide / particle.tides.parameters.internal.distance 
                                        * (particle.spin.x * particle.tides.coordinates.position.y  - particle.spin.y * particle.tides.coordinates.position.x - particle.tides.coordinates.velocity.z);
            ////////println!("factor3 {:e} {:e} {:e}", particle.tides.parameters.internal.radial_component_of_the_tidal_force, particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide, particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide);
            ////////println!("d {:e} vrad {:e}", particle.tides.parameters.internal., particle.tides.parameters.internal.radial_velocity);
            ////////println!("total {:e} {:e} {:e}", total_tidal_force_x, total_tidal_force_y, total_tidal_force_z);
            // //println!("The componant of force radial [{}] ortho [{}]", particle.tides.parameters.internal.radial_component_of_the_tidal_force, particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide );
            
            sum_total_tidal_force.x += total_tidal_force_x;
            sum_total_tidal_force.y += total_tidal_force_y;
            sum_total_tidal_force.z += total_tidal_force_z;

            // - Equation 19 from Bolmont et al. 2015 (first term)
            particle.tides.parameters.output.acceleration.x = factor1 * total_tidal_force_x; 
            particle.tides.parameters.output.acceleration.y = factor1 * total_tidal_force_y;
            particle.tides.parameters.output.acceleration.z = factor1 * total_tidal_force_z;
            //println!("\nIn CTL\nThe total tidal force = \tX:{} \tY{} \tZ{}",particle.tides.parameters.output.acceleration.x, particle.tides.parameters.output.acceleration.z, particle.tides.parameters.output.acceleration.z );
        }
        if let TidesEffect::KaulaCoplanarOrbitingBody = particle.tides.effect {
            let factor1 = 1. / particle.mass;
            let factor4 = (particle.mass+tidal_host_particle.mass)/(particle.mass*tidal_host_particle.mass);

            let mut radial_component_of_the_tidal_force = particle.tides.parameters.internal.radial_component_of_the_tidal_force;
            let mut orthogonal_component_of_the_tidal_force = particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide;
            let local_calculation = false;
            if local_calculation {
                // ---- Implement directly into calculate tidal acceleration ---- //
                // let mut radial_component_of_the_tidal_force = 0.;//particle.tides.parameters.internal.radial_component_of_the_tidal_force;
                // let mut orthogonal_component_of_the_tidal_force = 0.;//particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide;

                // let orbital_elements = tools::calculate_keplerian_orbital_elements_2(G*(tidal_host_particle.mass+particle.mass), particle.heliocentric_position, particle.heliocentric_velocity);
                // let semi_major_axis:f64                 = orbital_elements.0; //in AU
                // let mut eccentricity:f64                = orbital_elements.2;
                // if eccentricity<1e-7
                // {
                //     eccentricity                        = 0.;   
                // }
                // let mut inclination:f64                 = orbital_elements.3;
                // if inclination<1e-7{
                //     inclination                         = 0.;
                // }
                // let longitude_perihelion:f64            = orbital_elements.4;
                // let longitude_of_ascending_node:f64     = orbital_elements.5;
                // let argument_perihelion:f64             = longitude_of_ascending_node -longitude_perihelion;
                // let mean_anomaly:f64                    = orbital_elements.6;
                // let orbital_period:f64                  = orbital_elements.7;

                let gm = G*(tidal_host_particle.mass+particle.mass);
                let ( semi_major_axis, _perihelion_distance, eccentricity, inclination, longitude_perihelion, longitude_of_ascending_node, mean_anomaly, orbital_period 
                ) = tools::calculate_keplerian_orbital_elements_2( gm, particle.heliocentric_position, particle.heliocentric_velocity);

                let argument_perihelion:f64             = longitude_of_ascending_node -longitude_perihelion;
                let orbital_frequency:f64               = TWO_PI / (orbital_period *DAY); // in s^-1
                let spin:f64    = ( particle.spin.x.powi(2) + particle.spin.y.powi(2) + particle.spin.z.powi(2) ).sqrt() / DAY; //this is [rad.s^-1]
                // let phi:f64     = 0.; // atan(particle.spin.y / particle.spin.x ) ; // in rad
                let eccentricity_function_g_20q = eccentricty_function_g_20q(eccentricity);
                let eccentricity_function_g_21q = eccentricty_function_g_21q(eccentricity);
                let associated_legendre_p_2m    = associated_legendre_p_2m(0.); //0. in coplanar approximation
                let im_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.imaginary_part_love_number;
                let re_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.real_part_love_number;
                let w_lmpq      = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.love_number_excitation_frequency;
                let nm_data:f64 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.num_datapoints;
                let time:f64    = current_time*24.; // in seconds
                let phi:f64     = spin*time;
                let l:f64 = mean_anomaly;
                // ---------------------------------------- //
                // ---- Calculation of the tidal force ---- //
                let kaula_2 = false;
                if kaula_2 {
                    // panic!("At the Disco");
                    let mut sum_over_q:f64 = 0.;
                    let mut sum_over_j:f64 = 0.;
                    // ---- orthogonal part of tidal force ---- //
                    sum_over_q = 0.;
                    for q in -2..3{
                        sum_over_j = 0.;
                        for j in -2..3{
                            let phase_term:f64 = ((j as f64)-(q as f64))*l;
                            // let phase_term:f64 = 0.;
                            let frequence:f64 = sigma_2mpq(2., 0., (q as f64), spin, orbital_frequency);
                            sum_over_j = sum_over_j +eccentricity_function_g_20q[(j+2) as usize] *( phase_term.sin()*find_rek2(frequence, nm_data, re_k2, w_lmpq) +phase_term.cos()*find_imk2(frequence, nm_data, im_k2, w_lmpq)  );
                        }
                        sum_over_q = sum_over_q +eccentricity_function_g_20q[(q+2) as usize] *sum_over_j;
                    }
                    orthogonal_component_of_the_tidal_force = -(3./2.) *( (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(6))*particle.heliocentric_distance ) *sum_over_q;
                    particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = orthogonal_component_of_the_tidal_force;
                    
                    // ---- radial part of tidal force ---- //
                    sum_over_j = 0.;
                    sum_over_q = 0.;
                    for q in -2..3{
                        sum_over_j = 0.;
                        for j in -2..3{
                            let phase_term:f64 = ((j as f64)-(q as f64))*l;
                            // let phase_term:f64 = 0.;
                            let frequence:f64 = sigma_2mpq(0., 1., (q as f64), spin, orbital_frequency);
                            sum_over_j = sum_over_j +eccentricity_function_g_21q[(j+2) as usize]*( phase_term.cos()*find_rek2(frequence, nm_data, re_k2, w_lmpq) -phase_term.sin()*find_imk2(frequence, nm_data, im_k2, w_lmpq)  );
                        }
                        sum_over_q = sum_over_q +(3./4.)*eccentricity_function_g_21q[(q+2) as usize]*sum_over_j;
                    }
                    let term1:f64 = sum_over_q;
                    
                    sum_over_q = 0.;
                    for q in -2..3{
                        sum_over_j = 0.;
                        for j in -2..3{
                            let phase_term:f64 = ((j as f64)-(q as f64))*l;
                            // let phase_term:f64 = 0.;
                            let frequence:f64 = sigma_2mpq(2., 0., (q as f64), spin, orbital_frequency);
                            sum_over_j = sum_over_j +eccentricity_function_g_20q[(j+2) as usize]*( phase_term.cos()*find_rek2(frequence, nm_data, re_k2, w_lmpq) -phase_term.sin()*find_imk2(frequence, nm_data, im_k2, w_lmpq)  );
                        }
                        sum_over_q = sum_over_q +(9./4.)*eccentricity_function_g_20q[(q+2) as usize]*sum_over_j;
                    }
                    let term2:f64 = sum_over_q;

                    radial_component_of_the_tidal_force = -( (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(6) * particle.heliocentric_distance) ) *(term1+term2);
                    particle.tides.parameters.internal.radial_component_of_the_tidal_force = radial_component_of_the_tidal_force;
                }
                // --------------------------------------------------------
                let kaula_1 = false;
                if kaula_1 {
                    // panic!("At the disco");
                    let mut sum_over_q:f64=0.;

                    sum_over_q=0.;
                    for q in -2..3{
                        let phase_term_220q:f64 = (2.+(q as f64))*l -2.*spin*time +2.*argument_perihelion +2.*longitude_of_ascending_node -2.*phi;
                        // let phase_term_220q:f64 = 0.;
                        let frequence:f64=sigma_2mpq(2., 0., (q as f64), spin, orbital_frequency);
                        sum_over_q = sum_over_q +eccentricity_function_g_20q[(q+2) as usize] *( phase_term_220q.cos()*find_imk2(frequence, nm_data, im_k2, w_lmpq) -phase_term_220q.sin()*find_rek2(frequence, nm_data, re_k2, w_lmpq)  );
                        // sum_over_q = sum_over_q +eccentricity_function_g_20q[(q+2) as usize] *( phase_term_220q.cos()*find_imk2(frequence, nm_data, im_k2, w_lmpq) +phase_term_220q.sin()*find_rek2(frequence, nm_data, re_k2, w_lmpq)  ); // Gabriel sign
                    }
                    let sum_term:f64 = -(1./2.)*associated_legendre_p_2m[2]*sum_over_q;

                    orthogonal_component_of_the_tidal_force = ( (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(3))*particle.heliocentric_distance.powi(4) ) *sum_term;
                    particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = orthogonal_component_of_the_tidal_force;
                    // --------------------------------------------------------
                    sum_over_q=0.;
                    for q in -2..3{
                        let phase_term_201q:f64 = (q as f64)*l;
                        // let phase_term_201q:f64 = 0.;
                        let frequence:f64=sigma_2mpq(0., 1., (q as f64), spin, orbital_frequency);
                        sum_over_q = sum_over_q +eccentricity_function_g_21q[(q+2) as usize]*( phase_term_201q.cos()*find_rek2(frequence, nm_data, re_k2, w_lmpq) +phase_term_201q.sin()*find_imk2(frequence, nm_data, im_k2, w_lmpq)  );
                    }
                    let term1:f64 = (3./2.)*associated_legendre_p_2m[0]*sum_over_q;
                    
                    sum_over_q=0.;
                    for q in -2..3{
                        let phase_term_220q:f64 = (2.+(q as f64))*l -2.*spin*time +2.*argument_perihelion +2.*longitude_of_ascending_node -2.*phi;
                        // let phase_term_220q:f64 = 0.;
                        let frequence:f64=sigma_2mpq(2., 0., (q as f64), spin, orbital_frequency);
                        sum_over_q = sum_over_q +eccentricity_function_g_20q[(q+2) as usize]*( phase_term_220q.cos()*find_rek2(frequence, nm_data, re_k2, w_lmpq) +phase_term_220q.sin()*find_imk2(frequence, nm_data, im_k2, w_lmpq)  );
                    }
                    let term2:f64 = -(3./4.)*associated_legendre_p_2m[2]*sum_over_q;
                    radial_component_of_the_tidal_force = ( (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(3) * particle.heliocentric_distance.powi(4)) ) *(term1+term2);
                    particle.tides.parameters.internal.radial_component_of_the_tidal_force = radial_component_of_the_tidal_force;
                }
                // --------------------------------------------------------

            } // 
            
            let distance:f64 = ( particle.tides.coordinates.position.x.powi(2) +particle.tides.coordinates.position.y.powi(2) +particle.tides.coordinates.position.z.powi(2) ).sqrt();
            let total_tidal_force_x = radial_component_of_the_tidal_force *(particle.tides.coordinates.position.x / distance) - orthogonal_component_of_the_tidal_force *(particle.tides.coordinates.position.y / distance) ;
            let total_tidal_force_y = radial_component_of_the_tidal_force *(particle.tides.coordinates.position.y / distance) + orthogonal_component_of_the_tidal_force *(particle.tides.coordinates.position.x / distance) ;
            let total_tidal_force_z = 0.0; 
            // let total_tidal_force_x = radial_component_of_the_tidal_force *(particle.tides.coordinates.position.x / particle.tides.parameters.internal.distance) - orthogonal_component_of_the_tidal_force *(particle.tides.coordinates.position.y / particle.tides.parameters.internal.distance) ;
            // let total_tidal_force_y = radial_component_of_the_tidal_force *(particle.tides.coordinates.position.y / particle.tides.parameters.internal.distance) + orthogonal_component_of_the_tidal_force *(particle.tides.coordinates.position.x / particle.tides.parameters.internal.distance) ;
            // let total_tidal_force_z = 0.0; 

            sum_total_tidal_force.x += total_tidal_force_x;
            sum_total_tidal_force.y += total_tidal_force_y;
            sum_total_tidal_force.z += total_tidal_force_z;

            // - Equation 19 from Bolmont et al. 2015 (first term)
            particle.tides.parameters.output.acceleration.x = factor1 * total_tidal_force_x; 
            particle.tides.parameters.output.acceleration.y = factor1 * total_tidal_force_y;
            particle.tides.parameters.output.acceleration.z = factor1 * total_tidal_force_z;
            //println!("\nIn Kaula\nThe total tidal force = \tX:{} \tY{} \tZ{}",particle.tides.parameters.output.acceleration.x, particle.tides.parameters.output.acceleration.z, particle.tides.parameters.output.acceleration.z );
            
        }
        particle.tides.parameters.internal.force_by_tides = particle.tides.parameters.output.acceleration;
    }

    // - Equation 19 from Bolmont et al. 2015 (second term)
    //for particle in particles.iter_mut() {
        //particle.tides.parameters.output.acceleration.x += factor2 * sum_total_tidal_force.x;
        //particle.tides.parameters.output.acceleration.y += factor2 * sum_total_tidal_force.y;
        //particle.tides.parameters.output.acceleration.z += factor2 * sum_total_tidal_force.z;
    //}
    // Instead of the previous code, keep star tidal acceleration separated:
    tidal_host_particle.tides.parameters.output.acceleration.x = -1.0 * factor2 * sum_total_tidal_force.x;
    tidal_host_particle.tides.parameters.output.acceleration.y = -1.0 * factor2 * sum_total_tidal_force.y;
    tidal_host_particle.tides.parameters.output.acceleration.z = -1.0 * factor2 * sum_total_tidal_force.z;
    // println!("The sum_total_force {:?}", tidal_host_particle.tides.parameters.output.acceleration);

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Let's star the modifications:


pub fn eccentricty_function_g_20q(eccentricity: f64) -> [f64; 5]{

    let mut eccentricity_function = [0_f64;5];
    // let e_1 = eccentricity ;
    // let e_2 = e_1 ;
    // let e_3 = e_2 * e_1 ;
    // let e_4 = e_3 *e_1 ;

    eccentricity_function[0] = 0.;
    eccentricity_function[1] = -(1./2.)*eccentricity +(1./16.)*eccentricity.powi(3) - (15./384.)*eccentricity.powi(5) - (143./18432.)*eccentricity.powi(7);
    eccentricity_function[2] = 1. -(5./2.)*eccentricity.powi(2) +(13./16.)*eccentricity.powi(4) - (35./288.) * eccentricity.powi(6);
    eccentricity_function[3] = (7./2.)*eccentricity -(123./16.)*eccentricity.powi(3) + (489./128.) * eccentricity.powi(5) - (1763./2048.)*eccentricity.powi(7);
    eccentricity_function[4] = (17./2.)*eccentricity.powi(2) -(115./6.)*eccentricity.powi(4) + (601./48.)*eccentricity.powi(6);

    return eccentricity_function;
}
pub fn eccentricty_function_g_21q(eccentricity: f64) -> [f64; 5]{

    let mut eccentricity_function = [0_f64;5];
    // let e_1 = eccentricity ;
    // let e_2 = e_1 ;
    // let e_3 = e_2 * e_1 ;
    // let e_4 = e_3 *e_1 ;

    eccentricity_function[0] = (9./4.)*eccentricity.powi(2) +(7./4.)*eccentricity.powi(4);
    eccentricity_function[1] = (3./2.)*eccentricity +(27./16.)*eccentricity.powi(3);
    eccentricity_function[2] = (1. -eccentricity.powi(2)).powf(-(3./2.));
    eccentricity_function[3] = (3./2.)*eccentricity +(27./16.)*eccentricity.powi(3);
    eccentricity_function[4] = (9./4.)*eccentricity.powi(2) +(7./4.)*eccentricity.powi(4);

    return eccentricity_function;
}

pub fn associated_legendre_p_2m(x: f64) -> [f64;3]{

    let mut p_2m = [0.0;3];
    let x_2 = x*x ;

    p_2m[0] = 0.5*(3.*x_2 -1.);
    p_2m[1] = -3.*((1.-x_2).sqrt());
    p_2m[2] = 3.*(1.-x_2);

    return p_2m;
}

pub fn _eccentricty_function_g(eccentricity: f64) -> [[f64; 15]; 3]{

    let mut eccentricity_function = [[0.0;15];3];
    // in this notation:
    // [0] is -7
    // [1] is -6
    // [2] is -5
    // [3] is -4
    // [4] is -3
    // [5] is -2
    // [6] is -1
    // [7] is 0
    // [8] is 1
    // [9] is 2
    // [10] is 3   
    // [11] is 4
    // [12] is 5
    // [13] is 6
    // [14] is 7

    let ecc:f64=eccentricity;
    let ecc_2:f64=ecc*ecc;
    let ecc_3:f64=ecc_2*ecc;
    let ecc_4:f64=ecc_3*ecc;
    let ecc_5:f64=ecc_4*ecc;
    let ecc_6:f64=ecc_5*ecc;
    let ecc_7:f64=ecc_6*ecc;

    eccentricity_function[0][0] = (15625./129024.)*ecc_7 ;
    eccentricity_function[0][1] = (4./45.)*ecc_6 ;
    eccentricity_function[0][2] = (81./1280.)*ecc_5     +(81./2048.)*ecc_7 ;
    eccentricity_function[0][3] = (1./24.)*ecc_4        +(7./240.)*ecc_6 ;
    eccentricity_function[0][4] = (1./48.)*ecc_3        +(11./768.)*ecc_5       +(313./30720.)*ecc_7 ;
    eccentricity_function[0][5] = 0.;
    eccentricity_function[0][6] = -0.5*ecc              +(1./16.)*ecc_3         -(5./384.)*ecc_5        -(143./18432.)*ecc_7 ;
    eccentricity_function[0][7] = 1.                    -(5./2.)*ecc_2          +(13./16.)*ecc_4        -(35./288.)*ecc_6 ;
    eccentricity_function[0][8] = (7./2.)*ecc           -(123./16.)*ecc_3       +(489./128.)*ecc_5      -(1763./2048.)*ecc_7 ;
    eccentricity_function[0][9] = (17./2.)*ecc_2        -(115./16.)*ecc_4       +(601./48.)*ecc_6 ;
    eccentricity_function[0][10] = (845./48.)*ecc_3     -(32525./768.)*ecc_5    +(208225./6144.)*ecc_7 ;
    eccentricity_function[0][11] = (533./16.)*ecc_4     -(13827./160.)*ecc_6 ;
    eccentricity_function[0][12] = (228347./3840.)*ecc_5 -(3071075./18432.)*ecc_7 ;
    eccentricity_function[0][13] = (73369./720.)*ecc_6 ;
    eccentricity_function[0][14] = (12144273./71680.)*ecc_7 ;

    eccentricity_function[1][0] = (432091./30720.)*ecc_7 ;
    eccentricity_function[1][1] = (3167./320.)*ecc_6 ;
    eccentricity_function[1][2] = (1773./256.)*ecc_5    +(4987./6144.)*ecc_7 ;
    eccentricity_function[1][3] = (77./16.)*ecc_4       +(129./160.)*ecc_6 ;
    eccentricity_function[1][4] = (53./16.)*ecc_3       +(393./256.)*ecc_5      +(24753./10240.)*ecc_7;
    eccentricity_function[1][5] = (9./4.)*ecc_2         +(7./4.)*ecc_4          +(141./64.)*ecc_6 ;
    eccentricity_function[1][6] = (3./2.)*ecc           +(27./16.)*ecc_3        +(261./128.)*ecc_5      +(14309./6144.)*ecc_7 ;
    eccentricity_function[1][7] = (1.-ecc_2).powf(-3./2.) ;
    eccentricity_function[1][8] = eccentricity_function[1][6] ;
    eccentricity_function[1][9] = eccentricity_function[1][5] ;
    eccentricity_function[1][10] = eccentricity_function[1][4] ;
    eccentricity_function[1][11] = eccentricity_function[1][3] ;
    eccentricity_function[1][12] = eccentricity_function[1][2] ;
    eccentricity_function[1][13] = eccentricity_function[1][1] ;
    eccentricity_function[1][14] = eccentricity_function[1][0] ;

    eccentricity_function[2][0] = eccentricity_function[0][14];
    eccentricity_function[2][1] = eccentricity_function[0][13];
    eccentricity_function[2][2] = eccentricity_function[0][12];
    eccentricity_function[2][3] = eccentricity_function[0][11];
    eccentricity_function[2][4] = eccentricity_function[0][10];
    eccentricity_function[2][5] = eccentricity_function[0][9] ;
    eccentricity_function[2][6] = eccentricity_function[0][8] ;
    eccentricity_function[2][7] = eccentricity_function[0][7] ;
    eccentricity_function[2][8] = eccentricity_function[0][6] ;
    eccentricity_function[2][9] = eccentricity_function[0][5] ;
    eccentricity_function[2][10] =eccentricity_function[0][4] ;
    eccentricity_function[2][11] =eccentricity_function[0][3] ;
    eccentricity_function[2][12] =eccentricity_function[0][2] ;
    eccentricity_function[2][13] =eccentricity_function[0][1] ;
    eccentricity_function[2][14] =eccentricity_function[0][0] ;

    return eccentricity_function;
}

pub fn inclination_function_f2lm(inclination: f64) -> [[f64; 3]; 3]{

    let mut inclination_function = [[0.;3];3];
    // in this notation: 0 is q=-2, 1 is q=-1, 2 is q=0, 3 is q=1, 4 is q=2 ... cf Kaula 64 Table 3

    inclination_function[0][0] = -(3./8.)*(inclination.sin()).powi(2);
    inclination_function[0][1] = (3./4.)*(inclination.sin()).powi(2) -(1./2.);
    inclination_function[0][2] = -(3./8.)*(inclination.sin()).powi(2);

    inclination_function[1][0] = (3./4.)*inclination.sin()*(1. +inclination.cos());
    inclination_function[1][1] = -(3./2.)*inclination.sin()*inclination.cos();
    inclination_function[1][2] = (3./4.)*inclination.sin()*(inclination.cos() -1.);//(3./4.)*inclination.sin()*(1. +inclination.cos());

    inclination_function[2][0] = (3./4.)*(1. +inclination.cos()).powi(2);
    inclination_function[2][1] = (3./2.)*inclination.sin();
    inclination_function[2][2] = (3./4.)*(1. -inclination.cos());
    
    return inclination_function;
}


pub fn kaula_number(wk2:f64, nm_data:f64, real_part_love_number: [[f64;32];32], imaginary_part_love_number: [[f64;32];32], love_number_excitation_frequency: [[f64;32];32] ) -> (f64,f64){
    let mut w_k2 = wk2;
    let mut re_k2 = 0.;
    let mut im_k2 = 0.;
    let mut delta:f64 = 0.;
    let mut ctrl = true;
    let mut find_it = false; 
    let mut parity = false;
    let mut case_oneline_last_col = false;
    let mut more_than_one_col = true;
    let nm_data = nm_data as usize;

    let mut last_line:usize;
    let mut last_column:usize;

    let mut line:usize;
    let mut column:usize;

    // let im_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.imaginary_part_love_number;
    // let re_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.real_part_love_number;
    // let w_lmpq      = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.love_number_excitation_frequency;
    // let nm_data:f64 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.num_datapoints;

    if nm_data%32==0{ // If there is 32 line in the last column
        last_line=31;
        last_column = nm_data/32-1;
        more_than_one_col = false;
    } else {
        last_line = nm_data%32-1;
        last_column = nm_data/32;
        if last_line==0{  // If there is only one line in the last column
            case_oneline_last_col = true;
        }
    }
    
    // println!("\n===============================");
    // println!("\n\tNmbre ={} de points", nm_data);
    // println!("last COL ={} last LINE ={}",last_column,last_line );
    // println!("Excitative Last ={} Lsat-1 ={}  Firts ={} Second ={}",love_number_excitation_frequency[last_column][last_line], love_number_excitation_frequency[last_column][last_line-1], love_number_excitation_frequency[0][0],love_number_excitation_frequency[0][1]);
    // println!("Excitative Lastcol-1 last lign ={} Lastcol first lign = {} ",love_number_excitation_frequency[last_column-1][31], love_number_excitation_frequency[last_column][0]);
    // println!("Love number Last ={} Lsat-1 ={} Firts ={} Second ={}",imaginary_part_love_number[last_column][last_line], imaginary_part_love_number[last_column][last_line-1], imaginary_part_love_number[0][0],imaginary_part_love_number[0][1]);
    // println!("The Lens {} de ses morts", love_number_excitation_frequency.len());
    // println!("Search K2 for w_k2 ={} \n\t",w_k2);

    if w_k2 < 0.0 {
        w_k2 = w_k2.abs();
        parity = true;
    }

    if w_k2 > love_number_excitation_frequency[last_column][last_line-1] { //==================================================================Check beyond the last value of the 2D table
        if case_oneline_last_col {
            re_k2 = real_part_love_number[last_column][last_line]       + (w_k2 - love_number_excitation_frequency[last_column][last_line]) * (real_part_love_number[last_column][last_line] - real_part_love_number[last_column-1][31])/(love_number_excitation_frequency[last_column][last_line]-love_number_excitation_frequency[last_column-1][31]) ;
            im_k2 = imaginary_part_love_number[last_column][last_line]  + (w_k2 - love_number_excitation_frequency[last_column][last_line]) * (imaginary_part_love_number[last_column][last_line] - imaginary_part_love_number[last_column-1][31])/(love_number_excitation_frequency[last_column][last_line]-love_number_excitation_frequency[last_column-1][31]);
        } else {
            delta = (w_k2 - love_number_excitation_frequency[last_column][last_line-1]) / (love_number_excitation_frequency[last_column][last_line-1] - love_number_excitation_frequency[last_column][last_line-2]);
            im_k2 = delta*( imaginary_part_love_number[last_column][last_line-1] -imaginary_part_love_number[last_column][last_line-2] ) +imaginary_part_love_number[last_column][last_line-1];
            re_k2 = delta*( real_part_love_number[last_column][last_line-1] -real_part_love_number[last_column][last_line-2] )           +real_part_love_number[last_column][last_line-1]; 
            // re_k2 = real_part_love_number[last_column][last_line]       + (w_k2 - love_number_excitation_frequency[last_column][last_line]) * (real_part_love_number[last_column][last_line] - real_part_love_number[last_column][last_line-1])/(love_number_excitation_frequency[last_column][last_line]-love_number_excitation_frequency[last_column][last_line-1]) ;
            // im_k2 = imaginary_part_love_number[last_column][last_line]  + (w_k2 - love_number_excitation_frequency[last_column][last_line]) * (imaginary_part_love_number[last_column][last_line] - imaginary_part_love_number[last_column][last_line-1])/(love_number_excitation_frequency[last_column][last_line]-love_number_excitation_frequency[last_column][last_line-1]);
            if im_k2<0. || re_k2<0. { // HAVE TO CHECK
                re_k2 = real_part_love_number[last_column][last_line]       + (real_part_love_number[last_column][last_line] - real_part_love_number[last_column][last_line-1])/2.;
                im_k2 = imaginary_part_love_number[last_column][last_line]  + (imaginary_part_love_number[last_column][last_line] - imaginary_part_love_number[last_column][last_line-1])/2.;
            }
            ctrl = false;
        }
    } else if w_k2 > love_number_excitation_frequency[last_column-1][31] { //================================================================If is in the last column
        if ctrl && w_k2 < love_number_excitation_frequency[last_column][0] { //==============================================================Check is the value is btw the two last column
            // delta = (frequency-W_lmpq(i-1)) / (W_lmpq(i)-W_lmpq(i-1))  realK2 = (1.-delta)*REK2(i-1) + delta*REK2(i)    imaginaryK2 = (1.-delta)*IMK2(i-1) + delta*IMK2(i)
            delta = (w_k2 - love_number_excitation_frequency[last_column-1][31]) / (love_number_excitation_frequency[last_column][0] - love_number_excitation_frequency[last_column-1][31]);
            re_k2 = (1.-delta)* real_part_love_number[last_column-1][31]        +delta*(real_part_love_number[last_column][0]);
            im_k2 = (1.-delta)* imaginary_part_love_number[last_column-1][31]   +delta*(imaginary_part_love_number[last_column][0]);
            // re_k2 = real_part_love_number[last_column][0]        + (real_part_love_number[last_column][0]         - real_part_love_number[last_column-1][31] )     /2.;
            // im_k2 = imaginary_part_love_number[last_column][0]   + (imaginary_part_love_number[last_column][0]    - imaginary_part_love_number[last_column-1][31] )/2.;
            // println!("\n\tBtw the last-1 and the last\t");
            ctrl = false;
        } else if ctrl { //======================================================================================================Check into the last column
            for frequency2 in 0..last_line+1{ //====Because for 00 to last line-1
                if ctrl && love_number_excitation_frequency[last_column][frequency2] >= w_k2 {
                    if  love_number_excitation_frequency[last_column][frequency2] == w_k2 {
                        re_k2 = real_part_love_number[last_column][frequency2];
                        im_k2 = imaginary_part_love_number[last_column][frequency2];
                        ctrl = false;
                        break;
                    } else if ctrl {
                        // delta = (frequency-W_lmpq(i-1)) / (W_lmpq(i)-W_lmpq(i-1))  realK2 = (1.-delta)*REK2(i-1) + delta*REK2(i)    imaginaryK2 = (1.-delta)*IMK2(i-1) + delta*IMK2(i)
                        delta = (w_k2 - love_number_excitation_frequency[last_column][frequency2-1]) / (love_number_excitation_frequency[last_column][frequency2] -love_number_excitation_frequency[last_column][frequency2-1]);
                        re_k2 = (1.-delta)* real_part_love_number[last_column][frequency2-1]        +delta*(real_part_love_number[last_column][frequency2]);
                        im_k2 = (1.-delta)* imaginary_part_love_number[last_column][frequency2-1]   +delta*(imaginary_part_love_number[last_column][frequency2]);
                        // re_k2 = real_part_love_number[last_column][frequency2-1]         + (real_part_love_number[last_column][frequency2]        - real_part_love_number[last_column][frequency2-1])      /2.;
                        // im_k2 = imaginary_part_love_number[last_column][frequency2-1]    + (imaginary_part_love_number[last_column][frequency2]   - imaginary_part_love_number[last_column][frequency2-1]) /2.;
                        // re_k2 = real_part_love_number[last_column][frequency2-1] + (real_part_love_number[last_column][2] - real_part_love_number[last_column][1])/2.;
                        // im_k2 = imaginary_part_love_number[last_column][1] + (imaginary_part_love_number[last_column][2] - imaginary_part_love_number[last_column][1])/2.;
                        // println!("\n\t Into the last col");
                        // println!("Find {}, with {} <-> {}  ", im_k2, imaginary_part_love_number[last_column][frequency2], imaginary_part_love_number[last_column][frequency2-1]);
                        // //println!("Stuff {}", w_k2);
                        ctrl = false;
                        break;
                    }
                    if !ctrl {break;}
                }
                if !ctrl {break;}
            }
        } //-----------------------------End check last column
    } else { //--------------------------Check All other column
        for frequency1 in 0..last_column{ // Loop into the column without the last one (loop from 0 to last_column-1)
            if ctrl && love_number_excitation_frequency[frequency1][31] >= w_k2 { //If the last value of the col frequency1 is GT w2, w2 is in this col
                if ctrl && love_number_excitation_frequency[frequency1][0] > w_k2 { // Check is the fisrt value is GT w2, then the value is btw the two lat colunm
                    // delta = (frequency-W_lmpq(i-1)) / (W_lmpq(i)-W_lmpq(i-1))  realK2 = (1.-delta)*REK2(i-1) + delta*REK2(i)    imaginaryK2 = (1.-delta)*IMK2(i-1) + delta*IMK2(i)
                    delta = (w_k2 - love_number_excitation_frequency[frequency1-1][31]) / (love_number_excitation_frequency[frequency1][0] - love_number_excitation_frequency[frequency1-1][31]);
                    re_k2 = (1.-delta)* real_part_love_number[frequency1-1][31]        +delta*(real_part_love_number[frequency1][0]);
                    im_k2 = (1.-delta)* imaginary_part_love_number[frequency1-1][31]   +delta*(imaginary_part_love_number[frequency1][0]);
                    // re_k2 = real_part_love_number[frequency1][0]        + (real_part_love_number[frequency1][0]         - real_part_love_number[frequency1 -1][31] )        /2.;
                    // im_k2 = imaginary_part_love_number[frequency1][0]   + (imaginary_part_love_number[frequency1][0]    - imaginary_part_love_number[frequency1 -1][31] )   /2.;
                    ctrl = false;
                } else if ctrl {
                    for frequency2 in 0..love_number_excitation_frequency.len(){ //loop into the 32 elements of each line
                        if ctrl && love_number_excitation_frequency[frequency1][frequency2] >= w_k2 { 
                            if  love_number_excitation_frequency[frequency1][frequency2] == w_k2 {
                                re_k2 = real_part_love_number[frequency1][frequency2];
                                im_k2 = imaginary_part_love_number[frequency1][frequency2];
                                ctrl = false;
                                break;
                            }
                            else if ctrl { 
                                delta = (w_k2 - love_number_excitation_frequency[frequency1][frequency2-1]) / (love_number_excitation_frequency[frequency1][frequency2] -love_number_excitation_frequency[frequency1][frequency2-1]);
                                re_k2 = (1.-delta)* real_part_love_number[frequency1][frequency2-1]        +delta*(real_part_love_number[frequency1][frequency2]);
                                im_k2 = (1.-delta)* imaginary_part_love_number[frequency1][frequency2-1]   +delta*(imaginary_part_love_number[frequency1][frequency2]);
                                // re_k2 = real_part_love_number[frequency1][frequency2-1]         + (real_part_love_number[frequency1][frequency2]        - real_part_love_number[frequency1][frequency2-1])      /2.;
                                // im_k2 = imaginary_part_love_number[frequency1][frequency2-1]    + (imaginary_part_love_number[frequency1][frequency2]   - imaginary_part_love_number[frequency1][frequency2-1]) /2.;
                                ctrl = false;
                                break;
                            }
                        }
                        if !ctrl {break;}
                    }
                }
            }
            if !ctrl {break;}
        }
    }
    if parity {
        im_k2 = -im_k2;
        // re_k2 = -re_k2;
    }

    if (w_k2!=0. && im_k2==0.) || (w_k2!=0. && re_k2==0.) {
        println!("Frequ {} \t ImK2 {} \t ReK2 {}", w_k2, im_k2, re_k2);
        panic!("=== Error  in tide.rs fail interpolation k2 ===")
    }

    return (re_k2, im_k2)
}

pub fn sigma_2mpq(m:f64, p:f64, q:f64, spin:f64, orbital_frequency: f64) -> f64{
    // //////println!("\tsigma_calulation: the orbital_frequ {} and spin {}  give the sigma excit_frequ {}, the integer m {}, p {}, q {},  ", orbital_frequency, spin, (2. -2.*p + q)*orbital_frequency -m*spin, m, p, q );
    return (2.0 -2.0*p +q)*orbital_frequency -m*spin;
}

pub fn find_imk2(wk2:f64, nm_data:f64, imaginary_part_love_number: [[f64;32];32], love_number_excitation_frequency: [[f64;32];32] ) -> f64{
    let mut w_k2 = wk2;
    let mut im_k2 = 0.;
    let mut delta:f64 = 0.;
    let mut ctrl = true;
    let mut find_it = false; 
    let mut parity = false;
    let mut case_oneline_last_col = false;
    let mut more_than_one_col = true;
    let nm_data = nm_data as usize;

    let mut last_line:usize;
    let mut last_column:usize;

    let mut line_up:usize=0;
    let mut column_up:usize=0;
    let mut line_down:usize=0;
    let mut column_down:usize=0;

    if nm_data%32==0{ // If there is 32 line in the last column
        last_line=31;
        last_column = nm_data/32-1;
        more_than_one_col = false;
    } else {
        last_line = nm_data%32-1;
        last_column = nm_data/32;
        if last_line==0{  // If there is only one line in the last column
            case_oneline_last_col = true;
        }
    }
    
    if w_k2 < 0.0 {
        w_k2 = w_k2.abs();
        parity = true;
    }

    if w_k2 > love_number_excitation_frequency[last_column][last_line-1] { //==================================================================Check beyond the last value of the 2D table
        if case_oneline_last_col {
            im_k2 = imaginary_part_love_number[last_column][last_line]  + (w_k2 - love_number_excitation_frequency[last_column][last_line]) * (imaginary_part_love_number[last_column][last_line] - imaginary_part_love_number[last_column-1][31])/(love_number_excitation_frequency[last_column][last_line]-love_number_excitation_frequency[last_column-1][31]);
        } else {
            delta = (w_k2 - love_number_excitation_frequency[last_column][last_line-1]) / (love_number_excitation_frequency[last_column][last_line-1] - love_number_excitation_frequency[last_column][last_line-2]);
            im_k2 = delta*( imaginary_part_love_number[last_column][last_line-1] -imaginary_part_love_number[last_column][last_line-2] ) +imaginary_part_love_number[last_column][last_line-1];
            if im_k2<0. { // HAVE TO CHECK
                im_k2 = imaginary_part_love_number[last_column][last_line]  + (imaginary_part_love_number[last_column][last_line] - imaginary_part_love_number[last_column][last_line-1])/2.;
            }
            column_up = last_column;
            line_up = last_line;
            column_down = last_column;
            line_down = last_line;
            ctrl = false;
        }
    } else if w_k2 > love_number_excitation_frequency[last_column-1][31] { //==================================================================If is in the last column
        if ctrl && w_k2 < love_number_excitation_frequency[last_column][0] { //==============================================================Check is the value is btw the two last column
            // delta = (frequency-W_lmpq(i-1)) / (W_lmpq(i)-W_lmpq(i-1))  realK2 = (1.-delta)*REK2(i-1) + delta*REK2(i)    imaginaryK2 = (1.-delta)*IMK2(i-1) + delta*IMK2(i)
            delta = (w_k2 - love_number_excitation_frequency[last_column-1][31]) / (love_number_excitation_frequency[last_column][0] - love_number_excitation_frequency[last_column-1][31]);
            im_k2 = (1.-delta)* imaginary_part_love_number[last_column-1][31]   +delta*(imaginary_part_love_number[last_column][0]);
            column_up  = last_column;
            line_up    = 0;
            column_down  = last_column-1;
            line_down    = 31;
            ctrl = false;
        } else if ctrl { //======================================================================================================Check into the last column
            for frequency2 in 0..last_line+1{ //====Because for 00 to last line-1
                if ctrl && love_number_excitation_frequency[last_column][frequency2] >= w_k2 {
                    if  love_number_excitation_frequency[last_column][frequency2] == w_k2 {
                        im_k2 = imaginary_part_love_number[last_column][frequency2];
                        column_up  = last_column;
                        line_up    = frequency2;
                        column_down  = last_column;
                        line_down    = frequency2;
                        ctrl = false;
                        break;
                    } else if ctrl {
                        // delta = (frequency-W_lmpq(i-1)) / (W_lmpq(i)-W_lmpq(i-1))  realK2 = (1.-delta)*REK2(i-1) + delta*REK2(i)    imaginaryK2 = (1.-delta)*IMK2(i-1) + delta*IMK2(i)
                        delta = (w_k2 - love_number_excitation_frequency[last_column][frequency2-1]) / (love_number_excitation_frequency[last_column][frequency2] -love_number_excitation_frequency[last_column][frequency2-1]);
                        im_k2 = (1.-delta)* imaginary_part_love_number[last_column][frequency2-1]   +delta*(imaginary_part_love_number[last_column][frequency2]);
                        column_up  = last_column;
                        line_up    = frequency2;
                        column_down  = last_column;
                        line_down    = frequency2-1;
                        ctrl = false;
                        break;
                    }
                    if !ctrl {break;}
                }
                if !ctrl {break;}
            }
        } //-----------------------------End check last column
    } else { //================================================================================================================================Check All other column
        for frequency1 in 0..last_column{ // Loop into the column without the last one (loop from 0 to last_column-1)
            if ctrl && love_number_excitation_frequency[frequency1][31] >= w_k2 { //If the last value of the col frequency1 is GT w2, w2 is in this col
                if ctrl && love_number_excitation_frequency[frequency1][0] > w_k2 { // Check is the fisrt value is GT w2, then the value is btw the two lat colunm
                    // delta = (frequency-W_lmpq(i-1)) / (W_lmpq(i)-W_lmpq(i-1))  realK2 = (1.-delta)*REK2(i-1) + delta*REK2(i)    imaginaryK2 = (1.-delta)*IMK2(i-1) + delta*IMK2(i)
                    delta = (w_k2 - love_number_excitation_frequency[frequency1-1][31]) / (love_number_excitation_frequency[frequency1][0] - love_number_excitation_frequency[frequency1-1][31]);
                    im_k2 = (1.-delta)* imaginary_part_love_number[frequency1-1][31]   +delta*(imaginary_part_love_number[frequency1][0]);
                    column_up  = frequency1;
                    line_up    = 0;
                    column_down  = frequency1-1;
                    line_down    = 31;
                    ctrl = false;
                } else if ctrl {
                    for frequency2 in 0..love_number_excitation_frequency.len(){ //loop into the 32 elements of each line
                        if ctrl && love_number_excitation_frequency[frequency1][frequency2] >= w_k2 { 
                            if  love_number_excitation_frequency[frequency1][frequency2] == w_k2 {
                                im_k2 = imaginary_part_love_number[frequency1][frequency2];
                                column_up  = frequency1;
                                line_up    = frequency2;
                                column_down  = frequency1;
                                line_down    = frequency2;
                                ctrl = false;
                                break;
                            }
                            else if ctrl { 
                                delta = (w_k2 - love_number_excitation_frequency[frequency1][frequency2-1]) / (love_number_excitation_frequency[frequency1][frequency2] -love_number_excitation_frequency[frequency1][frequency2-1]);
                                im_k2 = (1.-delta)* imaginary_part_love_number[frequency1][frequency2-1]   +delta*(imaginary_part_love_number[frequency1][frequency2]);
                                column_up  = frequency1;
                                line_up    = frequency2;
                                column_down  = frequency1;
                                line_down    = frequency2-1;
                                ctrl = false;
                                break;
                            }
                        }
                        if !ctrl {break;}
                    }
                }
            }
            if !ctrl {break;}
        }
    }
    if parity {
        im_k2 = -im_k2;
    }

    if (w_k2!=0. && im_k2==0.) || w_k2>love_number_excitation_frequency[column_up][line_up] || w_k2<love_number_excitation_frequency[column_down][line_down] {
        panic!("=== Error  in tide.rs :: fail interpolation k2 === found btw column&line {} {} < {} {} in frequency: {} < {} < {} \t the ImK2 {}", 
        column_down, line_down, column_up, line_up, love_number_excitation_frequency[column_down][line_down], w_k2, love_number_excitation_frequency[column_up][line_up], im_k2);
    }
    return im_k2
}

pub fn find_rek2(wk2:f64, nm_data:f64, real_part_love_number: [[f64;32];32], love_number_excitation_frequency: [[f64;32];32] ) -> f64{
    let mut w_k2 = wk2;
    let mut re_k2 = 0.;
    let mut delta:f64 = 0.;
    let mut ctrl = true;
    let mut find_it = false; 
    let mut parity = false;
    let mut case_oneline_last_col = false;
    let mut more_than_one_col = true;
    let nm_data = nm_data as usize;

    let mut last_line:usize;
    let mut last_column:usize;

    let mut line_up:usize=0;
    let mut column_up:usize=0;
    let mut line_down:usize=0;
    let mut column_down:usize=0;

    if nm_data%32==0{ // If there is 32 line in the last column
        last_line=31;
        last_column = nm_data/32-1;
        more_than_one_col = false;
    } else {
        last_line = nm_data%32-1;
        last_column = nm_data/32;
        if last_line==0{  // If there is only one line in the last column
            case_oneline_last_col = true;
        }
    }
    if w_k2 < 0.0 {
        w_k2 = w_k2.abs(); 
    }

    if w_k2 > love_number_excitation_frequency[last_column][last_line-1] { //==================================================================Check beyond the last value of the 2D table
        if case_oneline_last_col {
            re_k2 = real_part_love_number[last_column][last_line]       + (w_k2 - love_number_excitation_frequency[last_column][last_line]) * (real_part_love_number[last_column][last_line] - real_part_love_number[last_column-1][31])/(love_number_excitation_frequency[last_column][last_line]-love_number_excitation_frequency[last_column-1][31]) ;
        } else {
            delta = (w_k2 - love_number_excitation_frequency[last_column][last_line-1]) / (love_number_excitation_frequency[last_column][last_line-1] - love_number_excitation_frequency[last_column][last_line-2]);
            re_k2 = delta*( real_part_love_number[last_column][last_line-1] -real_part_love_number[last_column][last_line-2] )           +real_part_love_number[last_column][last_line-1];
            if re_k2<0. { // HAVE TO CHECK
                re_k2 = real_part_love_number[last_column][last_line]       + (real_part_love_number[last_column][last_line] - real_part_love_number[last_column][last_line-1])/2.;
            }
            line_up = last_line;
            column_up = last_column;
            line_down = last_line;
            column_down = last_column;
            ctrl = false;
        }
    } else if w_k2 > love_number_excitation_frequency[last_column-1][31] { //==================================================================If is in the last column
        if ctrl && w_k2 < love_number_excitation_frequency[last_column][0] { //==============================================================Check is the value is btw the two last column
            // delta = (frequency-W_lmpq(i-1)) / (W_lmpq(i)-W_lmpq(i-1))  realK2 = (1.-delta)*REK2(i-1) + delta*REK2(i)    imaginaryK2 = (1.-delta)*IMK2(i-1) + delta*IMK2(i)
            delta = (w_k2 - love_number_excitation_frequency[last_column-1][31]) / (love_number_excitation_frequency[last_column][0] - love_number_excitation_frequency[last_column-1][31]);
            re_k2 = (1.-delta)* real_part_love_number[last_column-1][31]        +delta*(real_part_love_number[last_column][0]);
            line_up    = 0;
            column_up  = last_column;
            line_down    = 31;
            column_down  = last_column-1;
            ctrl = false;
        } else if ctrl { //======================================================================================================Check into the last column
            for frequency2 in 0..last_line+1{ //====Because for 00 to last line-1
                if ctrl && love_number_excitation_frequency[last_column][frequency2] >= w_k2 {
                    if  love_number_excitation_frequency[last_column][frequency2] == w_k2 {
                        re_k2 = real_part_love_number[last_column][frequency2];
                        line_up    = frequency2;
                        column_up  = last_column;
                        line_down    = frequency2;
                        column_down  = last_column;
                        ctrl = false;
                        break;
                    } else if ctrl {
                        // delta = (frequency-W_lmpq(i-1)) / (W_lmpq(i)-W_lmpq(i-1))  realK2 = (1.-delta)*REK2(i-1) + delta*REK2(i)    imaginaryK2 = (1.-delta)*IMK2(i-1) + delta*IMK2(i)
                        delta = (w_k2 - love_number_excitation_frequency[last_column][frequency2-1]) / (love_number_excitation_frequency[last_column][frequency2] -love_number_excitation_frequency[last_column][frequency2-1]);
                        re_k2 = (1.-delta)* real_part_love_number[last_column][frequency2-1]        +delta*(real_part_love_number[last_column][frequency2]);
                        line_up    = frequency2;
                        column_up  = last_column;
                        line_up    = frequency2-1;
                        column_up  = last_column;
                        ctrl = false;
                        break;
                    }
                    if !ctrl {break;}
                }
                if !ctrl {break;}
            }
        } //-----------------------------End check last column
    } else { //================================================================================================================================Check All other column
        for frequency1 in 0..last_column{ // Loop into the column without the last one (loop from 0 to last_column-1)
            if ctrl && love_number_excitation_frequency[frequency1][31] >= w_k2 { //If the last value of the col frequency1 is GT w2, w2 is in this col
                if ctrl && love_number_excitation_frequency[frequency1][0] > w_k2 { // Check is the fisrt value is GT w2, then the value is btw the two lat colunm
                    // delta = (frequency-W_lmpq(i-1)) / (W_lmpq(i)-W_lmpq(i-1))  realK2 = (1.-delta)*REK2(i-1) + delta*REK2(i)    imaginaryK2 = (1.-delta)*IMK2(i-1) + delta*IMK2(i)
                    // if frequency1 == 0 { println!("WTF dude why {} is .GT. {}  ", love_number_excitation_frequency[frequency1][0], w_k2);}
                    delta = (w_k2 - love_number_excitation_frequency[frequency1-1][31]) / (love_number_excitation_frequency[frequency1][0] - love_number_excitation_frequency[frequency1-1][31]);
                    re_k2 = (1.-delta)* real_part_love_number[frequency1-1][31]        +delta*(real_part_love_number[frequency1][0]);
                    line_up    = 0;
                    column_up  = frequency1;
                    line_down    = 31;
                    column_down  = frequency1-1;
                    ctrl = false;
                } else if ctrl {
                    for frequency2 in 0..love_number_excitation_frequency.len(){ //loop into the 32 elements of each line
                        if ctrl && love_number_excitation_frequency[frequency1][frequency2] >= w_k2 { 
                            if  love_number_excitation_frequency[frequency1][frequency2] == w_k2 {
                                re_k2 = real_part_love_number[frequency1][frequency2];
                                line_up    = frequency2;
                                column_up  = frequency1;
                                line_down    = frequency2;
                                column_down  = frequency1;
                                ctrl = false;
                                break;
                            }
                            else if ctrl { 
                                delta = (w_k2 - love_number_excitation_frequency[frequency1][frequency2-1]) / (love_number_excitation_frequency[frequency1][frequency2] -love_number_excitation_frequency[frequency1][frequency2-1]);
                                re_k2 = (1.-delta)* real_part_love_number[frequency1][frequency2-1]        +delta*(real_part_love_number[frequency1][frequency2]);
                                line_up    = frequency2;
                                column_up  = frequency1;
                                line_down    = frequency2-1;
                                column_down  = frequency1;
                                ctrl = false;
                                break;
                            }
                        }
                        if !ctrl {break;}
                    }
                }
            }
            if !ctrl {break;}
        }
    }

    if w_k2>love_number_excitation_frequency[column_up][line_up] || w_k2<love_number_excitation_frequency[column_down][line_down] {
        panic!("=== Error  in tide.rs :: fail interpolation k2 === found col {} line {} < {} {}  in frequency: {} < {} < {} \t ReK2 {}", 
        column_down, line_down, column_up, line_up, love_number_excitation_frequency[column_down][line_down], w_k2, love_number_excitation_frequency[column_up][line_up], re_k2);
    }

    return re_k2
}
