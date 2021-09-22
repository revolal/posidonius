use std::collections::HashMap;
use super::super::tools;
use super::super::constants::{K2};
use super::super::{Particle};
use super::super::{Axes};
use super::{EvolutionType};

use crate::constants::{G, TWO_PI,PI, DAY, DEG2RAD};//, AU, M_SUN, M2AU,R_SUN,R_EARTH};//, G_SI, AU, M_SUN};

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
    pub radial_component_of_the_tidal_force_2_kaula: f64,
    pub radial_component_of_the_tidal_force_1_kaula: f64,
    pub radial_component_of_the_tidal_force_1_kaula_nophases: f64,
    pub radial_component_of_the_tidal_force_1_kaula_averaged: f64,
    pub radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass: f64, // Needed to compute denergy_dt
    pub normal_component_of_the_tidal_force: f64,
    pub orthogonal_tidal_force_1_kaula: f64,
    pub orthogonal_tidal_force_1_kaula_averaged: f64,
    pub orthogonal_tidal_force_1_kaula_nophases: f64,
    pub orthogonal_tidal_force_2_kaula: f64,
    pub denergy_dt: f64, // Only for history output
    pub lag_angle: f64, // Used by EvolutionType::BolmontMathis2016, EvolutionType::GalletBolmont2017 and EvolutionType::LeconteChabrier2013(true)

    pub a: f64,
    pub spin: f64,
    pub orbital_frequency: f64,

    pub tidal_force_x: f64,
    pub tidal_force_y: f64,
    pub tidal_force_z: f64,

    pub tidal_torque_x: f64,
    pub tidal_torque_y: f64,
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
                    radial_component_of_the_tidal_force_2_kaula: 0.,
                    radial_component_of_the_tidal_force_1_kaula: 0.,
                    radial_component_of_the_tidal_force_1_kaula_nophases: 0.,
                    radial_component_of_the_tidal_force_1_kaula_averaged: 0.,
                    radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass: 0.,
                    normal_component_of_the_tidal_force: 0., // the theta component of the tidal force, only in the 3D inclined case

                    orthogonal_tidal_force_1_kaula: 0.,
                    orthogonal_tidal_force_1_kaula_averaged: 0.,
                    orthogonal_tidal_force_1_kaula_nophases: 0.,
                    orthogonal_tidal_force_2_kaula: 0.,

                    denergy_dt: 0., // Only for history output
                    lag_angle: 0., // It will be initialized the first time the evolver is called
                
                    a: 0.,
                    spin: 0.,
                    orbital_frequency: 0.,

                    tidal_force_x:0.,
                    tidal_force_y:0.,
                    tidal_force_z:0.,

                    tidal_torque_x: 0.,
                    tidal_torque_y: 0.,
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
            let by_calculation:bool=false;
            let mut torque_due_to_tides_x: f64 = 0.;
            let mut torque_due_to_tides_y: f64 = 0.;
            let mut torque_due_to_tides_z: f64 = 0.;
            if by_calculation {
                // --- The Love numbers data --- //
                let im_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.imaginary_part_love_number;
                let re_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.real_part_love_number;
                let w_lmpq      = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.love_number_excitation_frequency;
                let nm_data:f64 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.num_datapoints;
                // --- The keplerian elements calculation --- //
                let gm = G*(tidal_host_particle.mass+particle.mass);
                let ( semi_major_axis, perihelion_distance, eccentricity, inclination, longitude_perihelion, longitude_of_ascending_node, mean_anomaly, orbital_period, true_anomaly 
                ) = tools::calculate_keplerian_orbital_elements_2( gm, particle.heliocentric_position, particle.heliocentric_velocity);
                let argument_perihelion:f64             = longitude_perihelion -longitude_of_ascending_node ;
                let orbital_frequency:f64               = TWO_PI / (orbital_period *DAY); // in s^-1
                let spin:f64 = ( particle.spin.x.powi(2) + particle.spin.y.powi(2) + particle.spin.z.powi(2) ).sqrt() / DAY;
                let _t = current_time*24.*60.;
                let n = orbital_frequency;
                let l = mean_anomaly;

                // --- The cartesian angles --- //
                // Phi is calculated from the planet centred referecne frame, because of the pi factor
                let phi_planet:f64 = PI+(particle.tides.coordinates.position.y).atan2(particle.tides.coordinates.position.x);
                let cos_phi_planet:f64=phi_planet.cos();
                let sin_phi_planet:f64=phi_planet.sin();
                let cos_theta:f64=0.0_f64;//(0.5*PI).cos();
                let sin_theta:f64=1.0_f64;//(0.5*PI).sin();//(theta is the angle btw the spin axis and the orbital plane)

                
                particle.tides.parameters.internal.a = semi_major_axis;
                particle.tides.parameters.internal.spin = spin;
                particle.tides.parameters.internal.orbital_frequency = orbital_frequency;

                // --- Calculate the principal frequencies for output --- // 
                let mut do_calculation:bool=false;
                if (1.- particle.tides.parameters.internal.check_excitative_frequ / sigma_2mpq(2., 0., 0., spin, orbital_frequency)).abs() > 0.001 {
                    // //println!("Excitative frequency Variation {}", (1.- particle.tides.parameters.internal.check_excitative_frequ / sigma_2mpq(2., 0., 0., spin, orbital_frequency)).abs() );
                    particle.tides.parameters.internal.check_excitative_frequ = sigma_2mpq(2., 0., 0., spin, orbital_frequency);
                    do_calculation = true;
                }
                do_calculation = false;
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
                // --- Eq Torque 1ere implementation --- //
                // let mut sum_g_im_k2_over_q:f64=0.;
                // for x in 0..5{
                //     sum_g_im_k2_over_q = sum_g_im_k2_over_q + ( eccentricity_function[x].powi(2) * particle.tides.parameters.internal.im_love_number_sigma220q[x] ) ;
                // }
                // tidal_torque_kaula_coplanar[1] = - (3./2.)*( ( G* tidal_host_particle.mass.powi(2)* particle.radius.powi(5) ) / semi_major_axis.powi(6) ) *sum_g_im_k2_over_q;

                let mut tidal_torque_kaula_coplanar:f64=0.;
                let mut torque_tangential:f64=0.;
                let mut torque_normal:f64=0.;

                // --- Eq Torque 2em implementation --- //
                let torque_2d:bool=true;
                let torque_3d:bool=false;
                if torque_2d {
                    let eccentricity_function_2pq = eccentricty_function_g(eccentricity);
                    let eccentricity_function = eccentricty_function_g(eccentricity);
                    let associated_legendre_p_2m = associated_legendre_p_2m(0.);

                    let kaula_2_phased:bool=true;
                    let kaula_2_no_phased:bool = false;
                    let kaula_1_phased:bool=false;
                    let kaula_1_averaged_phased:bool=false;
                    let kaula_1_no_phased:bool=false;

                    let mut sum_over_q:f64 = 0.;
                    let mut sum_over_j:f64 = 0.;
                    let mut sum_over_p:f64 = 0.;
                    if kaula_2_phased { // --- Torque 2D  No simplification --- //
                        sum_over_q=0.;
                        for q in -7..8{
                            sum_over_j=0.;
                            for j in -7..8{
                                let phase_term:f64 = ((j as f64)-(q as f64))*l;
                                let frequence:f64 = sigma_2mpq(2., 0., q as f64, spin, orbital_frequency);
                                sum_over_j = sum_over_j +eccentricity_function_2pq[0][(j+7) as usize]*( //sum_over_j = sum_over_j +eccentricity_function_g_20q[(j+2) as usize]*(
                                    phase_term.sin()*find_rek2(frequence, nm_data, re_k2, w_lmpq) +phase_term.cos()*find_imk2(frequence, nm_data, im_k2, w_lmpq) 
                                );
                            }
                            sum_over_q = sum_over_q +eccentricity_function_2pq[0][(q+7) as usize]*sum_over_j; //sum_over_q = sum_over_q +eccentricity_function_g_20q[(q+2) as usize]*sum_over_j;
                        }
                        tidal_torque_kaula_coplanar = -(3./2.)*( (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5) ) / semi_major_axis.powi(6) ) *sum_over_q;
                        torque_normal = -tidal_torque_kaula_coplanar;
                    }

                    if kaula_2_no_phased { // --- Torque 2D  Simplification j = q --- //
                        sum_over_q=0.;
                        for q in -7..8{
                            let frequence:f64 = sigma_2mpq(2., 0., q as f64, spin, orbital_frequency);
                            sum_over_q = sum_over_q +eccentricity_function_2pq[0][(q+7) as usize]*eccentricity_function_2pq[0][(q+7) as usize]*( // eccentricity_function_g_20q[(q+7) as usize]*eccentricity_function_g_20q[(q+7) as usize]*(
                                find_imk2(frequence, nm_data, im_k2, w_lmpq)
                            );
                        }
                        tidal_torque_kaula_coplanar = -(3./2.)*( (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5) ) / semi_major_axis.powi(6) ) *sum_over_q;
                        torque_normal = -tidal_torque_kaula_coplanar;
                    }

                    if kaula_1_phased {
                        let constante:f64 = (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(3)*particle.heliocentric_distance.powi(3)) ;
                        sum_over_p=0.;
                        for p in 0..3{
                            sum_over_q=0.;
                            for q in -7..8{
                                let lambda_22pq = phase_2mpq(2., p as f64, q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_of_ascending_node);
                                let phase_22pq = lambda_22pq -2.*phi_planet;
                                let frequ_22pq = sigma_2mpq(2., p as f64, q as f64, spin, orbital_frequency);

                                sum_over_q = sum_over_q +eccentricity_function[p as usize][(q+7) as usize] * associated_legendre_p_2m[2] *(
                                        phase_22pq.sin()*find_rek2(frequ_22pq, nm_data, re_k2, w_lmpq)
                                        -phase_22pq.cos()*find_imk2(frequ_22pq, nm_data, im_k2, w_lmpq)
                                );
                            }
                            sum_over_p = sum_over_p +sum_over_q;
                        }
                        torque_normal = -(1./2.)*constante *sum_over_p;
                        tidal_torque_kaula_coplanar = -torque_normal
                    }

                    if kaula_1_averaged_phased {
                        let constante:f64 = (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(3)*particle.heliocentric_distance.powi(3)) ;
                        sum_over_p=0.;
                        for p in 0..3{
                            sum_over_q=0.;
                            for q in -7..8{
                                let lambda_22pq = phase_2mpq(2., p as f64, q as f64, 0., 0., argument_perihelion, longitude_of_ascending_node);
                                let phase_22pq = lambda_22pq ;
                                let frequ_22pq = sigma_2mpq(2., p as f64, q as f64, spin, orbital_frequency);

                                sum_over_q = sum_over_q +eccentricity_function[p as usize][(q+7) as usize] * associated_legendre_p_2m[2] *(
                                        phase_22pq.sin()*find_rek2(frequ_22pq, nm_data, re_k2, w_lmpq)
                                        -phase_22pq.cos()*find_imk2(frequ_22pq, nm_data, im_k2, w_lmpq)
                                );
                            }
                            sum_over_p = sum_over_p +sum_over_q;
                        }
                        torque_normal = -(1./2.)*constante *sum_over_p;
                        tidal_torque_kaula_coplanar = -torque_normal
                    }

                    if kaula_1_no_phased {
                        let constante:f64 = (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(3)*particle.heliocentric_distance.powi(3)) ;

                        sum_over_q=0.;
                        for q in -7..8{
                            let frequ_22pq = sigma_2mpq(2., 0., q as f64, spin, orbital_frequency);

                            sum_over_q = sum_over_q +eccentricity_function[0][(q+7) as usize] * associated_legendre_p_2m[2] *(
                                    -find_imk2(frequ_22pq, nm_data, im_k2, w_lmpq)
                            );
                        }

                        torque_normal = -(1./2.)*constante *sum_over_q;
                        tidal_torque_kaula_coplanar = -torque_normal
                        // println!("the normal torque {}", torque_normal);
                    }

                    particle.tides.parameters.internal.tidal_torque_z = tidal_torque_kaula_coplanar;
                    torque_tangential = 0.;
                    // println!("2D \t The normal torque {} \t\t The tangent torque {} \n ", torque_normal, torque_tangential);
                    // torque_normal = -tidal_torque_kaula_coplanar;
                }

                if torque_3d {
                    let constante:f64 = (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(3)*particle.heliocentric_distance.powi(3)) ; //is the theta one 
                    let eccentricity_function = eccentricty_function_g(eccentricity);
                    let inclination_function_f2mp = inclination_function_f2mp(0.); // inclination null
                    let associated_legendre_p_2m = associated_legendre_p_2m(0.); // legendre (cos(theta=\pi/2)=0)) . in coplanar approximation ht declination is null
                    let mut sum_over_p:f64 = 0.;
                    let mut sum_over_q:f64 = 0.;
                    
                    let no_phases:bool=false;
                    let averaged_phases:bool=false;
                    let phased:bool=false;
                    let kaula_2_no_phased:bool=false;
                    let kaula_2_phased:bool=true;

                    if kaula_2_no_phased{
                        let constante_2:f64 = (G*tidal_host_particle.mass.powi(2) *particle.radius.powi(5)) / semi_major_axis.powi(6) ;
                        let mut sum_over_m:f64 = 0.;

                        sum_over_m=0.;
                        for m in 0..3 {
                            // println!("m = {}", m);
                            sum_over_p=0.;
                            for p in 0..3{
                                // println!("\t p = {}", p);
                                sum_over_q=0.;
                                for q in -7..8{
                                    let frequ_2mpq = sigma_2mpq(m as f64, p as f64, q as f64, spin, orbital_frequency);
                                    sum_over_q = sum_over_q +eccentricity_function[p as usize][(q+7) as usize].powi(2) *find_imk2(frequ_2mpq, nm_data, im_k2, w_lmpq) ;
                                    // println!("SUM Q{} = {} ", q, sum_over_q);
                                }
                                sum_over_p = sum_over_p +inclination_function_f2mp[m as usize][p as usize].powi(2) *sum_over_q;
                                // println!("SUM P{} = {} ", p, sum_over_p);
                            }
                            sum_over_m = sum_over_m +sum_over_p *(m as f64) *(2.0_f64 -delta_kron(m)) *factorial(2-m)/factorial(2+m);
                            // println!("SUM M{} = {} ", m, sum_over_m);
                        }
                        torque_normal = constante_2 *sum_over_m;
                        torque_tangential = 0.;
                        // println!(" Normal {} \n", torque_normal);
                    }
                    if kaula_2_phased{
                        let constante_2:f64 = -(G*tidal_host_particle.mass.powi(2) *particle.radius.powi(5)) / semi_major_axis.powi(6) ;
                        let mut sum_over_m:f64 = 0.;
                        let mut sum_over_h:f64 = 0.;
                        let mut sum_over_j:f64 = 0.;

                        sum_over_m=0.;
                        for m in 0..3 {
                            sum_over_p=0.;
                            for p in 0..3{
                                sum_over_q=0.;
                                for q in -7..8{
                                    let frequ_2mpq = sigma_2mpq(m as f64, p as f64, q as f64, spin, orbital_frequency);
                                    sum_over_h=0.;
                                    for h in 0..3{
                                        sum_over_j=0.;
                                        for j in -7..8{
                                            let alpha = (2*(h-p) +(2*h-2*p +q -j)) as f64 ;
                                            sum_over_j = sum_over_j +eccentricity_function[h as usize][(j+7) as usize]*(
                                                alpha.sin() *find_rek2(frequ_2mpq, nm_data, re_k2, w_lmpq)
                                                -alpha.cos() *find_imk2(frequ_2mpq, nm_data, im_k2, w_lmpq)
                                            );
                                        }
                                        sum_over_h = sum_over_h +inclination_function_f2mp[m][h as usize] *sum_over_j ;
                                    }
                                    sum_over_q = sum_over_q +eccentricity_function[p as usize][(q+7) as usize] *sum_over_h ;
                                }
                                sum_over_p = sum_over_p +inclination_function_f2mp[m as usize][p as usize] *sum_over_q;
                            }
                            sum_over_m = sum_over_m +sum_over_p *(m as f64) *(2.0_f64 -delta_kron(m as u64)) *factorial(2-m as u64)/factorial(2+m as u64);
                        }
                        torque_normal = constante_2 *sum_over_m;
                        torque_tangential = 0.;
                    }
                    if no_phases {
                        sum_over_p=0.;
                        for p in 0..3{
                            sum_over_q=0.;
                            for q in -7..8{
                                let frequ_21pq = sigma_2mpq(1., p as f64, q as f64, spin, orbital_frequency);
                                let frequ_22pq = sigma_2mpq(2., p as f64, q as f64, spin, orbital_frequency);
                                // println!("\t the inclination function 21{} = {} \t the inclination function 22{} = {}", p, inclination_function_f2mp[1][p as usize],p, inclination_function_f2mp[2][p as usize] );
                                sum_over_q = sum_over_q +eccentricity_function[p as usize][(q+7) as usize] *(
                                    -(1./3.)*inclination_function_f2mp[1][p as usize] *associated_legendre_p_2m[1] *(
                                        find_rek2(frequ_21pq, nm_data, re_k2, w_lmpq)
                                    )
                                    +(1./6.)*inclination_function_f2mp[2][p as usize] *associated_legendre_p_2m[2] *( 
                                        -find_imk2(frequ_22pq, nm_data, im_k2, w_lmpq)
                                    )
                                );
                            }
                            sum_over_p = sum_over_p +sum_over_q;
                        }
                        torque_normal = -(1./sin_theta) *constante *sum_over_p; 

                        sum_over_p=0.;
                        for p in 0..3{
                            sum_over_q=0.;
                            for q in -7..8{
                                let frequ_20pq = sigma_2mpq(0., p as f64, q as f64, spin, orbital_frequency);
                                let frequ_21pq = sigma_2mpq(1., p as f64, q as f64, spin, orbital_frequency);
                                let frequ_22pq = sigma_2mpq(2., p as f64, q as f64, spin, orbital_frequency);
                                sum_over_q = sum_over_q +eccentricity_function[p as usize][(q+7) as usize] *(
                                    inclination_function_f2mp[0][p as usize] *associated_legendre_p_2m[1] *(
                                        find_rek2(frequ_20pq, nm_data, re_k2, w_lmpq)
                                    )
                                    +(1./6.)  *inclination_function_f2mp[1][p as usize] *(associated_legendre_p_2m[2] -6.*associated_legendre_p_2m[0]) *(
                                        -find_imk2(frequ_21pq, nm_data, im_k2, w_lmpq)
                                    )
                                    -(1./6.) *inclination_function_f2mp[2][p as usize] *associated_legendre_p_2m[1] *(
                                        find_rek2(frequ_22pq, nm_data, re_k2, w_lmpq)
                                    )
                                );
                            }
                            sum_over_p = sum_over_p +sum_over_q;
                        }
                        torque_tangential = constante *sum_over_p;// is phi_planet one

                        // println!("3D \t The normal torque {} \t\t The tangent torque {} \n\n ", torque_normal, torque_tangential);
                    }
                    if averaged_phases {
                        sum_over_p=0.;
                        for p in 0..3{
                            sum_over_q=0.;
                            for q in -7..8{
                                let lambda_21pq = phase_2mpq(1., p as f64, q as f64, 0., 0., argument_perihelion, longitude_of_ascending_node);
                                let lambda_22pq = phase_2mpq(2., p as f64, q as f64, 0., 0., argument_perihelion, longitude_of_ascending_node);
                                let phase_21pq = lambda_21pq ;
                                let phase_22pq = lambda_22pq ;
                                let frequ_21pq = sigma_2mpq(1., p as f64, q as f64, spin, orbital_frequency);
                                let frequ_22pq = sigma_2mpq(2., p as f64, q as f64, spin, orbital_frequency);

                                sum_over_q = sum_over_q +eccentricity_function[p as usize][(q+7) as usize] *(
                                    -(1./3.)*inclination_function_f2mp[1][p as usize] *associated_legendre_p_2m[1] *(
                                        phase_21pq.cos()*find_rek2(frequ_21pq, nm_data, re_k2, w_lmpq)
                                        +phase_21pq.sin()*find_imk2(frequ_21pq, nm_data, im_k2, w_lmpq)
                                    )
                                    +(1./6.)*inclination_function_f2mp[2][p as usize] *associated_legendre_p_2m[2] *( 
                                        phase_22pq.sin()*find_rek2(frequ_22pq, nm_data, re_k2, w_lmpq)
                                        -phase_22pq.cos()*find_imk2(frequ_22pq, nm_data, im_k2, w_lmpq)
                                    )
                                );
                            }
                            sum_over_p = sum_over_p +sum_over_q;
                        }
                        torque_normal = -(1./sin_theta) *constante *sum_over_p; 

                        sum_over_p=0.;
                        for p in 0..3{
                            sum_over_q=0.;
                            for q in -7..8{
                                let lambda_20pq = phase_2mpq(0., p as f64, q as f64, 0., 0., argument_perihelion, longitude_perihelion);
                                let lambda_21pq = phase_2mpq(1., p as f64, q as f64, 0., 0., argument_perihelion, longitude_perihelion);
                                let lambda_22pq = phase_2mpq(2., p as f64, q as f64, 0., 0., argument_perihelion, longitude_perihelion);
                                let phase_20pq = lambda_20pq ;
                                let phase_21pq = lambda_21pq ;
                                let phase_22pq = lambda_22pq ;
                                let frequ_20pq = sigma_2mpq(0., p as f64, q as f64, spin, orbital_frequency);
                                let frequ_21pq = sigma_2mpq(1., p as f64, q as f64, spin, orbital_frequency);
                                let frequ_22pq = sigma_2mpq(2., p as f64, q as f64, spin, orbital_frequency);
                                sum_over_q = sum_over_q +eccentricity_function[p as usize][(q+7) as usize] *(
                                    inclination_function_f2mp[0][p as usize] *associated_legendre_p_2m[1] *(
                                        phase_20pq.cos()*find_rek2(frequ_20pq, nm_data, re_k2, w_lmpq)
                                        +phase_20pq.sin()*find_imk2(frequ_20pq, nm_data, im_k2, w_lmpq)
                                    )
                                    +(1./6.)  *inclination_function_f2mp[1][p as usize] *(associated_legendre_p_2m[2] -6.*associated_legendre_p_2m[0]) *(
                                        phase_21pq.sin()*find_rek2(frequ_21pq, nm_data, re_k2, w_lmpq)
                                        -phase_21pq.cos()*find_imk2(frequ_21pq, nm_data, im_k2, w_lmpq)
                                    )
                                    -(1./6.) *inclination_function_f2mp[2][p as usize] *associated_legendre_p_2m[1] *(
                                        phase_22pq.cos()*find_rek2(frequ_22pq, nm_data, re_k2, w_lmpq)
                                        +phase_22pq.sin()*find_imk2(frequ_22pq, nm_data, im_k2, w_lmpq)
                                    )
                                );
                            }
                            sum_over_p = sum_over_p +sum_over_q;
                        }
                        torque_tangential = constante *sum_over_p;// is phi_planet one
                    }
                    if phased {
                        sum_over_p=0.;
                        for p in 0..3{
                            sum_over_q=0.;
                            for q in -7..8{
                                let lambda_21pq = phase_2mpq(1., p as f64, q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_of_ascending_node);
                                let lambda_22pq = phase_2mpq(2., p as f64, q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_of_ascending_node);
                                let phase_21pq = lambda_21pq -phi_planet;
                                let phase_22pq = lambda_22pq -2.*phi_planet;
                                let frequ_21pq = sigma_2mpq(1., p as f64, q as f64, spin, orbital_frequency);
                                let frequ_22pq = sigma_2mpq(2., p as f64, q as f64, spin, orbital_frequency);

                                sum_over_q = sum_over_q +eccentricity_function[p as usize][(q+7) as usize] *(
                                    -(1./3.)*inclination_function_f2mp[1][p as usize] *associated_legendre_p_2m[1] *(
                                        phase_21pq.cos()*find_rek2(frequ_21pq, nm_data, re_k2, w_lmpq)
                                        +phase_21pq.sin()*find_imk2(frequ_21pq, nm_data, im_k2, w_lmpq)
                                    )
                                    +(1./6.)*inclination_function_f2mp[2][p as usize] *associated_legendre_p_2m[2] *( 
                                        phase_22pq.sin()*find_rek2(frequ_22pq, nm_data, re_k2, w_lmpq)
                                        -phase_22pq.cos()*find_imk2(frequ_22pq, nm_data, im_k2, w_lmpq)
                                    )
                                );
                            }
                            sum_over_p = sum_over_p +sum_over_q;
                        }
                        torque_normal = -(1./sin_theta) *constante *sum_over_p; 

                        sum_over_p=0.;
                        for p in 0..3{
                            sum_over_q=0.;
                            for q in -7..8{
                                let lambda_20pq = phase_2mpq(0., p as f64, q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_perihelion);
                                let lambda_21pq = phase_2mpq(1., p as f64, q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_perihelion);
                                let lambda_22pq = phase_2mpq(2., p as f64, q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_perihelion);
                                let phase_20pq = lambda_20pq;
                                let phase_21pq = lambda_21pq -phi_planet;
                                let phase_22pq = lambda_22pq -2.*phi_planet;
                                let frequ_20pq = sigma_2mpq(0., p as f64, q as f64, spin, orbital_frequency);
                                let frequ_21pq = sigma_2mpq(1., p as f64, q as f64, spin, orbital_frequency);
                                let frequ_22pq = sigma_2mpq(2., p as f64, q as f64, spin, orbital_frequency);
                                sum_over_q = sum_over_q +eccentricity_function[p as usize][(q+7) as usize] *(
                                    inclination_function_f2mp[0][p as usize] *associated_legendre_p_2m[1] *(
                                        phase_20pq.cos()*find_rek2(frequ_20pq, nm_data, re_k2, w_lmpq)
                                        +phase_20pq.sin()*find_imk2(frequ_20pq, nm_data, im_k2, w_lmpq)
                                    )
                                    +(1./6.)  *inclination_function_f2mp[1][p as usize] *(associated_legendre_p_2m[2] -6.*associated_legendre_p_2m[0]) *(
                                        phase_21pq.sin()*find_rek2(frequ_21pq, nm_data, re_k2, w_lmpq)
                                        -phase_21pq.cos()*find_imk2(frequ_21pq, nm_data, im_k2, w_lmpq)
                                    )
                                    -(1./6.) *inclination_function_f2mp[2][p as usize] *associated_legendre_p_2m[1] *(
                                        phase_22pq.cos()*find_rek2(frequ_22pq, nm_data, re_k2, w_lmpq)
                                        +phase_22pq.sin()*find_imk2(frequ_22pq, nm_data, im_k2, w_lmpq)
                                    )
                                );
                            }
                            sum_over_p = sum_over_p +sum_over_q;
                        }
                        torque_tangential = constante *sum_over_p;// is phi_planet one
                    }
                }
                // ---------------------------- //
                // --- The cartesian angles --- //
                let distance:f64 = ( particle.tides.coordinates.position.x.powi(2) +particle.tides.coordinates.position.y.powi(2) +particle.tides.coordinates.position.z.powi(2) ).sqrt();
                let coplanar_distance:f64 = ( particle.tides.coordinates.position.x.powi(2) +particle.tides.coordinates.position.y.powi(2) ).sqrt();
                let cos_theta:f64 = particle.tides.coordinates.position.z / distance;
                let sin_theta:f64 = coplanar_distance / distance;
                let cos_phi:f64 = particle.tides.coordinates.position.x / coplanar_distance;
                let sin_phi:f64 = particle.tides.coordinates.position.y / coplanar_distance;

                // let sin_theta:f64 = (argument_perihelion+true_anomaly).sin() *inclination.sin() ;
                // let tan_phi:f64 = (argument_perihelion+true_anomaly).tan() *inclination.cos();

                // println!("\n Through Cartesian angle {} {} {} {}", cos_theta.acos(), sin_theta.asin(), cos_phi.sin(), sin_phi.asin() );
                // println!(" Through orbtial angle {} {} ",  );

                torque_due_to_tides_x= torque_normal*cos_theta*cos_phi -torque_tangential*sin_phi;
                torque_due_to_tides_y= torque_normal*cos_theta*sin_phi +torque_tangential*cos_phi;
                torque_due_to_tides_z= -torque_normal*sin_theta;
            } else {
                torque_due_to_tides_x = particle.tides.coordinates.position.y *particle.tides.parameters.internal.tidal_force_z 
                                                    -particle.tides.coordinates.position.z *particle.tides.parameters.internal.tidal_force_y  ;
                torque_due_to_tides_y = particle.tides.coordinates.position.z *particle.tides.parameters.internal.tidal_force_x 
                                                    -particle.tides.coordinates.position.x *particle.tides.parameters.internal.tidal_force_z  ;
                torque_due_to_tides_z = particle.tides.coordinates.position.x *particle.tides.parameters.internal.tidal_force_y 
                                                    -particle.tides.coordinates.position.y *particle.tides.parameters.internal.tidal_force_x  ;
            }

            let factor = -1.0;
            if central_body {
                // Integration of the spin (total torque tides):
                dangular_momentum_dt.x += factor * torque_due_to_tides_x;
                dangular_momentum_dt.y += factor * torque_due_to_tides_y;
                dangular_momentum_dt.z += factor * torque_due_to_tides_z;
            } else {
                particle.tides.parameters.output.dangular_momentum_dt.x = factor * torque_due_to_tides_x;
                particle.tides.parameters.output.dangular_momentum_dt.y = factor * torque_due_to_tides_y;
                particle.tides.parameters.output.dangular_momentum_dt.z = factor * torque_due_to_tides_z;
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
                // --- The Love numbers data --- //
                let im_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.imaginary_part_love_number;
                let re_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.real_part_love_number;
                let w_lmpq      = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.love_number_excitation_frequency;
                let nm_data:f64 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.num_datapoints;
                // --- The keplerian elements calculation --- //
                let gm = G*(tidal_host_particle.mass+particle.mass);
                let ( semi_major_axis, _perihelion_distance, eccentricity, inclination, longitude_perihelion, longitude_of_ascending_node, mean_anomaly, orbital_period, true_anomaly
                ) = tools::calculate_keplerian_orbital_elements_2( gm, particle.heliocentric_position, particle.heliocentric_velocity);
                let argument_perihelion:f64             = longitude_perihelion -longitude_of_ascending_node ;
                let orbital_frequency:f64               = TWO_PI / (orbital_period *DAY); // in s^-1
                let spin:f64 = ( particle.spin.x.powi(2) + particle.spin.y.powi(2) + particle.spin.z.powi(2) ).sqrt() / DAY; //this is [rad.s^-1]
                let _n:f64 = orbital_frequency;
                let l:f64 = mean_anomaly;
                // --- The cartesian angle elements --- //
                let phi:f64 = PI+(particle.tides.coordinates.position.y/particle.tides.coordinates.position.x).atan();

                let force_2d:bool=true;
                let force_3d:bool=false;
                if force_2d{
                    // --- Local variable --- //
                    let eccentricity_function_g_20q = eccentricty_function_g_20q(eccentricity);
                    let eccentricity_function_2pq   = eccentricty_function_g(eccentricity);
                    let inclination_function_f2mp   = inclination_function_f2mp(inclination);
                    let associated_legendre_p_2m    = associated_legendre_p_2m(0.); // legendre (cos(theta=\pi/2)=0)) . in coplanar approximation ht declination is null
                    // ---  Is 1_kaula or 2 kaula transformation ? -- //
                    let kaula_1:bool=false;
                    let kaula_1_nophases:bool=false;
                    let kaula_2:bool=true;
                    let kaula_2_no_phased:bool=false;

                    if kaula_1 {
                        // -------------------------------------------------------------------------------------------------------------------------------Eq48 Ortho-Force 1Kaula 2nd implementation
                        let constante:f64= (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(3)*particle.heliocentric_distance.powi(4) );
                        let mut sum_over_q:f64 = 0.;

                        sum_over_q=0.;
                        for q in -2..3{
                            let lambda_220q = phase_2mpq(2., 0., q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_of_ascending_node);
                            let phase_220q = lambda_220q -2.*phi ;
                            let frequ_220q:f64=sigma_2mpq(2., 0., q as f64, spin, orbital_frequency);
                            sum_over_q = sum_over_q +eccentricity_function_g_20q[(q+2) as usize] *( 
                                phase_220q.cos()*find_imk2(frequ_220q, nm_data, im_k2, w_lmpq) 
                                -phase_220q.sin()*find_rek2(frequ_220q, nm_data, re_k2, w_lmpq)  
                            );
                        }
                        let sum_term:f64 = -(1./2.)*associated_legendre_p_2m[2]*sum_over_q;
                        // particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = constante *sum_term;
                        particle.tides.parameters.internal.orthogonal_tidal_force_1_kaula = constante *sum_term;
                    }
                    if kaula_1_nophases {
                        // -------------------------------------------------------------------------------------------------------------------------------Eq48 Ortho-Force 1Kaula 2nd implementation
                        let constante:f64= (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(3)*particle.heliocentric_distance.powi(4) );
                        let mut sum_over_q:f64 = 0.;

                        sum_over_q=0.;
                        for q in -2..3{
                            let phase_220q = 0.0_f64;
                            let frequ_220q:f64=sigma_2mpq(2., 0., q as f64, spin, orbital_frequency);
                            sum_over_q = sum_over_q +eccentricity_function_g_20q[(q+2) as usize] *( 
                                phase_220q.cos()*find_imk2(frequ_220q, nm_data, im_k2, w_lmpq) 
                                -phase_220q.sin()*find_rek2(frequ_220q, nm_data, re_k2, w_lmpq)  
                            );
                        }
                        let sum_term:f64 = -(1./2.)*associated_legendre_p_2m[2]*sum_over_q;
                        // particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = constante *sum_term;
                        particle.tides.parameters.internal.orthogonal_tidal_force_1_kaula_nophases = constante *sum_term;
                    }
                    if kaula_2 {
                        // -------------------------------------------------------------------------------------------------------------------------------Eq48 Ortho-Force 2Kaula 2nd implementation
                        let mut sum_over_q:f64 = 0.;
                        let mut sum_over_j:f64 = 0.;
                        let constante = -(3./2.) *( (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(6)*particle.heliocentric_distance ));

                        sum_over_q = 0.;
                        for q in -7..8{
                            sum_over_j = 0.;
                            for j in -7..8{
                                let phase_term:f64 = ((j as f64)-(q as f64))*l;
                                let frequ_220q:f64 = sigma_2mpq(2., 0., q as f64, spin, orbital_frequency);
                                let (rek2, imk2) = kaula_number(frequ_220q, nm_data, re_k2, im_k2, w_lmpq );
                                sum_over_j = sum_over_j +eccentricity_function_2pq[0][(j+7) as usize] *(
                                    phase_term.sin()*rek2
                                    +phase_term.cos()*imk2  
                                );
                                // sum_over_j = sum_over_j +eccentricity_function_2pq[0][(j+7) as usize] *(
                                //     phase_term.sin()*find_rek2(frequ_220q, nm_data, re_k2, w_lmpq) 
                                //     +phase_term.cos()*find_imk2(frequ_220q, nm_data, im_k2, w_lmpq)  
                                // );
                            }
                            sum_over_q = sum_over_q +eccentricity_function_2pq[0][(q+7) as usize] *sum_over_j;
                        }
                        particle.tides.parameters.internal.orthogonal_tidal_force_2_kaula = constante *sum_over_q;

                    }
                    if kaula_2_no_phased {
                        // -------------------------------------------------------------------------------------------------------------------------------Eq48 Ortho-Force 2Kaula 1er implementation
                        let mut sum_over_q:f64 = 0.;
                        let mut sum_over_j:f64 = 0.;
                        let constante = -(3./2.) *( (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(6)*particle.heliocentric_distance ));
                        sum_over_q = 0.;
                        for q in -7..8{
                            let frequ_220q:f64 = sigma_2mpq(2., 0., q as f64, spin, orbital_frequency);
                            sum_over_q = sum_over_q +eccentricity_function_2pq[0][(q+7) as usize]*eccentricity_function_2pq[0][(q+7) as usize] *(
                                find_imk2(frequ_220q, nm_data, im_k2, w_lmpq)  
                            );
                        }
                        particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = constante *sum_over_q;
                        particle.tides.parameters.internal.orthogonal_tidal_force_2_kaula = constante *sum_over_q;
                    }
                }
                if force_3d{
                    // --- local quantities needed --- //
                    let constante:f64 = (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(3)*particle.heliocentric_distance.powi(4)) ;
                    let inclination_function_f2mp = inclination_function_f2mp(0.); // inclination null
                    let associated_legendre_p_2m  = associated_legendre_p_2m(0.); // legendre (cos(theta=\pi/2)=0)) . in coplanar approximation ht declination is null
                    let eccentricity_function_2pq = eccentricty_function_g(eccentricity);

                    let sin_theta:f64=(0.5*PI).sin();//(theta is the angle btw the spin axis and the orbital plane)

                    // --- summation terms --- //
                    let mut sum_over_p:f64 = 0.;
                    let mut sum_over_q:f64 = 0.;
                    let kaula_1:bool=true;
                    let kaula_1_averaged_phases:bool=true;
                    let kaula_1_nophases:bool=true;
                    let kaula_2:bool=false;
                    if kaula_1{
                        sum_over_p=0.;
                        for p in 0..3{
                            sum_over_q=0.;
                            for q in -7..8{ //for q in -2..3{
                                let lambda_21pq = phase_2mpq(1., p as f64, q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_of_ascending_node);
                                let lambda_22pq = phase_2mpq(2., p as f64, q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_of_ascending_node);
                                let phase_21pq = lambda_21pq -phi ;
                                let phase_22pq = lambda_22pq -2.*phi ;
                                let frequ_21pq = sigma_2mpq(1., p as f64, q as f64, spin, orbital_frequency);
                                let frequ_22pq = sigma_2mpq(2., p as f64, q as f64, spin, orbital_frequency);

                                sum_over_q = sum_over_q +eccentricity_function_2pq[p as usize][(q+7) as usize] *(
                                    -(1./3.)*inclination_function_f2mp[1][p as usize] *associated_legendre_p_2m[1] *(
                                        phase_21pq.cos()*find_rek2(frequ_21pq, nm_data, re_k2, w_lmpq)
                                        +phase_21pq.sin()*find_imk2(frequ_21pq, nm_data, im_k2, w_lmpq)
                                    )
                                    +(1./6.)*inclination_function_f2mp[2][p as usize] *associated_legendre_p_2m[2] *( 
                                        phase_22pq.sin()*find_rek2(frequ_22pq, nm_data, re_k2, w_lmpq)
                                        -phase_22pq.cos()*find_imk2(frequ_22pq, nm_data, im_k2, w_lmpq)
                                    )
                                );
                            }
                            sum_over_p = sum_over_p +sum_over_q;
                        }
                        // particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = (1./sin_theta) *constante *sum_over_p;
                        particle.tides.parameters.internal.orthogonal_tidal_force_1_kaula =  (1./sin_theta) *constante *sum_over_p;
                    }
                    if kaula_1_averaged_phases{
                        sum_over_p=0.;
                        for p in 0..3{
                            sum_over_q=0.;
                            for q in -7..8{ //for q in -2..3{
                                let lambda_21pq = phase_2mpq(1., p as f64, q as f64, 0., 0., argument_perihelion, longitude_of_ascending_node);
                                let lambda_22pq = phase_2mpq(2., p as f64, q as f64, 0., 0., argument_perihelion, longitude_of_ascending_node);
                                let phase_21pq = lambda_21pq; // -phi ;
                                let phase_22pq = lambda_22pq; // -2.*phi ;
                                let frequ_21pq = sigma_2mpq(1., p as f64, q as f64, spin, orbital_frequency);
                                let frequ_22pq = sigma_2mpq(2., p as f64, q as f64, spin, orbital_frequency);

                                sum_over_q = sum_over_q +eccentricity_function_2pq[p as usize][(q+7) as usize] *(
                                    -(1./3.)*inclination_function_f2mp[1][p as usize] *associated_legendre_p_2m[1] *(
                                        phase_21pq.cos()*find_rek2(frequ_21pq, nm_data, re_k2, w_lmpq)
                                        +phase_21pq.sin()*find_imk2(frequ_21pq, nm_data, im_k2, w_lmpq)
                                    )
                                    +(1./6.)*inclination_function_f2mp[2][p as usize] *associated_legendre_p_2m[2] *( 
                                        phase_22pq.sin()*find_rek2(frequ_22pq, nm_data, re_k2, w_lmpq)
                                        -phase_22pq.cos()*find_imk2(frequ_22pq, nm_data, im_k2, w_lmpq)
                                    )
                                );
                            }
                            sum_over_p = sum_over_p +sum_over_q;
                        }
                        // particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = (1./sin_theta) *constante *sum_over_p;
                        particle.tides.parameters.internal.orthogonal_tidal_force_1_kaula_averaged =  (1./sin_theta) *constante *sum_over_p;
                    }
                    if kaula_1_nophases{
                        sum_over_p=0.;
                        for p in 0..3{
                            sum_over_q=0.;
                            for q in -7..8{ //for q in -2..3{
                                let phase_21pq = 0.0_f64;
                                let phase_22pq = 0.0_f64;
                                let frequ_21pq = sigma_2mpq(1., p as f64, q as f64, spin, orbital_frequency);
                                let frequ_22pq = sigma_2mpq(2., p as f64, q as f64, spin, orbital_frequency);

                                sum_over_q = sum_over_q +eccentricity_function_2pq[p as usize][(q+7) as usize] *(
                                    -(1./3.)*inclination_function_f2mp[1][p as usize] *associated_legendre_p_2m[1] *(
                                        phase_21pq.cos()*find_rek2(frequ_21pq, nm_data, re_k2, w_lmpq)
                                        +phase_21pq.sin()*find_imk2(frequ_21pq, nm_data, im_k2, w_lmpq)
                                    )
                                    +(1./6.)*inclination_function_f2mp[2][p as usize] *associated_legendre_p_2m[2] *( 
                                        phase_22pq.sin()*find_rek2(frequ_22pq, nm_data, re_k2, w_lmpq)
                                        -phase_22pq.cos()*find_imk2(frequ_22pq, nm_data, im_k2, w_lmpq)
                                    )
                                );
                            }
                            sum_over_p = sum_over_p +sum_over_q;
                        }
                        // particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = (1./sin_theta) *constante *sum_over_p;
                        particle.tides.parameters.internal.orthogonal_tidal_force_1_kaula_nophases = (1./sin_theta) *constante *sum_over_p;
                    }
                    if kaula_2{

                    }
                    // println!("Ortho force 3D {:0.1e}",particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide );
                }
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
            // --- The Love numbers data --- //
            let im_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.imaginary_part_love_number;
            let re_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.real_part_love_number;
            let w_lmpq      = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.love_number_excitation_frequency;
            let nm_data:f64 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.num_datapoints;
            // --- Calculate keplerian elements --- //
            let gm = G*(tidal_host_particle.mass+particle.mass);
            let ( semi_major_axis, _perihelion_distance, eccentricity, inclination, longitude_perihelion, longitude_of_ascending_node, mean_anomaly, orbital_period, true_anomaly 
            ) = tools::calculate_keplerian_orbital_elements_2( gm, particle.heliocentric_position, particle.heliocentric_velocity);
            let argument_perihelion:f64             = longitude_perihelion -longitude_of_ascending_node ;
            let orbital_frequency:f64               = TWO_PI / (orbital_period *DAY); // in s^-1
            let spin:f64                            = ( particle.spin.x.powi(2) + particle.spin.y.powi(2) + particle.spin.z.powi(2) ).sqrt() / DAY; //this is [rad.s^-1]
            let l:f64 = mean_anomaly;
            let _n:f64 = orbital_frequency;
            let time:f64 = current_time*24.*60.; // in seconds
            // --- The cartesian angle elements --- //
            let phi:f64 = PI+(particle.tides.coordinates.position.y).atan2(particle.tides.coordinates.position.x);

            let force_2d:bool=true;
            let force_3d:bool=false;
            if force_2d {
                // --- Local quantities --- //
                let eccentricity_function_g_20q = eccentricty_function_g_20q(eccentricity); //eccentricity_function_2pq[0][(q+7) as usize]
                let eccentricity_function_g_21q = eccentricty_function_g_21q(eccentricity); //eccentricity_function_2pq[1][(q+7) as usize]
                let eccentricity_function_2pq = eccentricty_function_g(eccentricity);
                let associated_legendre_p_2m    = associated_legendre_p_2m(0.); //0. in coplanar approximation
                // ---  Is 1_kaula or 2 kaula transformation ? -- //
                let kaula_1:bool=false;
                let kaula_1_nophases:bool=false;
                let kaula_2:bool=true;
                let kaula_2_no_phased:bool=false;
                if kaula_1 { 
                    // println("{} 1 Kaula", '\n');
                    // ------------------------------------------------------------------------------------------------------------------------------------Radial-Force 1Kaula 2nd implementation
                    let mut sum_over_q:f64=0.;
                    let constante = (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(3) *particle.heliocentric_distance.powi(4)) ;
                    sum_over_q=0.;
                    for q in -2..3{
                        let lambda_201q = phase_2mpq(0., 1., q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_of_ascending_node);
                        let phase_201q = lambda_201q;
                        let frequ_201q = sigma_2mpq(0., 1., q as f64, spin, orbital_frequency);
                        sum_over_q = sum_over_q +eccentricity_function_g_21q[(q+2) as usize]*( 
                            phase_201q.cos()*find_rek2(frequ_201q, nm_data, re_k2, w_lmpq) 
                            +phase_201q.sin()*find_imk2(frequ_201q, nm_data, im_k2, w_lmpq)  
                        );
                    }
                    let term1:f64 = (3./2.)*associated_legendre_p_2m[0]*sum_over_q;
                    // println!("1er term {} ", term1);

                    sum_over_q=0.;
                    for q in -2..3{
                        let lambda_220q = phase_2mpq(2., 0., q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_of_ascending_node);
                        let phase_220q = lambda_220q -2.*phi;
                        let frequ_220q = sigma_2mpq(2., 0., q as f64, spin, orbital_frequency);
                        sum_over_q = sum_over_q +eccentricity_function_g_20q[(q+2) as usize]*( 
                            phase_220q.cos()*find_rek2(frequ_220q, nm_data, re_k2, w_lmpq) 
                            +phase_220q.sin()*find_imk2(frequ_220q, nm_data, im_k2, w_lmpq)  
                        );
                    }
                    let term2:f64 = -(3./4.)*associated_legendre_p_2m[2]*sum_over_q;
                    // println!("2em term {} ", term2);

                    particle.tides.parameters.internal.radial_component_of_the_tidal_force_1_kaula = constante *(term1+term2);
                    // particle.tides.parameters.internal.radial_component_of_the_tidal_force = constante *(term1+term2);
                }
                if kaula_1_nophases { 
                    // println("{} 1 Kaula", '\n');
                    // ------------------------------------------------------------------------------------------------------------------------------------Radial-Force 1Kaula 2nd implementation
                    let mut sum_over_q:f64=0.;
                    let constante = (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(3) *particle.heliocentric_distance.powi(4)) ;
                    sum_over_q=0.;
                    for q in -2..3{
                        let _lambda_201q = phase_2mpq(0., 1., q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_of_ascending_node);
                        let phase_201q = 0.0_f64;//lambda_201q;
                        let frequ_201q = sigma_2mpq(0., 1., q as f64, spin, orbital_frequency);
                        sum_over_q = sum_over_q +eccentricity_function_g_21q[(q+2) as usize]*( 
                            phase_201q.cos()*find_rek2(frequ_201q, nm_data, re_k2, w_lmpq) 
                            +phase_201q.sin()*find_imk2(frequ_201q, nm_data, im_k2, w_lmpq)  
                        );
                    }
                    let term1:f64 = (3./2.)*associated_legendre_p_2m[0]*sum_over_q;
                    // println!("1er term {} ", term1);

                    sum_over_q=0.;
                    for q in -2..3{
                        let _lambda_220q = phase_2mpq(2., 0., q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_of_ascending_node);
                        let phase_220q = 0.0_f64;//lambda_220q -2.*phi;
                        let frequ_220q = sigma_2mpq(2., 0., q as f64, spin, orbital_frequency);
                        sum_over_q = sum_over_q +eccentricity_function_g_20q[(q+2) as usize]*( 
                            phase_220q.cos()*find_rek2(frequ_220q, nm_data, re_k2, w_lmpq) 
                            +phase_220q.sin()*find_imk2(frequ_220q, nm_data, im_k2, w_lmpq)  
                        );
                    }
                    let term2:f64 = -(3./4.)*associated_legendre_p_2m[2]*sum_over_q;
                    // println!("2em term {} ", term2);

                    particle.tides.parameters.internal.radial_component_of_the_tidal_force_1_kaula_nophases = constante *(term1+term2);
                    // particle.tides.parameters.internal.radial_component_of_the_tidal_force = constante *(term1+term2);
                }
                if kaula_2 {
                    // ------------------------------------------------------------------------------------------------------------------------------------Eq48 Radial-Force 2Kaula 2nd implementation
                    let mut sum_over_q:f64 = 0.;
                    let mut sum_over_j:f64 = 0.;
                    let constante = -( (G *tidal_host_particle.mass.powi(2) *particle.radius.powi(5)) / (semi_major_axis.powi(6) *particle.heliocentric_distance) );
                    sum_over_q = 0.;
                    for q in -7..8{
                        sum_over_j = 0.;
                        for j in -7..8{
                            let phase_term:f64 = ((j as f64)-(q as f64))*l;
                            let frequence:f64 = sigma_2mpq(0., 1., q as f64, spin, orbital_frequency);
                            let (rek2, imk2) = kaula_number(frequence, nm_data, re_k2, im_k2, w_lmpq );
                            sum_over_j = sum_over_j +eccentricity_function_2pq[1][(j+7) as usize]*( 
                                phase_term.cos()*rek2 
                                -phase_term.sin()*imk2  
                            );
                            // sum_over_j = sum_over_j +eccentricity_function_2pq[1][(j+7) as usize]*( 
                            //     phase_term.cos()*find_rek2(frequence, nm_data, re_k2, w_lmpq) 
                            //     -phase_term.sin()*find_imk2(frequence, nm_data, im_k2, w_lmpq)  
                            // );
                        }
                        sum_over_q = sum_over_q +(3./4.)*eccentricity_function_2pq[1][(q+7) as usize]*sum_over_j;
                    }
                    let term1:f64 = sum_over_q;
                    sum_over_q = 0.;
                    for q in -7..8{
                        sum_over_j = 0.;
                        for j in -7..8{
                            let phase_term:f64 = ((j as f64)-(q as f64))*l;
                            let frequence:f64 = sigma_2mpq(2., 0., q as f64, spin, orbital_frequency);
                            let (rek2, imk2) = kaula_number(frequence, nm_data, re_k2, im_k2, w_lmpq );
                            sum_over_j = sum_over_j +eccentricity_function_2pq[0][(j+7) as usize]*( 
                                phase_term.cos()*rek2
                                -phase_term.sin()*imk2
                            );
                            // sum_over_j = sum_over_j +eccentricity_function_2pq[0][(j+7) as usize]*( 
                            //     phase_term.cos()*find_rek2(frequence, nm_data, re_k2, w_lmpq)
                            //     -phase_term.sin()*find_imk2(frequence, nm_data, im_k2, w_lmpq)
                            // );
                        }
                        sum_over_q = sum_over_q +(9./4.)*eccentricity_function_2pq[0][(q+7) as usize]*sum_over_j;
                    }
                    let term2:f64 = sum_over_q;
                    particle.tides.parameters.internal.radial_component_of_the_tidal_force_2_kaula = constante *(term1+term2);

                }
                if kaula_2_no_phased {
                    // ------------------------------------------------------------------------------------------------------------------------------------Eq48 Radial-Force 2Kaula 1er implementation
                    let mut sum_over_q:f64 = 0.;
                    let mut sum_over_j:f64 = 0.;
                    let constante = -( (G *tidal_host_particle.mass.powi(2) *particle.radius.powi(5)) / (semi_major_axis.powi(6) *particle.heliocentric_distance) );
                    sum_over_q = 0.;
                    for q in -7..8{
                        let frequence:f64 = sigma_2mpq(0., 1., q as f64, spin, orbital_frequency);
                        sum_over_q = sum_over_q +(3./4.)*eccentricity_function_2pq[1][(q+7) as usize] *eccentricity_function_2pq[1][(q+7) as usize]*( 
                            find_rek2(frequence, nm_data, re_k2, w_lmpq) 
                        );
                    }
                    let term1:f64 = sum_over_q;
                    sum_over_q = 0.;
                    for q in -7..8{
                        let frequence:f64 = sigma_2mpq(2., 0., q as f64, spin, orbital_frequency);
                        sum_over_q = sum_over_q +(9./4.)*eccentricity_function_2pq[0][(q+7) as usize]*eccentricity_function_2pq[0][(q+7) as usize]*( 
                            find_rek2(frequence, nm_data, re_k2, w_lmpq)
                        );
                    }
                    let term2:f64 = sum_over_q;
                    particle.tides.parameters.internal.radial_component_of_the_tidal_force_2_kaula = constante *(term1+term2);
                }
            }
            if force_3d {
                // --- local quantities needed --- //
                let constante:f64 = (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(3)*particle.heliocentric_distance.powi(4)) ;
                let inclination_function_f2mp = inclination_function_f2mp(0.); // inclination null
                let associated_legendre_p_2m  = associated_legendre_p_2m(0.); // legendre (cos(theta=\pi/2)=0)) . in coplanar approximation ht declination is null
                let eccentricity_function_2pq = eccentricty_function_g(eccentricity);

                // --- summation terms --- //
                let mut sum_over_p:f64 = 0.;
                let mut sum_over_q:f64 = 0.;

                let kaula_1:bool = true;
                let kaula_averged_phases_1:bool = true;
                let kaula_no_phases_1:bool = true;

                if kaula_1{
                    // --- 1 Kaula --- //
                    sum_over_p=0.;
                    for p in 0..3{
                        sum_over_q=0.;
                        for q in -7..8{ //for q in -2..3{
                            let lambda_20pq = phase_2mpq(0., p as f64, q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_of_ascending_node);
                            let lambda_21pq = phase_2mpq(1., p as f64, q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_of_ascending_node);
                            let lambda_22pq = phase_2mpq(2., p as f64, q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_of_ascending_node);
                            let phase_20pq = lambda_20pq ;
                            let phase_21pq = lambda_21pq -phi;
                            let phase_22pq = lambda_22pq  -2.*phi;
                            let frequ_20pq = sigma_2mpq(0., p as f64, q as f64, spin, orbital_frequency);
                            let frequ_21pq = sigma_2mpq(1., p as f64, q as f64, spin, orbital_frequency);
                            let frequ_22pq = sigma_2mpq(2., p as f64, q as f64, spin, orbital_frequency);

                            sum_over_q = sum_over_q -3.0*eccentricity_function_2pq[p as usize][(q+7) as usize] *(
                                inclination_function_f2mp[0][p as usize] *associated_legendre_p_2m[0] *(
                                    phase_20pq.cos()*find_rek2(frequ_20pq, nm_data, re_k2, w_lmpq)
                                    +phase_20pq.sin()*find_imk2(frequ_20pq, nm_data, im_k2, w_lmpq)
                                )
                                +(1./3.)*inclination_function_f2mp[1][p as usize] *associated_legendre_p_2m[1] *( 
                                    phase_21pq.sin()*find_rek2(frequ_21pq, nm_data, re_k2, w_lmpq)
                                    -phase_21pq.cos()*find_imk2(frequ_21pq, nm_data, im_k2, w_lmpq)
                                )
                                +(1./12.)*inclination_function_f2mp[2][p as usize] *associated_legendre_p_2m[2] *(
                                    phase_22pq.cos()*find_rek2(frequ_22pq, nm_data, re_k2, w_lmpq)
                                    +phase_22pq.sin()*find_imk2(frequ_22pq, nm_data, im_k2, w_lmpq)
                                )
                            );
                        }
                        sum_over_p = sum_over_p +sum_over_q;
                    }
                    particle.tides.parameters.internal.radial_component_of_the_tidal_force_1_kaula = constante *sum_over_p;
                }
                if kaula_averged_phases_1 {
                    // --- 1 Kaula Averaged phase --- //
                    sum_over_p=0.;
                    for p in 0..3{
                        sum_over_q=0.;
                        for q in -7..8{ //for q in -2..3{
                            let lambda_20pq = phase_2mpq(0., p as f64, q as f64, 0., 0., argument_perihelion, longitude_of_ascending_node);
                            let lambda_21pq = phase_2mpq(1., p as f64, q as f64, 0., 0., argument_perihelion, longitude_of_ascending_node);
                            let lambda_22pq = phase_2mpq(2., p as f64, q as f64, 0., 0., argument_perihelion, longitude_of_ascending_node);
                            let phase_20pq = lambda_20pq ;
                            let phase_21pq = lambda_21pq ; // -phi;
                            let phase_22pq = lambda_22pq ; // -2.*phi;
                            let frequ_20pq = sigma_2mpq(0., p as f64, q as f64, spin, orbital_frequency);
                            let frequ_21pq = sigma_2mpq(1., p as f64, q as f64, spin, orbital_frequency);
                            let frequ_22pq = sigma_2mpq(2., p as f64, q as f64, spin, orbital_frequency);

                            sum_over_q = sum_over_q -3.0*eccentricity_function_2pq[p as usize][(q+7) as usize] *(
                                inclination_function_f2mp[0][p as usize] *associated_legendre_p_2m[0] *(
                                    phase_20pq.cos()*find_rek2(frequ_20pq, nm_data, re_k2, w_lmpq)
                                    +phase_20pq.sin()*find_imk2(frequ_20pq, nm_data, im_k2, w_lmpq)
                                )
                                +(1./3.)*inclination_function_f2mp[1][p as usize] *associated_legendre_p_2m[1] *( 
                                    phase_21pq.sin()*find_rek2(frequ_21pq, nm_data, re_k2, w_lmpq)
                                    -phase_21pq.cos()*find_imk2(frequ_21pq, nm_data, im_k2, w_lmpq)
                                )
                                +(1./12.)*inclination_function_f2mp[2][p as usize] *associated_legendre_p_2m[2] *(
                                    phase_22pq.cos()*find_rek2(frequ_22pq, nm_data, re_k2, w_lmpq)
                                    +phase_22pq.sin()*find_imk2(frequ_22pq, nm_data, im_k2, w_lmpq)
                                )
                            );
                        }
                        sum_over_p = sum_over_p +sum_over_q;
                    }
                    particle.tides.parameters.internal.radial_component_of_the_tidal_force_1_kaula_averaged = constante *sum_over_p;
                }

                if kaula_no_phases_1 {
                    // --- 1 Kaula No Phases --- //
                    sum_over_p=0.;
                    for p in 0..3{
                        sum_over_q=0.;
                        for q in -7..8{ //for q in -2..3{
                            let phase_20pq = 0.0_f64;
                            let phase_21pq = 0.0_f64;
                            let phase_22pq = 0.0_f64;
                            let frequ_20pq = sigma_2mpq(0., p as f64, q as f64, spin, orbital_frequency);
                            let frequ_21pq = sigma_2mpq(1., p as f64, q as f64, spin, orbital_frequency);
                            let frequ_22pq = sigma_2mpq(2., p as f64, q as f64, spin, orbital_frequency);

                            sum_over_q = sum_over_q -3.0*eccentricity_function_2pq[p as usize][(q+7) as usize] *(
                                inclination_function_f2mp[0][p as usize] *associated_legendre_p_2m[0] *(
                                    phase_20pq.cos()*find_rek2(frequ_20pq, nm_data, re_k2, w_lmpq)
                                    +phase_20pq.sin()*find_imk2(frequ_20pq, nm_data, im_k2, w_lmpq)
                                )
                                +(1./3.)*inclination_function_f2mp[1][p as usize] *associated_legendre_p_2m[1] *( 
                                    phase_21pq.sin()*find_rek2(frequ_21pq, nm_data, re_k2, w_lmpq)
                                    -phase_21pq.cos()*find_imk2(frequ_21pq, nm_data, im_k2, w_lmpq)
                                )
                                +(1./12.)*inclination_function_f2mp[2][p as usize] *associated_legendre_p_2m[2] *(
                                    phase_22pq.cos()*find_rek2(frequ_22pq, nm_data, re_k2, w_lmpq)
                                    +phase_22pq.sin()*find_imk2(frequ_22pq, nm_data, im_k2, w_lmpq)
                                )
                            );
                        }
                        sum_over_p = sum_over_p +sum_over_q;
                    }
                    particle.tides.parameters.internal.radial_component_of_the_tidal_force_1_kaula_nophases = constante *sum_over_p;
                }

            }
        }
    }
}

pub fn calculate_normal_component_of_the_tidal_force(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>, current_time:f64) {

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody = particle.tides.effect {
            let force_normal = 0.;
        }

        if let TidesEffect::KaulaCoplanarOrbitingBody = particle.tides.effect {

            // --- The computation of the normal force (the theta component) --- //
            let mut force_normal:f64=0.;
            let force_2d:bool=true;
            let force_3d:bool=false;
            if force_2d {
                force_normal =0.;
            }
            if force_3d {
                // --- The Love numbers data --- //
                let im_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.imaginary_part_love_number;
                let re_k2       = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.real_part_love_number;
                let w_lmpq      = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.love_number_excitation_frequency;
                let nm_data:f64 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.num_datapoints;
                
                // --- The keplerian elements --- //
                let gm = G*(tidal_host_particle.mass+particle.mass);
                let ( semi_major_axis, _perihelion_distance, eccentricity, inclination, longitude_perihelion, longitude_of_ascending_node, mean_anomaly, orbital_period, true_anomaly 
                ) = tools::calculate_keplerian_orbital_elements_2( gm, particle.heliocentric_position, particle.heliocentric_velocity);
                let argument_perihelion:f64 = longitude_perihelion -longitude_of_ascending_node ;
                let orbital_frequency:f64   = TWO_PI / (orbital_period *DAY); // in s^-1
                let spin:f64                = ( particle.spin.x.powi(2) + particle.spin.y.powi(2) + particle.spin.z.powi(2) ).sqrt() / DAY; //this is [rad.s^-1]
                let l:f64 = mean_anomaly;
                let n:f64 = orbital_frequency;
                let time:f64 = current_time*24.*60.; //days in seconds
                // --- cartesian angle calculation --- /
                // let declination:f64 = PI-(particle.tides.coordinates.position.y).atan2(particle.tides.coordinates.position.x);
                // let declination:f64 = -PI+(particle.tides.coordinates.position.y).atan2(particle.tides.coordinates.position.x);
                let phi:f64 = PI-(particle.tides.coordinates.position.y).atan2(particle.tides.coordinates.position.x);
                // let phi:f64 = declination - spin*current_time*24.*60.;
                let cos_phi:f64=phi.cos();
                let sin_phi:f64=phi.sin();
                let cos_theta:f64=(0.5*PI).cos();
                let sin_theta:f64=(0.5*PI).sin();//(theta is the angle btw the spin axis and the orbital plane)

                let constante:f64 = (G*tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(3)*particle.heliocentric_distance.powi(4)) ;
                let eccentricity_function_2pq = eccentricty_function_g(eccentricity);
                let inclination_function_f2mp = inclination_function_f2mp(0.); // inclination null
                let associated_legendre_p_2m  = associated_legendre_p_2m(0.); // legendre (cos(theta=\pi/2)=0)) . in coplanar approximation ht declination is null
                // --- summation terms --- //
                let mut sum_over_p:f64 = 0.;
                let mut sum_over_q:f64 = 0.;
                
                sum_over_p=0.;
                for p in 0..3{
                    sum_over_q=0.;
                    for q in -2..3{
                        let lambda_20pq = phase_2mpq(0., p as f64, q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_of_ascending_node);
                        let lambda_21pq = phase_2mpq(1., p as f64, q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_of_ascending_node);
                        let lambda_22pq = phase_2mpq(2., p as f64, q as f64, mean_anomaly, spin*current_time*24.*60., argument_perihelion, longitude_of_ascending_node);
                        let phase_20pq = 0.0_f64;//lambda_20pq ;
                        let phase_21pq = 0.0_f64;//lambda_21pq -phi ;
                        let phase_22pq = 0.0_f64;//lambda_22pq -2.*phi ;
                        let frequ_20pq = sigma_2mpq(0., p as f64, q as f64, spin, orbital_frequency);
                        let frequ_21pq = sigma_2mpq(1., p as f64, q as f64, spin, orbital_frequency);
                        let frequ_22pq = sigma_2mpq(2., p as f64, q as f64, spin, orbital_frequency);

                        sum_over_q = sum_over_q +eccentricity_function_2pq[p as usize][(q+7) as usize] *(
                            inclination_function_f2mp[0][p as usize] *associated_legendre_p_2m[1] *(
                                phase_20pq.cos()*find_rek2(frequ_20pq, nm_data, re_k2, w_lmpq)
                                +phase_20pq.sin()*find_imk2(frequ_20pq, nm_data, im_k2, w_lmpq)
                            )
                            +(1./6.)*inclination_function_f2mp[1][p as usize] *(associated_legendre_p_2m[2]-6.0*associated_legendre_p_2m[0]) *(
                                phase_21pq.sin()*find_rek2(frequ_21pq, nm_data, re_k2, w_lmpq)
                                -phase_21pq.cos()*find_imk2(frequ_21pq, nm_data, im_k2, w_lmpq)
                            )
                            -(1./6.)*inclination_function_f2mp[2][p as usize] *associated_legendre_p_2m[1] *( 
                                phase_22pq.cos()*find_rek2(frequ_22pq, nm_data, re_k2, w_lmpq)
                                +phase_22pq.sin()*find_imk2(frequ_22pq, nm_data, im_k2, w_lmpq)
                            )
                        );
                    }
                    sum_over_p = sum_over_p +sum_over_q;
                }
                force_normal = constante *sum_over_p;

            }
            particle.tides.parameters.internal.normal_component_of_the_tidal_force = force_normal;
            // panic!(" the normal force {} {}", force_normal, '\n');
            // particle.tides.parameters.internal.normal_component_of_the_tidal_force = 0.0; //NORMAL FORCE NULL
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

            let radial_component_of_the_tidal_force     = particle.tides.parameters.internal.radial_component_of_the_tidal_force_2_kaula;
            let orthogonal_component_of_the_tidal_force = particle.tides.parameters.internal.orthogonal_tidal_force_2_kaula;
            let normal_component_of_the_tidal_force     = particle.tides.parameters.internal.normal_component_of_the_tidal_force;
            
            let mut total_tidal_force_x:f64=0.;
            let mut total_tidal_force_y:f64=0.;
            let mut total_tidal_force_z:f64=0.;

            let force_2d:bool=true;
            let force_3d:bool=false;
            if force_2d {
                let distance:f64 = ( particle.tides.coordinates.position.x.powi(2) +particle.tides.coordinates.position.y.powi(2) +particle.tides.coordinates.position.z.powi(2) ).sqrt();
                total_tidal_force_x = radial_component_of_the_tidal_force *(particle.tides.coordinates.position.x / distance) - orthogonal_component_of_the_tidal_force *(particle.tides.coordinates.position.y / distance) ;
                total_tidal_force_y = radial_component_of_the_tidal_force *(particle.tides.coordinates.position.y / distance) + orthogonal_component_of_the_tidal_force *(particle.tides.coordinates.position.x / distance) ;
                total_tidal_force_z = 0.0;
            }
            if force_3d {
                // --- The cartesian angles --- //
                let distance:f64 = ( particle.tides.coordinates.position.x.powi(2) +particle.tides.coordinates.position.y.powi(2) +particle.tides.coordinates.position.z.powi(2) ).sqrt();
                let coplanar_distance:f64 = ( particle.tides.coordinates.position.x.powi(2) +particle.tides.coordinates.position.y.powi(2) ).sqrt();
                let cos_theta:f64 = particle.tides.coordinates.position.z / distance;
                let sin_theta:f64 = coplanar_distance / distance;
                let cos_phi:f64 = particle.tides.coordinates.position.x / coplanar_distance;
                let sin_phi:f64 = particle.tides.coordinates.position.y / coplanar_distance;
                
                total_tidal_force_x = radial_component_of_the_tidal_force
                                            *sin_theta *cos_phi 
                                        + normal_component_of_the_tidal_force
                                            *cos_theta *cos_phi 
                                        - orthogonal_component_of_the_tidal_force 
                                            *sin_phi ;
                total_tidal_force_y = radial_component_of_the_tidal_force
                                            *sin_theta *sin_phi
                                        + normal_component_of_the_tidal_force
                                            *cos_theta *sin_phi
                                        + orthogonal_component_of_the_tidal_force 
                                            *cos_phi ;
                total_tidal_force_z = radial_component_of_the_tidal_force
                                            *cos_theta
                                        - normal_component_of_the_tidal_force
                                            *sin_theta ;
            }

            sum_total_tidal_force.x += total_tidal_force_x;
            sum_total_tidal_force.y += total_tidal_force_y;
            sum_total_tidal_force.z += total_tidal_force_z;
            particle.tides.parameters.internal.tidal_force_x = total_tidal_force_x;
            particle.tides.parameters.internal.tidal_force_y = total_tidal_force_y;
            particle.tides.parameters.internal.tidal_force_z = total_tidal_force_z;

            // - Equation 19 from Bolmont et al. 2015 (first term)
            particle.tides.parameters.output.acceleration.x = factor1 * total_tidal_force_x; 
            particle.tides.parameters.output.acceleration.y = factor1 * total_tidal_force_y;
            particle.tides.parameters.output.acceleration.z = factor1 * total_tidal_force_z;
            
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

fn _calculate_gradient_factors(gm: f64, position: Axes, velocity: Axes) -> f64{
    let ( a, _perihelion_distance, e, i, longitude_perihelion, longitude_of_ascending_node, l, p, f
    ) = tools::calculate_keplerian_orbital_elements_2( gm, position, velocity);
    let varpi:f64 = longitude_perihelion -longitude_of_ascending_node ;
    let r:f64 = (position.x.powi(2) +position.y.powi(2) +position.z.powi(2)).sqrt();
    let coplanar_distance:f64 = ( position.x.powi(2) +position.y.powi(2) ).sqrt();
    let sin_theta:f64 = coplanar_distance / r;
    let cos_f:f64 = f.cos();
    let sin_f:f64 = f.sin();
    let sin_i = i.sin();
    let cos_i = i.cos();
    let e_2:f64 = e*e;
    let tmp:f64 = 1.0+ e*cos_f;
    let tmp_2:f64 = tmp*tmp;

    // let ha:f64 = 3.0 *(1.0-e_2)/tmp;
    // let he:f64 = -3.0 *a *(1.0+e_2*cos_f+2.0*e)/tmp_2;
    // let hf1:f64 = 3.0 *a*e*(1.0-e_2)*sin_f/tmp_2;
    // let hf2:f64 = -3.0*r *(-varpi+f).cos()*sin_i/sin_theta ;
    // let hf3:f64;
    // let hi1:f64;
    // let hi2:f64;
    // let hvarpi1:f64;
    // let hvarpi2:f64;
    return 0.;
}

fn _factorial(num: u64) -> u64 {
    match num {
        0 | 1 => 1,
        _ => _factorial(num - 1) * num,
    }
}

fn factorial(num :u64) -> f64{
    return _factorial(num) as f64;
}

fn delta_kron(num: u64) -> f64 {
    match num {
        0 | 0 => 1.0_f64,
        _ => 0_f64,
    }
}

fn phase_2mpq(m:f64, p:f64, q:f64, mean_anomaly:f64, theta:f64, argument_perihelion:f64, longitude_perihelion:f64 ) ->f64{
    return (2.-2.*p +q)*mean_anomaly -m*theta +(2.-2.*p)*argument_perihelion +m*longitude_perihelion;
}

fn eccentricty_function_g_20q(eccentricity: f64) -> [f64; 5]{

    let mut eccentricity_function = [0_f64;5];
    // in this notation:
    // [0] is -2
    // [1] is -1
    // [2] is 0
    // [3] is 1
    // [4] is 2

    let ecc:f64=eccentricity;
    let ecc_2:f64=ecc*ecc;
    let ecc_3:f64=ecc_2*ecc;
    let ecc_4:f64=ecc_3*ecc;
    let ecc_5:f64=ecc_4*ecc;
    let ecc_6:f64=ecc_5*ecc;
    let ecc_7:f64=ecc_6*ecc;

    eccentricity_function[0] = 0.;
    eccentricity_function[1] =  -0.5*ecc              +(1./16.)*ecc_3         -(5./384.)*ecc_5        -(143./18432.)*ecc_7 ;
    eccentricity_function[2] = 1.                    -(5./2.)*ecc_2          +(13./16.)*ecc_4        -(35./288.)*ecc_6 ;
    eccentricity_function[3] = (7./2.)*ecc           -(123./16.)*ecc_3       +(489./128.)*ecc_5      -(1763./2048.)*ecc_7 ;
    eccentricity_function[4] = (17./2.)*ecc_2        -(115./16.)*ecc_4       +(601./48.)*ecc_6 ;

    return eccentricity_function;
}
fn eccentricty_function_g_21q(eccentricity: f64) -> [f64; 5]{

    let mut eccentricity_function = [0_f64;5];
    // in this notation:
    // [0] is -2
    // [1] is -1
    // [2] is 0
    // [3] is 1
    // [4] is 2

    let ecc:f64=eccentricity;
    let ecc_2:f64=ecc*ecc;
    let ecc_3:f64=ecc_2*ecc;
    let ecc_4:f64=ecc_3*ecc;
    let ecc_5:f64=ecc_4*ecc;
    let ecc_6:f64=ecc_5*ecc;
    let ecc_7:f64=ecc_6*ecc;

    eccentricity_function[0] = (9./4.)*ecc_2         +(7./4.)*ecc_4          +(141./64.)*ecc_6 ;
    eccentricity_function[1] = (3./2.)*ecc           +(27./16.)*ecc_3        +(261./128.)*ecc_5      +(14309./6144.)*ecc_7 ;
    eccentricity_function[2] = (1.-ecc_2).powf(-3./2.) ;
    eccentricity_function[3] = eccentricity_function[1];
    eccentricity_function[4] = eccentricity_function[0];

    return eccentricity_function;
}

fn associated_legendre_p_2m(x: f64) -> [f64;3]{

    let mut p_2m = [0.0;3];
    let x_2 = x*x ;

    p_2m[0] = 0.5*(3.*x_2 -1.);
    p_2m[1] = -3.*x*((1.-x_2).sqrt()); //p_2m[1] = -3.*((1.-x_2).sqrt());
    p_2m[2] = 3.*(1.-x_2);

    return p_2m;
}

fn eccentricty_function_g(eccentricity: f64) -> [[f64; 15]; 3]{

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
    eccentricity_function[1][2] = (1773./256.)*ecc_5    -(4987./6144.)*ecc_7 ;
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

fn _derivative_of_eccentricity_function_g(eccentricity:f64) -> [[f64;15];3]{
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
    let mut dg_de = [[0.0;15];3];
    let ecc:f64=eccentricity;
    let ecc_2:f64=ecc*ecc;
    let ecc_3:f64=ecc_2*ecc;
    let ecc_4:f64=ecc_3*ecc;
    let ecc_5:f64=ecc_4*ecc;
    let ecc_6:f64=ecc_5*ecc;
    let ecc_7:f64=ecc_6*ecc;
    let ecc_8:f64=ecc_7*ecc;
    let ecc_9:f64=ecc_8*ecc;
    let ecc_11:f64=ecc_9*ecc_2;
    let ecc_13:f64=ecc_11*ecc_2;
 
    dg_de[0][0] = (3417968750./16647192600.)*ecc_13;
    dg_de[0][1] = (64./675.)*ecc_11;
    dg_de[0][2] = (6561./163840.)*ecc_9 +(19683./327680.)*ecc_11 +(45927./2097152.)*ecc_13;
    dg_de[0][3] = (1./72.)*ecc_7 +(7./288.)*ecc_9 +(49./4800.)*ecc_11;
    dg_de[0][4] = (1./384.)*ecc_5 +(11./2304.)*ecc_7 +(619./98304.)*ecc_9 +(3443./983040.)*ecc_11 +(685783./471859200.)*ecc_13;
    dg_de[0][5] = 0.;
    dg_de[0][6] = 0.5*ecc -(1./4.)*ecc_3 +(13./128.)*ecc_5 +(113./2304.)*ecc_7;
    dg_de[0][7] = -10.0*ecc +(63./2.)*ecc_3 -(155./6.)*ecc_5 -(2921./288.)*ecc_7;
    dg_de[0][8] = (49./2.)*ecc -(861./4.)*ecc_3 +(65925./128.)*ecc_5 -(132635./256.)*ecc_7;
    dg_de[0][9] = 289.*ecc_3 -(5865./8.)*ecc_5 +(203147./96.)*ecc_7;
    dg_de[0][10] = (714025./384.)*ecc_5 -(27483625./2304.)*ecc_7 +(2936126875./98304.)*ecc_9 -(6772518125./196608.)*ecc_11 +(303503554375./18874368.)*ecc_13;
    dg_de[0][11] = (284089./32.)*ecc_7 -(66328119./1280.)*ecc_8 +(191185929./2560.)*ecc_9;
    dg_de[0][12] = (52142352409./1474560.)*ecc_9 -(140254152605./589824.)*ecc_11 +(66020511589375./169869312.)*ecc_13;
    dg_de[0][13] = (5383010161./43200.)*ecc_11;
    dg_de[0][14] = (147483366698529./367001600.)*ecc_13;

    dg_de[1][0] = (1306918425967./471859200.)*ecc_13;
    dg_de[1][1] = (30089667./25600.)*ecc_11;
    dg_de[1][2] = (15717645./32768.)*ecc_9 -(8841951./65536.)*ecc_11 +(174091183./18874368.)*ecc_13;
    dg_de[1][3] = (5929./32.)*ecc_7 +(9933./128.)*ecc_9;
    dg_de[1][4] = (8427./128.)*ecc_5 +(20829./256.)*ecc_7 +(6019881./32768.)*ecc_9 +(29183787./327680.)*ecc_11  +(4288977063./52428800.)*ecc_13;
    dg_de[1][5] = (81./4.)*ecc_3 +(189./4.)*ecc_5 +(1661./16.)*ecc_7;
    dg_de[1][6] = (9./2.)*ecc +(81./4.)*ecc_3 +(6885./128.)*ecc_5 +(12123879./20608.)*ecc_7;
    dg_de[1][7] = 6.*ecc*(1.-ecc_2).powf(-4.);
    dg_de[1][8] = dg_de[1][6];
    dg_de[1][9] = dg_de[1][5];
    dg_de[1][10] = dg_de[1][4];
    dg_de[1][11] = dg_de[1][3];
    dg_de[1][12] = dg_de[1][2];
    dg_de[1][13] = dg_de[1][1];
    dg_de[1][14] = dg_de[1][0];

    dg_de[2][1] = dg_de[0][14];
    dg_de[2][2] = dg_de[0][13];
    dg_de[2][3] = dg_de[0][12];
    dg_de[2][4] = dg_de[0][11];
    dg_de[2][5] = dg_de[0][10];
    dg_de[2][6] = dg_de[0][9];
    dg_de[2][7] = dg_de[0][8];
    dg_de[2][8] = dg_de[0][7];
    dg_de[2][9] = dg_de[0][6];
    dg_de[2][10] = dg_de[0][5];
    dg_de[2][11] = dg_de[0][4];
    dg_de[2][12] = dg_de[0][3];
    dg_de[2][13] = dg_de[0][2];
    dg_de[2][14] = dg_de[0][1];
    dg_de[2][15] = dg_de[0][0];

    return dg_de;
}

fn inclination_function_f2mp(inclination: f64) -> [[f64; 3]; 3]{

    let mut inclination_function = [[0.;3];3];
    // in this notation: 0 is q=-2, 1 is q=-1, 2 is q=0, 3 is q=1, 4 is q=2 ... cf Kaula 64 Table 3
    let sin_incl = inclination.sin();
    let sin_incl_2 = sin_incl*sin_incl;
    let cos_incl = inclination.cos();

    inclination_function[0][0] = -(3./8.)*sin_incl_2;
    inclination_function[0][1] = (3./4.)*sin_incl_2 -(1./2.);
    inclination_function[0][2] = -(3./8.)*sin_incl_2;

    inclination_function[1][0] = (3./4.)*sin_incl*(1. +cos_incl);
    inclination_function[1][1] = -(3./2.)*sin_incl*cos_incl;
    inclination_function[1][2] = (3./4.)*sin_incl*(cos_incl -1.);

    inclination_function[2][0] = (3./4.)*(1. +cos_incl).powi(2);
    inclination_function[2][1] = (3./2.)*sin_incl.powi(2);
    inclination_function[2][2] = (3./4.)*(1. -cos_incl).powi(2);
    
    return inclination_function;
}

fn _derivative_of_inclination_function_f(inclination:f64) -> [[f64;3];3]{
    let mut df_di = [[0.;3];3];
    // in this notation: 0 is q=-2, 1 is q=-1, 2 is q=0, 3 is q=1, 4 is q=2 ... cf Kaula 64 Table 3
    let sin_incl = inclination.sin();
    let sin_incl_2 = sin_incl*sin_incl;
    let sin_incl_3 = sin_incl_2*sin_incl;
    let cos_incl = inclination.cos();
    let cos_incl_2 = cos_incl*cos_incl;

    df_di[0][0] = (9./16.)*sin_incl_3*cos_incl;
    df_di[0][1] = ((3./2.)*sin_incl_2-1.)*(3./2.)*cos_incl *sin_incl;
    df_di[0][2] = df_di[0][0];

    df_di[1][0] = (9./8.)*sin_incl*(1.+cos_incl)*(cos_incl*(1.+cos_incl) -sin_incl_2 );
    df_di[1][1] = (9./2.)*sin_incl*cos_incl*(cos_incl_2-sin_incl_2);
    df_di[1][2] = (9./8.)*sin_incl*cos_incl*(cos_incl-1.).powi(2)  -sin_incl_3*(cos_incl -1.);

    df_di[2][0] = -(9./4.)*sin_incl*(1.+cos_incl).powi(3);
    df_di[2][1] = 9.*sin_incl_3*cos_incl;
    df_di[2][2] = (9./4.)*sin_incl*(1.-cos_incl).powi(3);

    return df_di;
}

fn kaula_number(wk2:f64, nm_data:f64, real_part_love_number: [[f64;32];32], imaginary_part_love_number: [[f64;32];32], love_number_excitation_frequency: [[f64;32];32] ) -> (f64,f64){
    let mut w_k2 = wk2;
    let mut im_k2 = 0.;
    let mut re_k2 = 0.;
    let mut delta:f64 = 0.;
    let mut ctrl = true;
    // let mut find_it = false; 
    let mut parity = false;
    let mut case_oneline_last_col = false;
    let mut more_than_one_col = true;
    let nm_data = nm_data as usize;

    let last_line:usize;
    let last_column:usize;

    let mut line_up:usize=0;
    let mut column_up:usize=0;
    let mut line_down:usize=0;
    let mut column_down:usize=0;
    // println!("{} La frequence de ses morts {}{}", '\n','\t', w_k2 );
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
            im_k2 = imaginary_part_love_number[last_column][last_line]  
                + (w_k2 - love_number_excitation_frequency[last_column][last_line]) 
                * (imaginary_part_love_number[last_column][last_line] - imaginary_part_love_number[last_column-1][31])
                /(love_number_excitation_frequency[last_column][last_line]-love_number_excitation_frequency[last_column-1][31]);
            re_k2 = real_part_love_number[last_column][last_line]       
                + (w_k2 - love_number_excitation_frequency[last_column][last_line]) 
                * (real_part_love_number[last_column][last_line] - real_part_love_number[last_column-1][31])
                /(love_number_excitation_frequency[last_column][last_line]-love_number_excitation_frequency[last_column-1][31]) ;
        } else {
            delta = (w_k2 - love_number_excitation_frequency[last_column][last_line-1]) 
                / (love_number_excitation_frequency[last_column][last_line-1] - love_number_excitation_frequency[last_column][last_line-2]);
            im_k2 = delta
                *( imaginary_part_love_number[last_column][last_line-1] -imaginary_part_love_number[last_column][last_line-2] ) 
                +imaginary_part_love_number[last_column][last_line-1];
            re_k2 = delta
                *( real_part_love_number[last_column][last_line-1] -real_part_love_number[last_column][last_line-2] )
                +real_part_love_number[last_column][last_line-1];
            if im_k2<0. { // HAVE TO CHECK
                im_k2 = imaginary_part_love_number[last_column][last_line]
                    + (imaginary_part_love_number[last_column][last_line] - imaginary_part_love_number[last_column][last_line-1])/2.;
            }
            if re_k2<0. { // HAVE TO CHECK
                re_k2 = real_part_love_number[last_column][last_line] 
                    + (real_part_love_number[last_column][last_line] - real_part_love_number[last_column][last_line-1])/2.;
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
            im_k2 = (1.-delta)* imaginary_part_love_number[last_column-1][31]
                +delta*(imaginary_part_love_number[last_column][0]);
            re_k2 = (1.-delta)* real_part_love_number[last_column-1][31]
                +delta*(real_part_love_number[last_column][0]);
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
                        re_k2 = real_part_love_number[last_column][frequency2];
                        column_up  = last_column;
                        line_up    = frequency2;
                        column_down  = last_column;
                        line_down    = frequency2;
                        ctrl = false;
                        break;
                    } else if ctrl {
                        // delta = (frequency-W_lmpq(i-1)) / (W_lmpq(i)-W_lmpq(i-1))  realK2 = (1.-delta)*REK2(i-1) + delta*REK2(i)    imaginaryK2 = (1.-delta)*IMK2(i-1) + delta*IMK2(i)
                        delta = (w_k2 - love_number_excitation_frequency[last_column][frequency2-1]) / (love_number_excitation_frequency[last_column][frequency2] -love_number_excitation_frequency[last_column][frequency2-1]);
                        im_k2 = (1.-delta)* imaginary_part_love_number[last_column][frequency2-1]
                            +delta*(imaginary_part_love_number[last_column][frequency2]);
                        re_k2 = (1.-delta)* real_part_love_number[last_column][frequency2-1] 
                            +delta*(real_part_love_number[last_column][frequency2]);
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
                    im_k2 = (1.-delta)* imaginary_part_love_number[frequency1-1][31]
                        +delta*(imaginary_part_love_number[frequency1][0]);
                    re_k2 = (1.-delta)* real_part_love_number[frequency1-1][31]
                        +delta*(real_part_love_number[frequency1][0]);
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
                                re_k2 = real_part_love_number[frequency1][frequency2];
                                column_up  = frequency1;
                                line_up    = frequency2;
                                column_down  = frequency1;
                                line_down    = frequency2;
                                ctrl = false;
                                break;
                            }
                            else if ctrl { 
                                delta = (w_k2 - love_number_excitation_frequency[frequency1][frequency2-1]) / (love_number_excitation_frequency[frequency1][frequency2] -love_number_excitation_frequency[frequency1][frequency2-1]);
                                im_k2 = (1.-delta)* imaginary_part_love_number[frequency1][frequency2-1] 
                                    +delta*(imaginary_part_love_number[frequency1][frequency2]);
                                re_k2 = (1.-delta)* real_part_love_number[frequency1][frequency2-1]
                                    +delta*(real_part_love_number[frequency1][frequency2]);
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
        panic!("=== Error  in tide.rs :: fail interpolation Imk2 === found btw column&line {} {} < {} {} in frequency: {} < {} < {} \t the ImK2 {}", 
        column_down, line_down, column_up, line_up, love_number_excitation_frequency[column_down][line_down], w_k2, love_number_excitation_frequency[column_up][line_up], im_k2);
    }

    return (re_k2, im_k2)
}

fn sigma_2mpq(m:f64, p:f64, q:f64, spin:f64, orbital_frequency: f64) -> f64{
    // //////println!("\tsigma_calulation: the orbital_frequ {} and spin {}  give the sigma excit_frequ {}, the integer m {}, p {}, q {},  ", orbital_frequency, spin, (2. -2.*p + q)*orbital_frequency -m*spin, m, p, q );
    return (2.0_f64 -2.0_f64*p +q)*orbital_frequency -m*spin;
}

fn find_imk2(wk2:f64, nm_data:f64, imaginary_part_love_number: [[f64;32];32], love_number_excitation_frequency: [[f64;32];32] ) -> f64{
    let mut w_k2 = wk2;
    let mut im_k2 = 0.;
    let mut delta:f64 = 0.;
    let mut ctrl = true;
    // let mut find_it = false; 
    let mut parity = false;
    let mut case_oneline_last_col = false;
    let mut more_than_one_col = true;
    let nm_data = nm_data as usize;

    let last_line:usize;
    let last_column:usize;

    let mut line_up:usize=0;
    let mut column_up:usize=0;
    let mut line_down:usize=0;
    let mut column_down:usize=0;
    // println!("{} La frequence de ses morts {}{}", '\n','\t', w_k2 );
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
        panic!("=== Error  in tide.rs :: fail interpolation Imk2 === found btw column&line {} {} < {} {} in frequency: {} < {} < {} \t the ImK2 {}", 
        column_down, line_down, column_up, line_up, love_number_excitation_frequency[column_down][line_down], w_k2, love_number_excitation_frequency[column_up][line_up], im_k2);
    }
    return im_k2
}

fn find_rek2(wk2:f64, nm_data:f64, real_part_love_number: [[f64;32];32], love_number_excitation_frequency: [[f64;32];32] ) -> f64{
    let mut w_k2 = wk2;
    let mut re_k2 = 0.;
    let mut delta:f64 = 0.;
    let mut ctrl = true;
    // let mut find_it = false; 
    // let mut parity = false;
    let mut case_oneline_last_col = false;
    let mut more_than_one_col = true;
    let nm_data = nm_data as usize;

    let last_line:usize;
    let last_column:usize;

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
                        line_down    = frequency2-1;
                        column_down  = last_column;
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
    // println!("The last frequency {} ", love_number_excitation_frequency[last_column][last_line-1]);
    // println!("the before last {} ", love_number_excitation_frequency[last_column-1][31]);
    if (w_k2>love_number_excitation_frequency[column_up][line_up] && w_k2<love_number_excitation_frequency[last_column][last_line-1] ) || w_k2<love_number_excitation_frequency[column_down][line_down] {
        panic!("=== Error  in tide.rs :: fail interpolation Rek2 === found col {} line {} < {} {}  in frequency: {} < {} < {} \t ReK2 {}", 
        column_down, line_down, column_up, line_up, love_number_excitation_frequency[column_down][line_down], w_k2, love_number_excitation_frequency[column_up][line_up], re_k2);
    }

    return re_k2
}
