use std::collections::HashMap;
use super::super::tools;
use super::super::constants::{K2};
use super::super::{Particle};
use super::super::{Axes};
use super::{EvolutionType};

use crate::constants::G;
use crate::constants::PI;

//use crate::TidesTypes::{KaulaCoplanarModel, Nil};

// /////////////////////////////////////////////////////////////////////////////////////////////
// #[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
// //#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
// pub enum TidesTypes {
//     ConstTimeLagModel,
//     KaulaCoplanarModel,
//     NoTides,
// }

// #[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
// impl TidesTypes {
//     ConstTimeLagModel: Tides,
//     KaulaCoplanarModel : KaulaCoplanarTides,
// }

// #[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
// pub struct Machin {
//     pub lestypes: TidesTypes,
// }

// #[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
// //#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
// pub struct KaulaCoplanarTidesParticleParameters {
//     pub input: KaulaCoplanarTidesInputParticleParameters,
//     // pub internal: TidesParticleInternalParameters,
//     // pub output: TidesParticleOutputParameters,
// }

// #[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
// pub struct KaulaCoplanarTides {
//     pub parameters: KaulaCoplanarTidesParticleParameters,
// }

// impl KaulaCoplanarTides {
//     // pub fn new( love_number_excitation_frequency: Vec<f64>, real_part_love_number: Vec<f64>, imaginary_part_love_number: Vec<f64>) -> KaulaCoplanarTides {
//     pub fn new( love_number_excitation_frequency: [[f64;32];32], real_part_love_number: [[f64;32];32], imaginary_part_love_number: [[f64;32];32], num_datapoints:usize) -> KaulaCoplanarTides {
//         KaulaCoplanarTides {
//             parameters: KaulaCoplanarTidesParticleParameters {
//                 input: KaulaCoplanarTidesInputParticleParameters{
//                     love_number_excitation_frequency: love_number_excitation_frequency,
//                     real_part_love_number: real_part_love_number,
//                     imaginary_part_love_number: imaginary_part_love_number,
//                     num_datapoints: num_datapoints,
//                 },
//                 // internal: TidesParticleInternalParameters {
//                 //     distance: 0.,
//                 //     radial_velocity: 0.,
//                 //     scaled_dissipation_factor: dissipation_factor_scale*dissipation_factor,
//                 //     scalar_product_of_vector_position_with_stellar_spin: 0.,
//                 //     scalar_product_of_vector_position_with_planetary_spin: 0.,
//                 //     orthogonal_component_of_the_tidal_force_due_to_stellar_tide: 0.,
//                 //     orthogonal_component_of_the_tidal_force_due_to_planetary_tide: 0.,
//                 //     radial_component_of_the_tidal_force: 0.,
//                 //     radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass: 0.,
//                 //     denergy_dt: 0., // Only for history output
//                 //     lag_angle: 0., // It will be initialized the first time the evolver is called
//                 // },
//                 // output: TidesParticleOutputParameters {
//                 //     acceleration: Axes{x: 0., y: 0., z: 0.},
//                 //     dangular_momentum_dt: Axes{x: 0., y: 0., z: 0.},
//                 // },
//             },
//         }
//     }
// }
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
//#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TidesParticleInputParameters {
    pub dissipation_factor: f64,
    pub dissipation_factor_scale: f64, // to scale the dissipation factor (multiply)
    pub love_number: f64,   // Love number of degree 2 (i.e., k2). Dimensionless parameters that measure the rigidity of a planetary body and the 
                            // susceptibility of its shape to change in response to a tidal potential.
    pub kaula_coplanar_tides_input_parameters: KaulaCoplanarTidesInputParticleParameters,
    // pub love_number_excitation_frequency: [[f64;32];32],
    // pub real_part_love_number: [[f64;32];32],
    // pub imaginary_part_love_number: [[f64;32];32],
    // pub num_datapoints: i32,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
//#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct KaulaCoplanarTidesInputParticleParameters {
    pub love_number_excitation_frequency: [[f64;32];32],
    pub real_part_love_number: [[f64;32];32],
    pub imaginary_part_love_number: [[f64;32];32],
    pub num_datapoints: i32,
}
///////////////////////////////////////////////////////////////////////////////////////////////

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
    /////////////////////////////////////////////////////////////////
    //ConstTimeLagOrbitingBody,
    KaulaCoplanarOrbitingBody,
    /////////////////////////////////////////////////////////////////
    Disabled,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
//#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Tides {
    pub effect: TidesEffect,
    pub parameters: TidesParticleParameters,
    pub coordinates: TidesParticleCoordinates,
}

impl Tides {
    //pub fn new(effect: TidesEffect, dissipation_factor: f64, dissipation_factor_scale: f64, love_number: f64) -> Tides {
    pub fn new(effect: TidesEffect, dissipation_factor: f64, dissipation_factor_scale: f64, love_number: f64, love_number_excitation_frequency: [[f64;32];32], real_part_love_number: [[f64;32];32], imaginary_part_love_number: [[f64;32];32], num_datapoints: i32) -> Tides {
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
            //if let TidesEffect::ConstTimeLagOrbitingBody = particle.tides.effect {
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
            //if let TidesEffect::ConstTimeLagOrbitingBody = particle.tides.effect {
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
            //if let TidesEffect::ConstTimeLagOrbitingBody = particle.tides.effect {
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
                //println!("Insert {} in {}", planet_dependent_dissipation_factor, particle.id);
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
pub fn calculate_torque_due_to_tides(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], central_body:bool) {
    let mut dangular_momentum_dt = Axes{x: 0., y: 0., z:0.};
    let mut reference_spin = tidal_host_particle.spin.clone();
    let mut orthogonal_component_of_the_tidal_force: f64;
    let mut reference_rscalspin: f64;

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody = particle.tides.effect {
        //if let TidesEffect::ConstTimeLagOrbitingBody = particle.tides.effect {
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
            }
        }
    }

    if central_body {
        // - Equation 25 from Bolmont et al. 2015
        tidal_host_particle.tides.parameters.output.dangular_momentum_dt.x = dangular_momentum_dt.x;
        tidal_host_particle.tides.parameters.output.dangular_momentum_dt.y = dangular_momentum_dt.y;
        tidal_host_particle.tides.parameters.output.dangular_momentum_dt.z = dangular_momentum_dt.z;
    }
    //if mon stuff 
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody = particle.tides.effect {
        //if let TidesEffect::KaulaCoplanarOrbitingBody = particle.tides.effect {

            let orbital_elements = tools::calculate_keplerian_orbital_elements(G*(tidal_host_particle.mass+particle.mass), particle.heliocentric_position, particle.heliocentric_velocity);
            //it return (a, q, eccentricity, i, p, n, l, orbital_period)
            let semi_major_axis = orbital_elements.0;
            let eccentricity = orbital_elements.2;
            let orbital_period = orbital_elements.7;
            let orbital_frequency = (2.*PI) / (orbital_period *24.*60.*60.); // in 1/s
            let spin = particle.norm_spin_vector_2.sqrt();

            let mut tidal_torque_kaula_coplanar = [0.;3];

            let im_k2 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.imaginary_part_love_number;
            let re_k2 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.real_part_love_number;
            let w_lmpq = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.love_number_excitation_frequency;
            let nm_data = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.num_datapoints;

            //let mut e = 0.1;

            let eccentricity_function = eccentricty_function_g(eccentricity);

            let mut sum_g_im_k2_over_q = 0.;
            let mut q = -2.;

            for x in 0..5{
                let excitative_frequency = 2.*(spin - orbital_frequency) + (q)*spin ;
                let imaginary_kaula_number = kaula_number(excitative_frequency, nm_data, re_k2, im_k2, w_lmpq);
                sum_g_im_k2_over_q = sum_g_im_k2_over_q + ( eccentricity_function[0][x].powf(2.) * imaginary_kaula_number.1 );
                q+=1.;
            }

            //tidal_torque_Kaula_coplanar[2] = (3./2.)*(G* tidal_host_particle.mass.powf(2.)* particle.radius.powf(5.))*( eccentricity_function[0][0] +eccentricity_function[0][1] +eccentricity_function[0][2] +eccentricity_function[0][3] +eccentricity_function[0][4]);
            tidal_torque_kaula_coplanar[2] = (3./2.)*((G* tidal_host_particle.mass.powf(2.)* particle.radius.powf(5.))/semi_major_axis.powf(6.)) *sum_g_im_k2_over_q;
        }
    }
}

pub fn calculate_orthogonal_component_of_the_tidal_force(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>) {
    let mut tidal_host_particle = tidal_host_particle;
    let mut particles = particles;
    let mut more_particles = more_particles;
    let mut star_planet_dependent_dissipation_factors = star_planet_dependent_dissipation_factors;
    let central_body = true;
    calculate_orthogonal_component_of_the_tidal_force_for(central_body, &mut tidal_host_particle, &mut particles, &mut more_particles, &mut star_planet_dependent_dissipation_factors);
    calculate_orthogonal_component_of_the_tidal_force_for(!central_body, &mut tidal_host_particle, &mut particles, &mut more_particles, &mut star_planet_dependent_dissipation_factors);
}

fn calculate_orthogonal_component_of_the_tidal_force_for(central_body:bool, tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody = particle.tides.effect {
        //if let TidesEffect::ConstTimeLagOrbitingBody = particle.tides.effect {
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

                // SBC
                //println!("> {:e} {:e} {:e} {:e}", tidal_host_particle.mass_g, particle.radius.powi(10), particle.scaled_dissipation_factor, distance_7);
                //println!("> {:e} {:e} {:e}", particle.tides.coordinates.position.x, particle.tides.coordinates.position.y, particle.tides.coordinates.position.z);
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
            if !central_body{
                let orbital_elements = tools::calculate_keplerian_orbital_elements(G*(tidal_host_particle.mass+particle.mass), particle.heliocentric_position, particle.heliocentric_velocity);
                //it return (a, q, eccentricity, i, p, n, l, orbital_period)
                let semi_major_axis = orbital_elements.0;
                let eccentricity = orbital_elements.2;
                let orbital_period = orbital_elements.7;
                let orbital_frequency = (2.*PI) / (orbital_period *24.*60.*60.); // in 1/s
                let spin = particle.norm_spin_vector_2.sqrt();

                let imaginary_k2 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.imaginary_part_love_number;
                let real_k2 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.real_part_love_number;
                let w_lmpq = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.love_number_excitation_frequency;
                let nm_data = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.num_datapoints;

                let eccentricity_function = eccentricty_function_g(eccentricity);
                let mut excitative_frequency: f64;
                let mut sum_over_q = (3./2.)*(G* tidal_host_particle.mass.powf(2.)* particle.radius.powf(5.)) / (semi_major_axis.powf(6.) * particle.heliocentric_distance);
                let mut q = -2.;

                for x in 0..5{
                
                    // sum_over_q = sum_over_q -(3./4.) ( eccentricity_function[1][x].powf(2.) * real_part_love_number(nm_data, re_k2, w_lmpq, q, orbital_frequency, spin) );
                    excitative_frequency = sigma_2mpq(2., 0., q, spin, orbital_frequency);
                    let im_kaula_number = kaula_number(excitative_frequency, nm_data, real_k2, imaginary_k2, w_lmpq);
                    sum_over_q = sum_over_q -( eccentricity_function[0][x].powf(2.) * im_kaula_number.1 );
                    
                    q+=1.;
                }
                particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = sum_over_q

            }
        }
    }
}

pub fn calculate_radial_component_of_the_tidal_force(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>) {
    let star_mass_2 = tidal_host_particle.mass * tidal_host_particle.mass;

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody = particle.tides.effect {
        //if let TidesEffect::ConstTimeLagOrbitingBody = particle.tides.effect {
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
        }

        if let TidesEffect::KaulaCoplanarOrbitingBody = particle.tides.effect {
            
            let orbital_elements = tools::calculate_keplerian_orbital_elements(G*(tidal_host_particle.mass+particle.mass), particle.heliocentric_position, particle.heliocentric_velocity);
            //it return (a, q, eccentricity, i, p, n, l, orbital_period)
            let semi_major_axis = orbital_elements.0;
            let eccentricity = orbital_elements.2;
            let orbital_period = orbital_elements.7;
            let orbital_frequency = (2.*PI) / (orbital_period *24.*60.*60.); // in 1/s
            let spin = particle.norm_spin_vector_2.sqrt();

            let imaginary_k2 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.imaginary_part_love_number;
            let real_k2 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.real_part_love_number;
            let w_lmpq = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.love_number_excitation_frequency;
            let nm_data = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.num_datapoints;

            let eccentricity_function = eccentricty_function_g(eccentricity);
            let mut excitative_frequency: f64;
            let mut sum_over_q = 0.;
            let mut q = -2.;
            let mut _j = -2.;

            for x in 0..5{
                
                // sum_over_q = sum_over_q -(3./4.) ( eccentricity_function[1][x].powf(2.) * real_part_love_number(nm_data, re_k2, w_lmpq, q, orbital_frequency, spin) );
                excitative_frequency = sigma_2mpq(0., 1., q, spin, orbital_frequency);
                let mut real_kaula_number = kaula_number(excitative_frequency, nm_data, real_k2, imaginary_k2, w_lmpq);
                sum_over_q = sum_over_q -(3./4.)*( eccentricity_function[1][x].powf(2.) * real_kaula_number.0 );
                
                excitative_frequency = sigma_2mpq(2., 0., q, spin, orbital_frequency);
                real_kaula_number = kaula_number( excitative_frequency, nm_data, real_k2, imaginary_k2, w_lmpq );
                sum_over_q = sum_over_q -(9./4.)*( eccentricity_function[0][x].powf(2.) * real_kaula_number.0 );
                
                q+=1.;
            }
            particle.tides.parameters.internal.radial_component_of_the_tidal_force = (G* tidal_host_particle.mass.powf(2.)* particle.radius.powf(5.)) / (semi_major_axis.powf(6.) * particle.heliocentric_distance);
            
        }
    }
}

pub fn calculate_denergy_dt(particles: &mut [Particle], more_particles: &mut [Particle]) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody = particle.tides.effect {
        //if let TidesEffect::ConstTimeLagOrbitingBody = particle.tides.effect {
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


pub fn calculate_tidal_acceleration(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
    let factor2 = 1. / tidal_host_particle.mass;
    let mut sum_total_tidal_force = Axes{x:0., y:0., z:0.};

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        //*****************Here implement the foce due to tidal FORCE**********************
        //if CTL
        //elif CREEP
        //elif KAULA
        if let TidesEffect::OrbitingBody = particle.tides.effect {
        //if let TidesEffect::ConstTimeLagOrbitingBody = particle.tides.effect {
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
            //println!("factor3 {:e} {:e} {:e}", particle.tides.parameters.internal.radial_component_of_the_tidal_force, particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide, particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide);
            //println!("d {:e} vrad {:e}", particle.tides.parameters.internal., particle.tides.parameters.internal.radial_velocity);
            //println!("total {:e} {:e} {:e}", total_tidal_force_x, total_tidal_force_y, total_tidal_force_z);

            sum_total_tidal_force.x += total_tidal_force_x;
            sum_total_tidal_force.y += total_tidal_force_y;
            sum_total_tidal_force.z += total_tidal_force_z;

            // - Equation 19 from Bolmont et al. 2015 (first term)
            particle.tides.parameters.output.acceleration.x = factor1 * total_tidal_force_x; 
            particle.tides.parameters.output.acceleration.y = factor1 * total_tidal_force_y;
            particle.tides.parameters.output.acceleration.z = factor1 * total_tidal_force_z;
        }
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
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Let's star the modifications:

pub fn eccentricty_function_g(eccentricity: f64) -> [[f64; 5]; 3]{

    let mut eccentricity_function = [[0.9;5];3];
    // in this notation: 0 is q=-2, 1 is q=-1, 2 is q=0, 3 is q=1, 4 is q=2 ... cf Kaula 64 Table 3

    eccentricity_function[0][0] = 0.;
    eccentricity_function[0][1] = -(1./2.)*eccentricity +(1./16.)*eccentricity.powf(3.);
    eccentricity_function[0][2] = 1. -(5./2.)*eccentricity.powf(2.) +(13./16.)*eccentricity.powf(4.);
    eccentricity_function[0][3] = (7./2.)*eccentricity -(123./16.)*eccentricity.powf(3.);
    eccentricity_function[0][4] = (17./2.)*eccentricity.powf(2.) -(115./6.)*eccentricity.powf(4.);

    eccentricity_function[1][0] = (9./4.)*eccentricity.powf(2.) +(7./4.)*eccentricity.powf(4.);
    eccentricity_function[1][1] = (3./2.)*eccentricity +(27./16.)*eccentricity.powf(3.);
    eccentricity_function[1][2] = (1. -eccentricity.powf(2.)).powf(-(3./2.));
    eccentricity_function[1][3] = (3./2.)*eccentricity +(27./16.)*eccentricity.powf(3.);
    eccentricity_function[1][4] = (9./4.)*eccentricity.powf(2.) +(7./4.)*eccentricity.powf(4.);

    eccentricity_function[2][0] = (17./2.)*eccentricity.powf(2.) -(115./6.)*eccentricity.powf(4.);
    eccentricity_function[2][1] = (7./2.)*eccentricity -(123./16.)*eccentricity.powf(3.);
    eccentricity_function[2][2] = 1. -(5./2.)*eccentricity.powf(2.) +(13./16.)*eccentricity.powf(4.);
    eccentricity_function[2][3] = -(1./2.)*eccentricity + (1./16.)*eccentricity.powf(3.);
    eccentricity_function[2][4] = 0.;

    return eccentricity_function;
}

pub fn _inclination_function_f(inclination: f64) -> [[f64; 3]; 3]{

    let mut inclination_function = [[0.;3];3];
    // in this notation: 0 is q=-2, 1 is q=-1, 2 is q=0, 3 is q=1, 4 is q=2 ... cf Kaula 64 Table 3

    inclination_function[0][0] = -(3./8.)*inclination.sin().powf(2.);
    inclination_function[0][1] = (3./4.)*inclination.sin().powf(2.) -(1./2.);
    inclination_function[0][2] = -(3./8.)*inclination.sin().powf(2.);

    inclination_function[1][0] = (3./4.)*inclination.sin()*(1. +inclination.cos());
    inclination_function[1][1] = -(3./2.)*inclination.sin()*inclination.cos();
    inclination_function[1][2] = (3./4.)*inclination.sin()*(1. +inclination.cos());

    inclination_function[2][0] = (3./4.)*(1. +inclination.cos()).powf(2.);
    inclination_function[2][1] = (3./2.)*inclination.sin();
    inclination_function[2][2] = (3./4.)*(1. -inclination.cos());
    
    return inclination_function;
}

// pub fn imaginary_kaula_number(nm_data:i32, imaginary_part_love_number: [[f64;32];32], love_number_excitation_frequency: [[f64;32];32], integer_q: f64, orbital_frequency: f64, spin_frequency: f64) -> f64{

//     let mut w_k2 = 2.*(spin_frequency - orbital_frequency) + integer_q*spin_frequency ;
//     let mut im_k2 = 0.;
//     let mut love_frequency_low_value = love_number_excitation_frequency[0][0];
//     let mut love_number_low_value = imaginary_part_love_number[0][0];
//     let mut x = 0;

//     for frequency1 in 0..love_number_excitation_frequency.len(){
//         for frequency2 in 0..love_number_excitation_frequency.len(){
//             if love_number_excitation_frequency[frequency1][frequency2] > w_k2 {
//                 if love_number_excitation_frequency[frequency1][frequency2] == love_frequency_low_value {
//                     w_k2 = love_number_excitation_frequency[frequency1][frequency2] ;
//                     im_k2 = imaginary_part_love_number[frequency1][frequency2];
//                     break;
//                 }
//                 else { 
//                     w_k2 = ( love_number_excitation_frequency[frequency1][frequency2] - love_frequency_low_value )/ 2.;
//                     im_k2 = ( imaginary_part_love_number[frequency1][frequency2] - love_number_low_value )/2.;
//                     break;
//                 }
//             }
//             love_frequency_low_value = love_number_excitation_frequency[frequency1][frequency2];
//             love_number_low_value = imaginary_part_love_number[frequency1][frequency2];

//             x = x +1;
//             if x==nm_data{break;}

//             // match w_k2.cmp(&love_number_excitation_frequency[frequency]){
//             //     Ordering::Greater =>,
//             //     Ordering::Equal => ,
//             //     Ordering::Less =>,
//             // }
//         }
//     }
//     return im_k2
// }
// //Thoses two function must be just one!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// pub fn real_kaula_number( w_k2:f64, nm_data:i32, real_part_love_number: [[f64;32];32], love_number_excitation_frequency: [[f64;32];32], _integer_q: f64, _orbital_frequency: f64, _spin_frequency: f64) -> f64{

//     //let mut w_k2 = 2.*(spin_frequency - orbital_frequency) + integer_q*spin_frequency ;
//     let mut re_k2 = 0.;
//     let mut love_frequency_low_value = love_number_excitation_frequency[0][0];
//     let mut love_number_low_value = real_part_love_number[0][0];
//     let mut x = 0;

//     for frequency1 in 0..love_number_excitation_frequency.len(){
//         for frequency2 in 0..love_number_excitation_frequency.len(){
//             if love_number_excitation_frequency[frequency1][frequency2] > w_k2 {
//                 if love_number_excitation_frequency[frequency1][frequency2] == love_frequency_low_value {
//                     //w_k2 = love_number_excitation_frequency[frequency1][frequency2] ;
//                     re_k2 = real_part_love_number[frequency1][frequency2];
//                     break;
//                 }
//                 else { 
//                     //w_k2 = ( love_number_excitation_frequency[frequency1][frequency2] - love_frequency_low_value )/ 2.;
//                     re_k2 = ( real_part_love_number[frequency1][frequency2] - love_number_low_value )/2.;
//                     break;
//                 }
//             }
//             love_frequency_low_value = love_number_excitation_frequency[frequency1][frequency2];
//             love_number_low_value = real_part_love_number[frequency1][frequency2];

//             x = x +1;
//             if x==nm_data{break;}

//         }
//     }
//     return re_k2
// }

pub fn sigma_2mpq(m:f64, p:f64, q:f64, spin:f64, orbital_frequency: f64) -> f64{

    return (2. -2.*p + q)*spin -m*orbital_frequency;
}

pub fn kaula_number(w_k2:f64, nm_data:i32, real_part_love_number: [[f64;32];32], imaginary_part_love_number: [[f64;32];32], love_number_excitation_frequency: [[f64;32];32] ) -> (f64,f64){

    //let mut w_k2 = 2.*(spin_frequency - orbital_frequency) + integer_q*spin_frequency ;
    let mut re_k2 = 0.;
    let mut im_k2 = 0.;
    let mut love_frequency_low_value = love_number_excitation_frequency[0][0];
    let mut im_love_number_low_value = imaginary_part_love_number[0][0];
    let mut re_love_number_low_value = real_part_love_number[0][0];
    let mut x = 0;

    for frequency1 in 0..love_number_excitation_frequency.len(){
        for frequency2 in 0..love_number_excitation_frequency.len(){
            if love_number_excitation_frequency[frequency1][frequency2] > w_k2 {
                if love_number_excitation_frequency[frequency1][frequency2] == love_frequency_low_value {
                    //w_k2 = love_number_excitation_frequency[frequency1][frequency2];
                    re_k2 = real_part_love_number[frequency1][frequency2];
                    im_k2 = imaginary_part_love_number[frequency1][frequency2];
                    break;
                }
                else { 
                    //w_k2 = ( love_number_excitation_frequency[frequency1][frequency2] - love_frequency_low_value )/ 2.;
                    re_k2 = ( real_part_love_number[frequency1][frequency2] - re_love_number_low_value )/2.;
                    im_k2 = ( imaginary_part_love_number[frequency1][frequency2] - im_love_number_low_value )/2.;
                    break;
                }
            }
            love_frequency_low_value = love_number_excitation_frequency[frequency1][frequency2];
            re_love_number_low_value = real_part_love_number[frequency1][frequency2];
            im_love_number_low_value = imaginary_part_love_number[frequency1][frequency2];

            x = x +1;
            if x==nm_data{break;}

            // match w_k2.cmp(&love_number_excitation_frequency[frequency]){
            //     Ordering::Greater =>,
            //     Ordering::Equal =>,
            //     Ordering::Less =>,
            // }
        }
    }
    return (re_k2, im_k2)
}