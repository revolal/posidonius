use std::collections::HashMap;
use super::super::tools;
use super::super::constants::{K2};
use super::super::{Particle};
use super::super::{Axes};
use super::{EvolutionType};


use crate::constants::{G,TWO_PI, DAY, AU};//, G_SI, AU, M_SUN};

// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////MY MODIFICATION
// #[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
// pub enum Tides {
//     // CTLCentralBody(dissipation_factor: f64, dissipation_factor_scale: f64, love_number: f64),
//     // CTLOrbitingBody(dissipation_factor: f64, dissipation_factor_scale: f64, love_number: f64),
//     // KCOrbitingBody(love_number_excitation_frequency: f64, real_part_love_number: f64, imaginary_part_love_number: f64, num_datapoints: f64),
//     CTLCentralBody(f64, f64, f64),
//     CTLOrbitingBody( f64, f64, f64),
//     KaulaCoplanarOrbitingBody(f64, f64, f64, f64),
//     Disabled,
// }
// #[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
// pub struct CTLCentralBody{
//     pub effect: TidesEffect,
//     pub parameters: TidesParticleParameters,
//     pub coordinates: TidesParticleCoordinates,
// }
// #[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
// pub struct CTLOrbitingBody{
//     pub effect: TidesEffect,
//     pub parameters: TidesParticleParameters,
//     pub coordinates: TidesParticleCoordinates,
// }
// #[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
// pub struct KCOrbitingBody{
//     pub effect: TidesEffect,
//     pub parameters: KaulaCoplanarTidesInputParticleParameters,
//     pub coordinates: TidesParticleCoordinates,
// }

// ////
// #[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
// pub enum TidesEffect {
//     CentralBody,
//     OrbitingBody,
//     KaulaCoplanarOrbitingBody,
//     Disabled,
// }
// ///////
// // #[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
// // pub struct KaulaCoplanarTidesParticleParameters {
// //     pub input: KaulaCoplanarTidesInputParticleParameters,
// //     pub output: TidesParticleOutputParameters,
// // }
// #[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
// pub struct KaulaCoplanarTidesInputParticleParameters{
//     pub love_number_excitation_frequency: [[f64;32];32],
//     pub real_part_love_number: [[f64;32];32],
//     pub imaginary_part_love_number: [[f64;32];32],
//     // pub love_number_excitation_frequency: f64,
//     // pub real_part_love_number: f64,
//     // pub imaginary_part_love_number: f64,
//     pub num_datapoints: f64,
// }
// //
// #[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
// pub struct TidesParticleParameters {
//     pub input: TidesParticleInputParameters,
//     pub internal: TidesParticleInternalParameters,
//     pub output: TidesParticleOutputParameters,
// }
// #[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
// pub struct TidesParticleInputParameters {
//     pub dissipation_factor: f64,
//     pub dissipation_factor_scale: f64, // to scale the dissipation factor (multiply)
//     pub love_number: f64,   // Love number of degree 2 (i.e., k2). Dimensionless parameters that measure the rigidity of a planetary body and the 
//                             // susceptibility of its shape to change in response to a tidal potential.
// }
// #[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
// pub struct TidesParticleInternalParameters {
//     pub distance: f64,
//     pub radial_velocity: f64,
//     pub scaled_dissipation_factor: f64, // sigma (dissipation_factor_scale*dissipation_factor)
//     pub scalar_product_of_vector_position_with_stellar_spin: f64,
//     pub scalar_product_of_vector_position_with_planetary_spin: f64,
//     pub orthogonal_component_of_the_tidal_force_due_to_stellar_tide: f64,
//     pub orthogonal_component_of_the_tidal_force_due_to_planetary_tide: f64,
//     pub radial_component_of_the_tidal_force: f64,
//     pub radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass: f64, // Needed to compute denergy_dt
//     pub denergy_dt: f64, // Only for history output
//     pub lag_angle: f64, // Used by EvolutionType::BolmontMathis2016, EvolutionType::GalletBolmont2017 and EvolutionType::LeconteChabrier2013(true)
// }
// #[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
// pub struct TidesParticleOutputParameters {
//     pub acceleration: Axes,
//     pub dangular_momentum_dt: Axes, // Force
// }
// //
// #[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
// pub struct TidesParticleCoordinates {
//     // Positions/velocities in a heliocentric frame 
//     // (i.e., the host is at rest with respect to the origin of the coordinate system)
//     pub position: Axes,
//     pub velocity: Axes,
// }


// impl CTLCentralBody {
//     pub fn new(dissipation_factor: f64, dissipation_factor_scale: f64, love_number: f64) -> CTLCentralBody{
//         CTLCentralBody{
//             parameters: TidesParticleParameters {
//                 input: TidesParticleInputParameters {
//                     dissipation_factor: dissipation_factor,
//                     dissipation_factor_scale: dissipation_factor_scale,
//                     love_number: love_number,
//                 },
//                 internal: TidesParticleInternalParameters {
//                     distance: 0.,
//                     radial_velocity: 0.,
//                     scaled_dissipation_factor: dissipation_factor_scale*dissipation_factor,
//                     scalar_product_of_vector_position_with_stellar_spin: 0.,
//                     scalar_product_of_vector_position_with_planetary_spin: 0.,
//                     orthogonal_component_of_the_tidal_force_due_to_stellar_tide: 0.,
//                     orthogonal_component_of_the_tidal_force_due_to_planetary_tide: 0.,
//                     radial_component_of_the_tidal_force: 0.,
//                     radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass: 0.,
//                     denergy_dt: 0., // Only for history output
//                     lag_angle: 0., // It will be initialized the first time the evolver is called
//                 },
//                 output: TidesParticleOutputParameters {
//                     acceleration: Axes{x: 0., y: 0., z: 0.},
//                     dangular_momentum_dt: Axes{x: 0., y: 0., z: 0.},
//                 },
//             },
//             coordinates: TidesParticleCoordinates {
//                 position: Axes{x: 0., y: 0., z: 0.},
//                 velocity: Axes{x: 0., y: 0., z: 0.},
//             },
//         }
//     }
// }
// impl CTLOrbitingBody {
//     pub fn new(dissipation_factor: f64, dissipation_factor_scale: f64, love_number: f64) -> CTLCentralBody{
//         CTLCentralBody{
//             parameters: TidesParticleParameters {
//                 input: TidesParticleInputParameters {
//                     dissipation_factor: dissipation_factor,
//                     dissipation_factor_scale: dissipation_factor_scale,
//                     love_number: love_number,
//                 },
//                 internal: TidesParticleInternalParameters {
//                     distance: 0.,
//                     radial_velocity: 0.,
//                     scaled_dissipation_factor: dissipation_factor_scale*dissipation_factor,
//                     scalar_product_of_vector_position_with_stellar_spin: 0.,
//                     scalar_product_of_vector_position_with_planetary_spin: 0.,
//                     orthogonal_component_of_the_tidal_force_due_to_stellar_tide: 0.,
//                     orthogonal_component_of_the_tidal_force_due_to_planetary_tide: 0.,
//                     radial_component_of_the_tidal_force: 0.,
//                     radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass: 0.,
//                     denergy_dt: 0., // Only for history output
//                     lag_angle: 0., // It will be initialized the first time the evolver is called
//                 },
//                 output: TidesParticleOutputParameters {
//                     acceleration: Axes{x: 0., y: 0., z: 0.},
//                     dangular_momentum_dt: Axes{x: 0., y: 0., z: 0.},
//                 },
//             },
//             coordinates: TidesParticleCoordinates {
//                 position: Axes{x: 0., y: 0., z: 0.},
//                 velocity: Axes{x: 0., y: 0., z: 0.},
//             },
//         }
//     }
// }
// impl KCOrbitingBody {
//     pub fn new(love_number_excitation_frequency: [[f64;32];32], real_part_love_number: [[f64;32];32], imaginary_part_love_number: [[f64;32];32], num_datapoints: i32) -> KCOrbitingBody{
//         KCOrbitingBody{
//             parameters: KaulaCoplanarTidesInputParticleParameters {
//                 love_number_excitation_frequency: love_number_excitation_frequency,
//                 real_part_love_number: real_part_love_number,
//                 imaginary_part_love_number: imaginary_part_love_number,
//                 num_datapoints: num_datapoints,
//             },
//             coordinates: TidesParticleCoordinates {
//                 position: Axes{x: 0., y: 0., z: 0.},
//                 velocity: Axes{x: 0., y: 0., z: 0.},
//             },
//         }
//     }
// }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

    // pub im_love_number_sigma220_2: f64,
    // pub im_love_number_sigma220_1: f64,
    pub im_love_number_sigma2200: f64,
    // pub im_love_number_sigma2201: f64,
    // pub im_love_number_sigma2202: f64,

    // pub re_love_number_sigma220_2: f64,
    // pub re_love_number_sigma220_1: f64,
    // pub re_love_number_sigma2200: f64,
    // pub re_love_number_sigma2201: f64,
    // pub re_love_number_sigma2202: f64,

    // pub im_love_number:f64,
    // pub re_love_number:f64,

    // pub sigma220_2_excitative_frequency: f64,
    // pub sigma220_1_excitative_frequency: f64,
    pub sigma2200_excitative_frequency: f64,
    // pub sigma2201_excitative_frequency: f64,
    // pub sigma2202_excitative_frequency: f64,

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

                    // sigma220_2_excitative_frequency:0.,
                    // sigma220_1_excitative_frequency:0.,
                    sigma2200_excitative_frequency:0.,
                    // sigma2201_excitative_frequency:0.,
                    // sigma2202_excitative_frequency:0.,

                    // im_love_number_sigma220_2: 0.,
                    // im_love_number_sigma220_1: 0.,
                    im_love_number_sigma2200: 0.,
                    // im_love_number_sigma2201: 0.,
                    // im_love_number_sigma2202: 0.,
                
                    // re_love_number_sigma220_2: 0.,
                    // re_love_number_sigma220_1: 0.,
                    // re_love_number_sigma2200: 0.,
                    // re_love_number_sigma2201: 0.,
                    // re_love_number_sigma2202: 0.,

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
pub fn calculate_torque_due_to_tides(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle], central_body:bool, _current_time:f64) {
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
        //if mon stuff 
        
        if let TidesEffect::KaulaCoplanarOrbitingBody = particle.tides.effect {

            // let orbital_elements = tools::calculate_keplerian_orbital_elements(G*(tidal_host_particle.mass), particle.inertial_position, particle.inertial_velocity);
            let orbital_elements = tools::calculate_keplerian_orbital_elements(G*(tidal_host_particle.mass+particle.mass), particle.heliocentric_position, particle.heliocentric_velocity);
            //it return (a, q, eccentricity, i, p, n, l, orbital_period)
            let semi_major_axis = orbital_elements.0;
            let eccentricity = orbital_elements.2;

            let orbital_period = orbital_elements.7; //in [s]
            let orbital_frequency = TWO_PI / (orbital_period*DAY) ; // in [rad/s]
            //////println!("a {}\t period [s] {}\t frequency [s^-1] {}\t ", semi_major_axis, orbital_period, orbital_frequency);
            // let orbital_period = TWO_PI *( (semi_major_axis *AU).powi(3)/(G_SI*M_SUN*(tidal_host_particle.mass+particle.mass))).sqrt();
            // let orbital_frequency = TWO_PI / (orbital_period) ; // in [rad/s]
            // //////println!("a {}\t period [s] {}\t frequency [s^-1] {}\t ", semi_major_axis, orbital_period, orbital_frequency);

            // let spin = particle.norm_spin_vector_2.sqrt()/DAY; //this is [rad.s^-1]
            let spin = ( particle.spin.x.powi(2) + particle.spin.y.powi(2) + particle.spin.z.powi(2) ).sqrt() / DAY; //this is [rad.s^-1]
            //////println!("The Planetary Spin in [rad.s^-1]: {}\t Period in [s.rad^-1]: {}\t in [day]: {} Vector: {:?}", spin, 1.0/(spin), TWO_PI/(spin*DAY), particle.spin);
            let mut tidal_torque_kaula_coplanar = [0.;3];

            let im_k2 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.imaginary_part_love_number;
            let re_k2 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.real_part_love_number;
            let w_lmpq = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.love_number_excitation_frequency;
            let nm_data = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.num_datapoints;

            let eccentricity_function = eccentricty_function_g_20q(eccentricity);
            let mut _excitative_frequency: f64;
            // let mut sum_g_im_k2_over_q = 0.;
            let mut _q = -2.;

            particle.tides.parameters.internal.a = semi_major_axis;
            particle.tides.parameters.internal.spin = spin;
            particle.tides.parameters.internal.orbital_frequency = orbital_frequency;

            //////println!("\nTroque due to tides, KaulaCoplanarOrbitingBody:");
            // //////println!("\n\nWELLL {}\n", orbital_frequency);

            // if current_time > 10000.*orbital_elements.7 {
            //     panic!("Seem it's working.");
            // }

            let mut do_calculation:bool=false;
            if (1.- particle.tides.parameters.internal.check_excitative_frequ / sigma_2mpq(2., 0., -2., spin, orbital_frequency)).abs() > 0.001 {
                // //println!("Excitative frequency Variation {}", (1.- particle.tides.parameters.internal.check_excitative_frequ / sigma_2mpq(2., 0., -2., spin, orbital_frequency)).abs() );
                particle.tides.parameters.internal.check_excitative_frequ = sigma_2mpq(2., 0., -2., spin, orbital_frequency);
                do_calculation = true;
            }
            // do_calculation = true;
            if do_calculation {
                
                // //println!("\n\tDO_CALCULATION");
                // //println!("\nsigma220_2");
                particle.tides.parameters.internal.excitative_frequ_sigma220q[0] = sigma_2mpq(2., 0., -2., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma220q[0], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma220q[0] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma220q[0] = kaula.0 ;
                // particle.tides.parameters.internal.im_love_number_sigma220_2 = kaula.1;
                // //println!("\nsigma201_2");
                particle.tides.parameters.internal.excitative_frequ_sigma201q[0] = sigma_2mpq(0., 1., -2., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma201q[0], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma201q[0] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma201q[0] = kaula.0 ;

                // //println!("sigma220_2:{} \t Im:{} \t Re:{}", particle.tides.parameters.internal.sigma220_2_excitative_frequency, particle.tides.parameters.internal.im_love_number_sigma220_2, particle.tides.parameters.internal.re_love_number_sigma220_2);

                // //println!("sigma220_1");
                particle.tides.parameters.internal.excitative_frequ_sigma220q[1] = sigma_2mpq(2., 0., -1., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma220q[1], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma220q[1] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma220q[1] = kaula.0 ;
                // particle.tides.parameters.internal.im_love_number_sigma220_1 = kaula.1;
                // //println!("\nsigma201_1");
                particle.tides.parameters.internal.excitative_frequ_sigma201q[1] = sigma_2mpq(0., 1., -1., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma201q[1], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma201q[1] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma201q[1] = kaula.0 ;

                // println!("sigma2200");
                particle.tides.parameters.internal.excitative_frequ_sigma220q[2] = sigma_2mpq(2., 0., 0., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma220q[2], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma220q[2] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma220q[2] = kaula.0 ;

                particle.tides.parameters.internal.im_love_number_sigma2200 = kaula.1;
                particle.tides.parameters.internal.sigma2200_excitative_frequency = sigma_2mpq(2., 0., 0., spin, orbital_frequency);

                // println!("\nsigma2010");
                particle.tides.parameters.internal.excitative_frequ_sigma201q[2] = sigma_2mpq(0., 1., 0., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma201q[2], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma201q[2] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma201q[2] = kaula.0 ;

                // println!("sigma2201");
                particle.tides.parameters.internal.excitative_frequ_sigma220q[3] = sigma_2mpq(2., 0., 1., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma220q[3], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma220q[3] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma220q[3] = kaula.0 ;
                // particle.tides.parameters.internal.im_love_number_sigma2201 = kaula.1;
                // println!("\nsigma2011");
                particle.tides.parameters.internal.excitative_frequ_sigma201q[3] = sigma_2mpq(0., 1., 1., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma201q[3], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma201q[3] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma201q[3] = kaula.0 ;

                // println!("sigma2202");
                particle.tides.parameters.internal.excitative_frequ_sigma220q[4] = sigma_2mpq(2., 0., 2., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma220q[4], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma220q[4] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma220q[4] = kaula.0 ;
                // particle.tides.parameters.internal.im_love_number_sigma2202 = kaula.1;
                // println!("\nsigma2012");
                particle.tides.parameters.internal.excitative_frequ_sigma201q[4] = sigma_2mpq(0., 1., 2., spin, orbital_frequency);
                let kaula = kaula_number(particle.tides.parameters.internal.excitative_frequ_sigma201q[4], nm_data, re_k2, im_k2, w_lmpq );
                particle.tides.parameters.internal.im_love_number_sigma201q[4] = kaula.1 ;
                particle.tides.parameters.internal.re_love_number_sigma201q[4] = kaula.0 ;
                // for q in 0..5{
                    // println!(" sigma 220 {}     im {}       re {} ", q, particle.tides.parameters.internal.im_love_number_sigma220q[q], particle.tides.parameters.internal.re_love_number_sigma220q[q]);
                    // println!(" sigma 201 {}     im {}       re {} ", q, particle.tides.parameters.internal.im_love_number_sigma201q[q], particle.tides.parameters.internal.re_love_number_sigma201q[q]);
                // }
            }

            let mut sum_g_im_k2_over_q:f64=0.;
            for x in 0..5{
                sum_g_im_k2_over_q = sum_g_im_k2_over_q + eccentricity_function[x].powi(2) * particle.tides.parameters.internal.im_love_number_sigma220q[x] ;
            }

            // let sum_g_im_k2_over_q = eccentricity_function[0].powi(2) * particle.tides.parameters.internal.im_love_number_sigma220q[0] +
            //                         eccentricity_function[1].powi(2) * particle.tides.parameters.internal.im_love_number_sigma220q[1] +
            //                         eccentricity_function[2].powi(2) * particle.tides.parameters.internal.im_love_number_sigma220q[2] +
            //                         eccentricity_function[3].powi(2) * particle.tides.parameters.internal.im_love_number_sigma220q[3] +
            //                         eccentricity_function[4].powi(2) * particle.tides.parameters.internal.im_love_number_sigma220q[4];

            // for x in 0..5{
            //     excitative_frequency = sigma_2mpq(2., 0., q, spin, orbital_frequency);
            //     // excitative_frequency = 2.0* (orbital_frequency - spin) + q*orbital_frequency;

            //     let im_kaula_number = kaula_number(excitative_frequency, nm_data, re_k2, im_k2, w_lmpq);
            //     sum_g_im_k2_over_q = sum_g_im_k2_over_q + ( eccentricity_function[0][x].powi(2) * im_kaula_number.1 );

            //     // ////println!("The Spin {} and Orbital frequency {} ", spin, orbital_frequency);
            //     // ////println!("\tThe excitative frequ {} give the im_Kauma_num {}", excitative_frequency, im_kaula_number.1);
            //     // ////println!("\tThe sum_G*im_k2: {} \t the im_Kaula {}, \t the ecc_function_G {} \n", sum_g_im_k2_over_q, im_kaula_number.1, eccentricity_function[0][x]);
            //     q+=1.;
            // }

            //////println!("\t\nThe Total  sum_over_q into TidalTorque {} \n", sum_g_im_k2_over_q );

            tidal_torque_kaula_coplanar[1] = - (3./2.)*( ( G* tidal_host_particle.mass.powi(2)* particle.radius.powi(5) ) / semi_major_axis.powi(6) ) *sum_g_im_k2_over_q;

            // println!("\n\tInto tidal Torque:");
            // println!("Ecc {}, semi_maj_axis {} AU, orbital period {} Days", eccentricity, semi_major_axis, orbital_period);
            // println!("The Star Mass {} Msol, The planet radius {} Rearth, ", tidal_host_particle.mass, particle.radius*AU/6.3781e6 );
            // panic!("...");

            //////println!("Star mass {}, particle radius {} semi_maj_axis {}",tidal_host_particle.mass, particle.radius, semi_major_axis );
            //////println!("Tidal torque [radial, ortho, azimuth] {:?}",tidal_torque_kaula_coplanar );


            //// Torque calculation (star)
            let torque_due_to_tides_x: f64 = 0.;

            let torque_due_to_tides_y: f64 = 0.;

            let torque_due_to_tidez_z: f64 = tidal_torque_kaula_coplanar[1];

            let factor = -1.0;

            if central_body {
                // Integration of the spin (total torque tides):
                dangular_momentum_dt.x += factor * torque_due_to_tides_x;
                dangular_momentum_dt.y += factor * torque_due_to_tides_y;
                dangular_momentum_dt.z += factor * torque_due_to_tidez_z;
            } else {
                particle.tides.parameters.output.dangular_momentum_dt.x =  factor * torque_due_to_tides_x;
                particle.tides.parameters.output.dangular_momentum_dt.y =  factor * torque_due_to_tides_y;
                particle.tides.parameters.output.dangular_momentum_dt.z =  factor * torque_due_to_tidez_z; 
                //println!("dRotationnal_Angular_momentum dt in particle {:?}",particle.tides.parameters.output.dangular_momentum_dt );
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
    //////println!("dRotationnal_Angular_momentum dt in central body {:?}",tidal_host_particle.tides.parameters.output.dangular_momentum_dt );
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
           //////println!("calculate_orthogonal_component_of_the_tidal_force");
           if !central_body{
                let orbital_elements = tools::calculate_keplerian_orbital_elements(G*(tidal_host_particle.mass+particle.mass), particle.heliocentric_position, particle.heliocentric_velocity);
                //keplerian element return (a, q, eccentricity, i, p, n, l, orbital_period)
                let semi_major_axis = orbital_elements.0;
                let eccentricity = orbital_elements.2;
                let orbital_period = orbital_elements.7;
                let orbital_frequency = TWO_PI / (orbital_period *DAY) ;
            
                let _spin = particle.norm_spin_vector_2.sqrt()/DAY;
                let _imaginary_k2 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.imaginary_part_love_number;
                let _real_k2 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.real_part_love_number;
                let _w_lmpq = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.love_number_excitation_frequency;
                let _nm_data = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.num_datapoints;
                // num_datapoint is f64
                
                let eccentricity_function_g_20q = eccentricty_function_g_20q(eccentricity);

                let mut sum_over:f64 = 0.;
                let n = orbital_frequency;
                let t = current_time;
                let mut integer_q = -2.;
                let mut integer_j = -2.;
                // println!("\n\tInto Radial tidal force:");
                // println!("Ecc {}, semi_maj_axis {} AU, orbital period {} Days", eccentricity, semi_major_axis, orbital_period);
                // println!("The Star Mass {} Msol, The planet radius {} Rearth, ", tidal_host_particle.mass, particle.radius*AU/6.3781e6 );
                // panic!("...");
    
                for q in 0..5 {
                    for j in 0..5{
                        let phase_term = (integer_q - integer_j)*n*t;
                        sum_over = sum_over + eccentricity_function_g_20q[q] * eccentricity_function_g_20q[j] 
                                                * ( phase_term.sin() * particle.tides.parameters.internal.re_love_number_sigma220q[q] 
                                                - phase_term.cos() * particle.tides.parameters.internal.im_love_number_sigma220q[q] );

                        integer_j = integer_j +1.;
                        // println!("n° {}     real {}   imaginary {}", q, particle.tides.parameters.internal.re_love_number_sigma220q[q], particle.tides.parameters.internal.im_love_number_sigma220q[q] );
                    }
                    integer_q = integer_q +1.;
                }
                particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = - (3./2.) *( (G* tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(6) * particle.heliocentric_distance) ) * sum_over;
                // println!("Ortho force {}", particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide);
                // for q in 0..5 {
                //     for j in 0..5{
                //         let phase_term = (integer_q - integer_j)*n*t;
                //         sum_over = eccentricity_function_g_20q[q] * eccentricity_function_g_20q[j] * 
                //                                 ( -phase_term.sin() * particle.tides.parameters.internal.re_love_number_sigma220q[q] + phase_term.cos() * particle.tides.parameters.internal.im_love_number_sigma220q[q] );

                //         integer_j = integer_j +1.;
                //     }
                //     integer_q = integer_q +1.;
                // }
                // particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = - (3./2.) *( (G* tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(6) * particle.heliocentric_distance) ) * sum_over

                // let mut excitative_frequency: f64;
                // let mut sum_over_q = 0.;
                // let mut q = -2.;

                //////println!("\nInto Ortho tidal force:");

                // let sum_over_q = eccentricity_function_g_20q[0].powi(2) * particle.tides.parameters.internal.im_love_number_sigma220q[0] +
                //                     eccentricity_function_g_20q[1].powi(2) * particle.tides.parameters.internal.im_love_number_sigma220q[1] +
                //                     eccentricity_function_g_20q[2].powi(2) * particle.tides.parameters.internal.im_love_number_sigma220q[2]+
                //                     eccentricity_function_g_20q[3].powi(2) * particle.tides.parameters.internal.im_love_number_sigma220q[3]+
                //                     eccentricity_function_g_20q[4].powi(2) * particle.tides.parameters.internal.im_love_number_sigma220q[4];

                // for x in 0..5{
                //     excitative_frequency = sigma_2mpq(2., 0., q, spin, orbital_frequency);
                //     let im_kaula_number = kaula_number(excitative_frequency, nm_data, real_k2, imaginary_k2, w_lmpq);
                //     //////println!("\tThe excitative frequ {} give the im_Kauma_num {}", excitative_frequency, im_kaula_number.1);
                //     sum_over_q = sum_over_q +( eccentricity_function[0][x].powf(2.) * im_kaula_number.1 );
                //     //////println!("\tThe sum_over_q: {} \t the im_Kaula {}, \t the ecc_function_G {} \n", sum_over_q, im_kaula_number.1, eccentricity_function[0][x]);
                //     q+=1.;
                // }
                // //println!("\t\nThe Total  sum_over_q into Ortho tidal force {} \n",sum_over );

                // particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 0.0; //ORTHO FORCE NULL
            }
        }
    }
    //panic!("Does it reach this point?-calculate_orthogonal_component_of_the_tidal_force_for");

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
            //////println!("calculate_radial_component_of_the_tidal_force");

            let orbital_elements = tools::calculate_keplerian_orbital_elements(G*(tidal_host_particle.mass+particle.mass), particle.heliocentric_position, particle.heliocentric_velocity);
            //calculate keplerian el return (a, q, eccentricity, i, p, n, l, orbital_period)
            let semi_major_axis = orbital_elements.0;
            let eccentricity = orbital_elements.2;

            let orbital_period = orbital_elements.7;
            let orbital_frequency = TWO_PI / (orbital_period *DAY) ;

            // let _spin = particle.norm_spin_vector_2.sqrt()/DAY; // in rad/s
            // //////println!("The Planetary Spin in [s^-1]: {}\t Period in [s]: {}\t in [day]: {} Vector: {:?}", spin, 1./(spin), 1./(spin*DAY), particle.spin);
            // let _imaginary_k2 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.imaginary_part_love_number;
            // let _real_k2 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.real_part_love_number;
            // let _w_lmpq = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.love_number_excitation_frequency;
            // let _nm_data = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.num_datapoints;

            let eccentricity_function_g_20q = eccentricty_function_g_20q(eccentricity);
            let eccentricity_function_g_21q = eccentricty_function_g_21q(eccentricity);

            // let mut _excitative_frequency: f64;
            let mut sum_over:f64 = 0.;
            let n = orbital_frequency;
            let t = current_time;
            let mut integer_q:f64 = -2.;
            let mut integer_j:f64 = -2.;

            // let im_k2 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.imaginary_part_love_number;
            // let re_k2 = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.real_part_love_number;
            // let w_lmpq = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.love_number_excitation_frequency;
            // let nm_data = particle.tides.parameters.input.kaula_coplanar_tides_input_parameters.num_datapoints;

            // println!("\n\tInto Ortho tidal force:");
            // println!("Ecc {}, semi_maj_axis {} AU, orbital period {} Days", eccentricity, semi_major_axis, orbital_period);
            // println!("The Star Mass {} Msol, The planet radius {} Rearth, ", tidal_host_particle.mass, particle.radius*AU/6.3781e6 );
            // panic!("...");

            for q in 0..5{
                for j in 0..5{
                    let phase_term =  (integer_q - integer_j)*n*t;
                    sum_over = sum_over + (3./4.) * eccentricity_function_g_21q[q] * eccentricity_function_g_21q[j]
                                           * ( phase_term.cos() * particle.tides.parameters.internal.re_love_number_sigma201q[q] -  phase_term.sin() * particle.tides.parameters.internal.im_love_number_sigma201q[q] )
                                + (9./4.) * eccentricity_function_g_20q[q] * eccentricity_function_g_20q[j]
                                           * ( phase_term.cos()*particle.tides.parameters.internal.re_love_number_sigma220q[q] -  phase_term.sin()*particle.tides.parameters.internal.im_love_number_sigma220q[q] ) ;
                    integer_j = integer_j +1.;
                }
                integer_q = integer_q +1.;
            }
            particle.tides.parameters.internal.radial_component_of_the_tidal_force = - ( (G* tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(6) * particle.heliocentric_distance) ) * sum_over;
            // println!("Radial force {}", particle.tides.parameters.internal.radial_component_of_the_tidal_force);


            // for q in 0..5{
            //     for j in 0..5{
            //         let phase_term =  (integer_q - integer_j)*n*t;
            //         let kaula = kaula_number(integer_q*n, nm_data, re_k2, im_k2, w_lmpq );
            //         sum_over = -(3./4.) *(1./semi_major_axis.powi(3)) *eccentricity_function_g_21q[q]
            //                                * ( (integer_q*n*t).cos() * kaula.0 + (integer_q*n*t).sin() * kaula.1 )
            //                     - (9./4.) *(1./particle.heliocentric_distance.powi(3)) * eccentricity_function_g_20q[q] * eccentricity_function_g_20q[j]
            //                                * ( phase_term.cos()*particle.tides.parameters.internal.re_love_number_sigma220q[q] + phase_term.sin()*particle.tides.parameters.internal.im_love_number_sigma220q[q] ) ;
            //         integer_j = integer_j +1.;
            //     }
            //     integer_q = integer_q +1.;
            // }
            // particle.tides.parameters.internal.radial_component_of_the_tidal_force = - ( (G* tidal_host_particle.mass.powi(2)* particle.radius.powi(5)) / (semi_major_axis.powi(3) * particle.heliocentric_distance) ) * sum_over;

            // let sum_over_q = - (3./4.)*( eccentricity_function_g_21q[0].powi(2) * particle.tides.parameters.internal.re_love_number_sigma220_2 +
            //                                 eccentricity_function_g_21q[1].powi(2) * particle.tides.parameters.internal.re_love_number_sigma220_1 +
            //                                 eccentricity_function_g_21q[2].powi(2) * particle.tides.parameters.internal.re_love_number_sigma2200 +
            //                                 eccentricity_function_g_21q[3].powi(2) * particle.tides.parameters.internal.re_love_number_sigma2201 +
            //                                 eccentricity_function_g_21q[4].powi(2) * particle.tides.parameters.internal.re_love_number_sigma2202 ) 
            //                 + (9./4.)*( eccentricity_function_g_20q[0].powi(2) * particle.tides.parameters.internal.re_love_number_sigma220_2 +
            //                                 eccentricity_function_g_20q[1].powi(2) * particle.tides.parameters.internal.re_love_number_sigma220_1 +
            //                                 eccentricity_function_g_20q[2].powi(2) * particle.tides.parameters.internal.re_love_number_sigma2200 +
            //                                 eccentricity_function_g_20q[3].powi(2) * particle.tides.parameters.internal.re_love_number_sigma2201 +
            //                                 eccentricity_function_g_20q[4].powi(2) * particle.tides.parameters.internal.re_love_number_sigma2202 
            //                             );

            // for x in 0..5{ 
            //     excitative_frequency = sigma_2mpq(0., 1., q, spin, orbital_frequency);
            //     let mut real_kaula_number = kaula_number(excitative_frequency, nm_data, real_k2, imaginary_k2, w_lmpq);
            //     //////println!("\tThe sigma201q_excit_frequ {} give the real_Kauma_num {}", excitative_frequency, real_kaula_number.0);
            //     sum_over_q = sum_over_q -(3./4.)*( eccentricity_function[1][x].powf(2.) * real_kaula_number.0 );
            //     excitative_frequency = sigma_2mpq(2., 0., q, spin, orbital_frequency);
            //     real_kaula_number = kaula_number( excitative_frequency, nm_data, real_k2, imaginary_k2, w_lmpq );
            //     //////println!("\tThe sigma220q_excit_frequ {} give the real_Kaula_num {}", excitative_frequency, real_kaula_number.0);
            //     sum_over_q = sum_over_q +(9./4.)*( eccentricity_function[0][x].powf(2.) * real_kaula_number.0 );
            //     //////println!("\tThe sum_over_q {} \t the Re_Kaula {}, \t the ecc_function_g {} \n ", sum_over_q, real_kaula_number.0, eccentricity_function[0][x]);
            //     q+=1.;
            // }
            // //println!("\t\nThe Total sum_over_q into Radial tidal force {} \n",sum_over );

            
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


pub fn calculate_tidal_acceleration(tidal_host_particle: &mut Particle, particles: &mut [Particle], more_particles: &mut [Particle]) {
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

            let radial_component_of_the_tidal_force = particle.tides.parameters.internal.radial_component_of_the_tidal_force;
            let orthogonal_component_of_the_tidal_force = particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide;

            let total_tidal_force_x = radial_component_of_the_tidal_force *(particle.tides.coordinates.position.x / particle.tides.parameters.internal.distance) - orthogonal_component_of_the_tidal_force *(particle.tides.coordinates.position.y / particle.tides.parameters.internal.distance) ;
            let total_tidal_force_y = radial_component_of_the_tidal_force *(particle.tides.coordinates.position.y / particle.tides.parameters.internal.distance) + orthogonal_component_of_the_tidal_force *(particle.tides.coordinates.position.x / particle.tides.parameters.internal.distance) ;
            let total_tidal_force_z = 0.0; 

            // //println!("\nThe componant of force radial [{}] ortho [{}]\n", radial_component_of_the_tidal_force, orthogonal_component_of_the_tidal_force );
            ////////println!("factor3 {:e} {:e} {:e}", particle.tides.parameters.internal.radial_component_of_the_tidal_force, particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide, particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide);
            ////////println!("d {:e} vrad {:e}", particle.tides.parameters.internal., particle.tides.parameters.internal.radial_velocity);
            ////////println!("total {:e} {:e} {:e}", total_tidal_force_x, total_tidal_force_y, total_tidal_force_z);

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

// pub fn eccentricty_function_g(eccentricity: f64) -> [[f64; 5]; 3]{

//     let mut eccentricity_function = [[0.9;5];3];
//     // in this notation: 0 is q=-2, 1 is q=-1, 2 is q=0, 3 is q=1, 4 is q=2 ... cf Kaula 64 Table 3

//     eccentricity_function[0][0] = 0.;
//     eccentricity_function[0][1] = -(1./2.)*eccentricity +(1./16.)*eccentricity.powi(3);
//     eccentricity_function[0][2] = 1. -(5./2.)*eccentricity.powi(2) +(13./16.)*eccentricity.powi(4);
//     eccentricity_function[0][3] = (7./2.)*eccentricity -(123./16.)*eccentricity.powi(3);
//     eccentricity_function[0][4] = (17./2.)*eccentricity.powi(2) -(115./6.)*eccentricity.powi(4);

//     eccentricity_function[1][0] = (9./4.)*eccentricity.powi(2) +(7./4.)*eccentricity.powi(4);
//     eccentricity_function[1][1] = (3./2.)*eccentricity +(27./16.)*eccentricity.powi(3);
//     eccentricity_function[1][2] = (1. -eccentricity.powi(2)).powf(-(3./2.));
//     eccentricity_function[1][3] = (3./2.)*eccentricity +(27./16.)*eccentricity.powi(3);
//     eccentricity_function[1][4] = (9./4.)*eccentricity.powi(2) +(7./4.)*eccentricity.powi(4);

//     eccentricity_function[2][0] = (17./2.)*eccentricity.powi(2) -(115./6.)*eccentricity.powi(4);
//     eccentricity_function[2][1] = (7./2.)*eccentricity -(123./16.)*eccentricity.powi(3);
//     eccentricity_function[2][2] = 1. -(5./2.)*eccentricity.powi(2) +(13./16.)*eccentricity.powi(4);
//     eccentricity_function[2][3] = -(1./2.)*eccentricity + (1./16.)*eccentricity.powi(3);
//     eccentricity_function[2][4] = 0.;

//     return eccentricity_function;
// }
pub fn eccentricty_function_g_20q(eccentricity: f64) -> [f64; 5]{

    let mut eccentricity_function = [0.9;5];
    // let e_1 = eccentricity ;
    // let e_2 = e_1 ;
    // let e_3 = e_2 * e_1 ;
    // let e_4 = e_3 *e_1 ;

    eccentricity_function[0] = 0.;
    eccentricity_function[1] = -(1./2.)*eccentricity +(1./16.)*eccentricity.powi(3);
    eccentricity_function[2] = 1. -(5./2.)*eccentricity.powi(2) +(13./16.)*eccentricity.powi(4);
    eccentricity_function[3] = (7./2.)*eccentricity -(123./16.)*eccentricity.powi(3);
    eccentricity_function[4] = (17./2.)*eccentricity.powi(2) -(115./6.)*eccentricity.powi(4);

    return eccentricity_function;
}
pub fn eccentricty_function_g_21q(eccentricity: f64) -> [f64; 5]{

    let mut eccentricity_function = [0.9;5];
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


pub fn kaula_number(wk2:f64, _nm_data:f64, real_part_love_number: [[f64;32];32], imaginary_part_love_number: [[f64;32];32], love_number_excitation_frequency: [[f64;32];32] ) -> (f64,f64){
    let mut w_k2 = wk2;
    let mut re_k2 = 0.;
    let mut im_k2 = 0.;
    let mut _x = 0.;
    let mut _y = 0.;
    let mut ctrl = true; 
    let mut parity = false;

    // //println!("NUM {}", nm_data);
    // //println!("TEST{} {}", nm_data%32., nm_data - (nm_data %32.)*32.);
    // let A :i32 = (nm_data%32. -1.) as i32;
    // let B :i32 = (nm_data - (nm_data %32.)*32. -1.) as i32;
    // //println!("TEST {} {} GIVE {} , {}", A, B, love_number_excitation_frequency [14][13], love_number_excitation_frequency [14][12]);
    // ////println!("frequ one [0][1] {} frequ two [1][0]{} ", love_number_excitation_frequency[0][1], love_number_excitation_frequency[1][0]);
    // if  love_number_excitation_frequency[0][1] > love_number_excitation_frequency[1][0]{
    //     ////println!("frequ one [0][1] > [1][0]");
    // }
    // if  love_number_excitation_frequency[0][1] < love_number_excitation_frequency[1][0]{
    //     ////println!("frequ one [0][1] < [1][0]");
    // }

    // println!("Search K2 for w_k2 {}",w_k2);
    if w_k2 < 0.0 {
        w_k2 = w_k2.abs();
        parity = true;
        // //println!("\tThe new wk2 {}", w_k2);
    }

    // // ////println!("Ctrl de ses mort:\n IM  First {} End {} \n RE First {} End {} \n Frequ First {} End {}",imaginary_part_love_number[0][0], imaginary_part_love_number[14][11], real_part_love_number[0][0], real_part_love_number[14][11], love_number_excitation_frequency[0][0], love_number_excitation_frequency[14][11] );
    // if w_k2 < love_number_excitation_frequency[0][0] { //HAVE TO CHANGE THIS CONDITION
    //     re_k2 = real_part_love_number[0][0];
    //     im_k2 = imaginary_part_love_number[0][0];
    // // } else if w_k2 > love_number_excitation_frequency[31][31] {
    // //     re_k2 = real_part_love_number[31][31];
    // //     im_k2 = imaginary_part_love_number[31][31];
    // } else {


    if w_k2 > love_number_excitation_frequency[14][12] {
        re_k2 = real_part_love_number[14][12];
        im_k2 = imaginary_part_love_number[14][12];
    } else if w_k2 > love_number_excitation_frequency[13][31] {
        if ctrl && love_number_excitation_frequency[14][0] > w_k2 {
            re_k2 = real_part_love_number[14][0] + (real_part_love_number[14][0] - real_part_love_number[14-1][31] )/2.;
            im_k2 = imaginary_part_love_number[14][0] + (imaginary_part_love_number[14][0] - imaginary_part_love_number[14-1][31] )/2.;
        } else if ctrl {
            for frequency2 in 0..13{
                if ctrl && love_number_excitation_frequency[14][frequency2] >= w_k2 {
                    if  love_number_excitation_frequency[14][frequency2] == w_k2 {
                        re_k2 = real_part_love_number[14][frequency2];
                        im_k2 = imaginary_part_love_number[14][frequency2];
                        ctrl = false;
                        break;
                    } else if ctrl {
                        re_k2 = real_part_love_number[14][frequency2-1] + (real_part_love_number[14][frequency2] - real_part_love_number[14][frequency2-1])/2.;
                        im_k2 = imaginary_part_love_number[14][frequency2-1] + (imaginary_part_love_number[14][frequency2] - imaginary_part_love_number[14][frequency2-1])/2.;
                        // re_k2 = real_part_love_number[14][frequency2-1] + (real_part_love_number[14][2] - real_part_love_number[14][1])/2.;
                        // im_k2 = imaginary_part_love_number[14][1] + (imaginary_part_love_number[14][2] - imaginary_part_love_number[14][1])/2.;
                        // //println!("\n \n Find {}, with {} <-> {}  ", im_k2, imaginary_part_love_number[14][frequency2], imaginary_part_love_number[14][frequency2-1]);
                        // //println!("Stuff {}", w_k2);
                        ctrl = false;
                        break;
                    }
                    if !ctrl {break;}
                }
                if !ctrl {break;}
            }
        }
    } else {
        for frequency1 in 0..love_number_excitation_frequency.len()-1{
            // //println!("here?");
            if ctrl && love_number_excitation_frequency[frequency1][31] > w_k2 {
                // //println!("the firts {}", love_number_excitation_frequency[frequency1][31]);
                if ctrl && love_number_excitation_frequency[frequency1][0] > w_k2 {
                    re_k2 = real_part_love_number[frequency1][0] + (real_part_love_number[frequency1][0] - real_part_love_number[frequency1 -1][31] )/2.;
                    im_k2 = imaginary_part_love_number[frequency1][0] + (imaginary_part_love_number[frequency1][0] - imaginary_part_love_number[frequency1 -1][31] )/2.;
                } else if ctrl {
                    for frequency2 in 0..love_number_excitation_frequency.len(){
                        if ctrl && love_number_excitation_frequency[frequency1][frequency2] >= w_k2 { 
                            if  love_number_excitation_frequency[frequency1][frequency2] == w_k2 {
                                re_k2 = real_part_love_number[frequency1][frequency2];
                                im_k2 = imaginary_part_love_number[frequency1][frequency2];
                                ctrl = false;
                                break;
                            }
                            else if ctrl { 
                                re_k2 = real_part_love_number[frequency1][frequency2-1] + (real_part_love_number[frequency1][frequency2] - real_part_love_number[frequency1][frequency2-1])/2.;
                                im_k2 = imaginary_part_love_number[frequency1][frequency2-1] + (imaginary_part_love_number[frequency1][frequency2] - imaginary_part_love_number[frequency1][frequency2-1])/2.;
                                // //println!("\n \n Find {}, with {} <-> {}  ", im_k2, imaginary_part_love_number[frequency1][frequency2], imaginary_part_love_number[frequency1][frequency2-1]);
                                ctrl = false;
                                break;
                            }
                        }
                        if !ctrl {break;}
                    }
                }
            }
            // ////println!("\nEND THE BOUCLE Y -{}", y);
            // y = y +1.;
            if !ctrl {break;}
        }
    }
    if parity {
        im_k2 = -im_k2;
    }
    // println!("The two number RE:{}  Im:{}\n", re_k2, im_k2);
    return (re_k2*1., im_k2*1.)
}


pub fn sigma_2mpq(m:f64, p:f64, q:f64, spin:f64, orbital_frequency: f64) -> f64{
    // //////println!("\tsigma_calulation: the orbital_frequ {} and spin {}  give the sigma excit_frequ {}, the integer m {}, p {}, q {},  ", orbital_frequency, spin, (2. -2.*p + q)*orbital_frequency -m*spin, m, p, q );
    return (2.0 -2.0*p +q)*orbital_frequency -m*spin;
}

