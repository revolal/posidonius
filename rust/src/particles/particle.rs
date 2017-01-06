use super::super::constants::{K2, SUN_DYN_FREQ, TWO_PI, R_SUN, R_EARTH, M2EARTH};
use super::super::rand::random;
use super::{Evolver, EvolutionType, SolarEvolutionType};
use super::{Axes};
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Particle {
    pub id: String, // Unique internal identifier
    pub mass: f64,
    pub mass_g: f64,
    pub radius: f64,
    pub scaled_dissipation_factor: f64, // sigma
    pub dissipation_factor_scale: f64, // to scale the dissipation factor (multiply)
    pub radius_of_gyration_2: f64,  // radius of gyration square can be computed in terms of the mass moment of inertia, which 
                                    // depends on the shape of the body and determines the torque needed for a desired angular acceleration
    pub love_number: f64,   // k. Dimensionless parameters that measure the rigidity of a planetary body and the 
                            // susceptibility of its shape to change in response to a tidal potential.
    pub fluid_love_number: f64,   // love number for a completely fluid planet (used for rotational flattening effects)
    //
    pub position: Axes,
    pub velocity: Axes,
    pub acceleration: Axes,
    //
    pub radial_velocity: f64,
    pub norm_velocity_vector: f64,
    pub norm_velocity_vector_2: f64,
    pub distance: f64,
    // Tides
    pub orthogonal_component_of_the_tidal_force_due_to_stellar_tide: f64,
    pub orthogonal_component_of_the_tidal_force_due_to_planetary_tide: f64,
    pub radial_component_of_the_tidal_force: f64,
    pub denergy_dt: f64,
    pub torque: Axes, // Force
    pub spin: Axes,
    pub dspin_dt: Axes,
    pub tidal_acceleration: Axes,
    // Rotational flattening
    pub acceleration_induced_by_rotational_flattering: Axes,
    // General Relativity
    pub general_relativity_factor: f64,
    pub general_relativity_acceleration: Axes,
    // Evolution
    pub evolver: Evolver,
    pub lag_angle: f64, // MathisSolarLike
    pub planet_dependent_dissipation_factors : HashMap<String, f64>,
}

impl Particle {
    pub fn new(mass: f64, radius: f64, dissipation_factor: f64, dissipation_factor_scale: f64, radius_of_gyration_2: f64, love_number: f64, fluid_love_number: f64, position: Axes, velocity: Axes, acceleration: Axes, spin: Axes, evolution_type: EvolutionType) -> Particle {
        let torque = Axes{x: 0., y: 0., z: 0.};
        let dspin_dt = Axes{x: 0., y: 0., z: 0.};
        let tidal_acceleration = Axes{x: 0., y: 0., z: 0.};
        let acceleration_induced_by_rotational_flattering = Axes{x: 0., y: 0., z: 0.};
        let general_relativity_acceleration = Axes{x: 0., y: 0., z: 0.};
        let evolver = Evolver::new(evolution_type);
        let id = Particle::rand_string(16);
        match evolution_type {
            EvolutionType::BrownDwarf(_) => println!("WARNING: Bodies with BrownDwarf evolution will ignore initial radius and radius of gyration."),
            EvolutionType::MDwarf => println!("WARNING: Bodies with MDwarf evolution will ignore initial radius."),
            EvolutionType::SolarLike(model) => {
                match model {
                    SolarEvolutionType::ConstantDissipation => println!("WARNING: Bodies with SolarLike evolution will ignore initial radius."),
                    SolarEvolutionType::EvolvingDissipation(_) => println!("WARNING: Bodies with MathisSolarLike evolution will ignore initial radius and dissipation factor."),
                }
            },
            EvolutionType::Jupiter => println!("WARNING: Bodies with Jupiter evolution will ignore initial radius, radius of gyration and love number."),
            EvolutionType::NonEvolving => {},
        }
        Particle { id:id, mass:mass, mass_g: mass*K2, radius: radius, 
                    scaled_dissipation_factor:dissipation_factor_scale*dissipation_factor, 
                    dissipation_factor_scale:dissipation_factor_scale,
                    radius_of_gyration_2:radius_of_gyration_2, 
                    love_number:love_number, fluid_love_number:fluid_love_number,
                    position:position, velocity:velocity, acceleration:acceleration, spin:spin,
                    radial_velocity: 0., norm_velocity_vector:0., norm_velocity_vector_2:0., distance:0.,
                    orthogonal_component_of_the_tidal_force_due_to_stellar_tide:0., orthogonal_component_of_the_tidal_force_due_to_planetary_tide:0., radial_component_of_the_tidal_force:0., 
                    denergy_dt:0., torque:torque, dspin_dt:dspin_dt, 
                    tidal_acceleration:tidal_acceleration, 
                    acceleration_induced_by_rotational_flattering:acceleration_induced_by_rotational_flattering,
                    general_relativity_acceleration:general_relativity_acceleration,
                    general_relativity_factor: 0.,
                    evolver:evolver,
                    lag_angle:0., // It will be initialized the first time evolve is called
                    planet_dependent_dissipation_factors:HashMap::new(),
        }
    }

    pub fn new_brown_dwarf(mass: f64, dissipation_factor_scale: f64, position: Axes, velocity: Axes, acceleration: Axes, evolution_type: EvolutionType) -> Particle {
        let (rotation_period, love_number) = match evolution_type {
            EvolutionType::NonEvolving => { 
                let rotation_period: f64 = 70.0; // hours
                let love_number: f64 = 0.307; // BrownDwarf
                (rotation_period, love_number)
            },
            EvolutionType::BrownDwarf(mass) => { 
                let rotation_period; // hours
                let love_number;
                if mass <= 0.0101 && mass >= 0.0099 {
                    rotation_period = 8.0;
                    love_number = 0.3790;
                } else if mass <= 0.0121 && mass >= 0.0119 {
                    rotation_period = 13.0;
                    love_number = 0.3780;
                } else if mass <= 0.0151 && mass >= 0.0149 {
                    rotation_period = 19.0;
                    love_number = 0.3760;
                } else if mass <= 0.0201 && mass >= 0.0199 {
                    rotation_period = 24.0;
                    love_number = 0.3690;
                } else if mass <= 0.0301 && mass >= 0.0299 {
                    rotation_period = 30.0;
                    love_number = 0.3550;
                } else if mass <= 0.0401 && mass >= 0.0399 {
                    rotation_period = 36.0;
                    love_number = 0.3420;
                } else if mass <= 0.0501 && mass >= 0.0499 {
                    rotation_period = 41.0;
                    love_number = 0.3330;
                } else if mass <= 0.0601 && mass >= 0.0599 {
                    rotation_period = 47.0;
                    love_number = 0.3250;
                } else if mass <= 0.0701 && mass >= 0.0699 {
                    rotation_period = 53.0;
                    love_number = 0.3110;
                } else if mass <= 0.0721 && mass >= 0.0719 {
                    rotation_period = 58.0;
                    love_number = 0.3080;
                } else if mass <= 0.0751 && mass >= 0.0749 {
                    rotation_period = 64.0;
                    love_number = 0.3070;
                } else if mass <= 0.0801 && mass >= 0.0799 {
                    rotation_period = 70.0;
                    love_number = 0.3070;
                } else {
                    panic!("The evolution type BrownDwarf does not support a mass of {} Msun!", mass);
                }
                (rotation_period, love_number)
            },
            _ => { panic!("Evolution type should be solar like or non evolving to create a solar like body!"); }
        };

        let spin0 = TWO_PI/(rotation_period/24.); // days^-1
        let spin = Axes{x:0., y:0., z:spin0 };

        let fluid_love_number = love_number;
        // BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
        let dissipation_factor: f64 = 2.006*3.845764e4; // -60+64

        let radius_factor: f64 = 0.845649342247916;
        let radius: f64 = radius_factor * R_SUN;
        let radius_of_gyration_2: f64 = 1.94e-1; // Brown dwarf
        let mut brown_dwarf = Particle::new(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number,
                                                position, velocity, acceleration, spin, evolution_type);
        let current_time = 0.;
        brown_dwarf.evolve(current_time); // If it is NonEvolving type, this method will do nothing
        brown_dwarf
    }

    pub fn new_solar_like(mass: f64, dissipation_factor_scale: f64, position: Axes, velocity: Axes, acceleration: Axes, evolution_type: EvolutionType) -> Particle {
        match evolution_type {
            EvolutionType::SolarLike(_) => { },
            EvolutionType::NonEvolving => { },
            _ => { panic!("Evolution type should be solar like or non evolving to create a solar like body!"); }
        }

        let rotation_period = 8.; // hours
        let love_number = 0.03; // SolarLike
        let spin0 = TWO_PI/(rotation_period/24.); // days^-1
        let spin = Axes{x:0., y:0., z:spin0 };

        let fluid_love_number = love_number;
        // Sun-like-star: sigmast = 4.992e-66 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
        let dissipation_factor: f64 = 4.992*3.845764e-2; // -66+64

        let radius_factor: f64 = 1.;
        let radius: f64 = radius_factor * R_SUN;
        let radius_of_gyration_2: f64 = 5.9e-2; // Sun
        let mut solarlike_star = Particle::new(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number,
                                                position, velocity, acceleration, spin, evolution_type);
        let current_time = 0.;
        solarlike_star.evolve(current_time); // If it is NonEvolving type, this method will do nothing
        solarlike_star
    }

    pub fn new_m_dwarf(mass: f64, dissipation_factor_scale: f64, position: Axes, velocity: Axes, acceleration: Axes, evolution_type: EvolutionType) -> Particle {
        match evolution_type {
            EvolutionType::MDwarf => { },
            EvolutionType::NonEvolving => { },
            _ => { panic!("Evolution type should be M dwarf or non evolving to create a M dwarf body!"); }
        }

        let rotation_period = 70.; // hours
        let love_number: f64 = 0.307; // M Dwarf
        let fluid_love_number = love_number;

        let spin0 = TWO_PI/(rotation_period/24.); // days^-1
        let spin = Axes{x:0., y:0., z:spin0 };

        // BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
        let dissipation_factor: f64 = 2.006*3.845764e4; // -60+64

        let radius_factor: f64 = 0.845649342247916;
        let radius: f64 = radius_factor * R_SUN;
        let radius_of_gyration_2: f64 = 2.0e-1; // M-dwarf
        let mut m_dwarf = Particle::new(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number,
                                                position, velocity, acceleration, spin, evolution_type);
        let current_time = 0.;
        m_dwarf.evolve(current_time); // If it is NonEvolving type, this method will do nothing
        m_dwarf
    }

    pub fn new_jupiter_like(mass: f64, dissipation_factor_scale: f64, position: Axes, velocity: Axes, acceleration: Axes, evolution_type: EvolutionType) -> Particle {
        match evolution_type {
            EvolutionType::Jupiter => { },
            EvolutionType::NonEvolving => { },
            _ => { panic!("Evolution type should be jupyter or non evolving to create a jupiter body!"); }
        }

        let rotation_period = 9.8; // hours
        let love_number: f64 = 0.299; // Gas giant
        let fluid_love_number = love_number;

        let spin0 = TWO_PI/(rotation_period/24.); // days^-1
        let spin = Axes{x:0., y:0., z:spin0 };

        let radius_factor: f64 = 10.9; // Jupiter in R_EARTH
        let radius: f64 = radius_factor * R_EARTH;

        // k2delta_t for Jupiter: 2-3d-2 s, here in day (Leconte)
        let k2pdelta : f64 = 2.893519e-7;
        let dissipation_factor: f64 = 2.0 * K2 * k2pdelta / (3.0 * radius.powi(5));

        let radius_of_gyration_2: f64 = 3.308e-1; // Gas giant

        let mut jupiter = Particle::new(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number,
                                                position, velocity, acceleration, spin, evolution_type);
        let current_time = 0.;
        jupiter.evolve(current_time); // If it is NonEvolving type, this method will do nothing
        jupiter
    }

    pub fn new_earth_like(mass: f64, dissipation_factor_scale: f64, position: Axes, velocity: Axes, acceleration: Axes) -> Particle {
        let evolution_type = EvolutionType::NonEvolving;

        let rotation_period = 24.; // hours
        let love_number: f64 = 0.305; // Earth
        let fluid_love_number = 0.9532; // Earth

        let spin0 = TWO_PI/(rotation_period/24.); // days^-1
        let spin = Axes{x:0., y:0., z:spin0 };

        // Earth-like => mass-radius relationship from Fortney 2007
        let radius_factor : f64 = (0.0592*0.7+0.0975) * (mass.log10() + M2EARTH.log10() - K2.log10()).powi(2) 
                                 + (0.2337*0.7+0.4938) * (mass.log10() + M2EARTH.log10() - K2.log10()) 
                                 + 0.3102*0.7+0.7932;
        let radius: f64 = radius_factor * R_EARTH;
        let radius_of_gyration_2: f64 = 3.308e-1; // Earth type planet
        let k2pdelta: f64 = 2.465278e-3; // Terrestrial planets
        let dissipation_factor: f64 = 2. * K2 * k2pdelta/(3. * radius.powi(5));

        let mut earth_like = Particle::new(mass, radius, dissipation_factor, dissipation_factor_scale, radius_of_gyration_2, love_number, fluid_love_number,
                                                position, velocity, acceleration, spin, evolution_type);
        let current_time = 0.;
        earth_like.evolve(current_time); // If it is NonEvolving type, this method will do nothing
        earth_like
    }

    pub fn new_terrestrial(mass: f64, dissipation_factor_scale: f64, position: Axes, velocity: Axes, acceleration: Axes) -> Particle {
        let mut terrestrial = Particle::new_earth_like(mass, dissipation_factor_scale, position, velocity, acceleration);
        let radius_factor = 1.;
        terrestrial.radius = radius_factor * R_EARTH;
        terrestrial
    }
    
    pub fn new_gas_giant(mass: f64, dissipation_factor_scale: f64, position: Axes, velocity: Axes, acceleration: Axes) -> Particle {
        let evolution_type = EvolutionType::NonEvolving;
        let mut gas_giant = Particle::new_jupiter_like(mass, dissipation_factor_scale, position, velocity, acceleration, evolution_type);
        let dissipation_factor: f64 = 2.006*3.845764e4; // Gas giant
        gas_giant.scaled_dissipation_factor = gas_giant.dissipation_factor_scale * dissipation_factor;
        gas_giant
    }

    fn rand_string(n_characters: usize) -> String {
        (0..n_characters).map(|_| (0x20u8 + (random::<f32>() * 96.0) as u8) as char).collect()
    }

    pub fn evolve(&mut self, current_time: f64) {
        self.radius = self.evolver.radius(current_time, self.radius);
        self.radius_of_gyration_2 = self.evolver.radius_of_gyration_2(current_time, self.radius_of_gyration_2);
        self.love_number = self.evolver.love_number(current_time, self.love_number);
        self.lag_angle = match self.evolver.evolution_type {
                EvolutionType::SolarLike(model) => {
                    match model {
                        SolarEvolutionType::EvolvingDissipation(_) => {
                            let inverse_tidal_q_factor = self.evolver.inverse_tidal_q_factor(current_time, 0.);
                            //
                            // Calculation of the norm square of the spin for the star
                            let normspin_2 = self.spin.x.powi(2) + self.spin.y.powi(2) + self.spin.z.powi(2);
                            let epsilon_squared = normspin_2/SUN_DYN_FREQ;
                            // Normal formula = 3.d0*epsilon_squared*Q_str_inv/(4.d0*k2s)
                            // but as for sigma it is necessary to divide by k2s, we do not divide here
                            let lag_angle = 3.0*epsilon_squared*inverse_tidal_q_factor/4.0;
                            lag_angle
                        },
                        _ => 0.,
                    }
                },
                _ => 0.,
        };
        //println!("[{}] Evolve Radius {:e} Gyration {:e} Love {:e} Lag {:e}", current_time, self.radius, self.radius_of_gyration_2, self.love_number, self.lag_angle);
    }

    pub fn planet_dependent_dissipation_factor(&self, id: &String) -> f64 {
        match self.evolver.evolution_type {
            EvolutionType::SolarLike(model) => {
                match model {
                    SolarEvolutionType::EvolvingDissipation(_) => {
                        match self.planet_dependent_dissipation_factors.get(id) {
                            Some(&value) => value,
                            _ => self.scaled_dissipation_factor // This should not happen
                        }
                    },
                    _ => self.scaled_dissipation_factor,
                }
            },
            _ => self.scaled_dissipation_factor,
        }
    }


}


