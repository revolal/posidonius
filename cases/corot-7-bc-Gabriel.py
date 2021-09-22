import posidonius
import numpy as np
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('output_filename', action='store', help='Filename where the initial snapshot will be stored (e.g., universe_integrator.json)')

    args = parser.parse_args()
    filename = args.output_filename
    #filename = posidonius.constants.BASE_DIR+"target/case7.json"

    #initial_time = 4.5e6*365.25 # time [days] where simulation starts
    initial_time = 0.0 # time [days] where simulation starts
    time_step = 0.03 # days
    #time_step = 0.05 # days
    #time_limit   = 4*time_step # days
    #time_limit   = 365.25 * 1.0e8 # days
    time_limit   = 365.25*1.0e8 # days
    historic_snapshot_period = 100.0*1.0e4*time_step # days
    recovery_snapshot_period = 1000.*historic_snapshot_period # days
    consider_effects = posidonius.ConsiderEffects({
        "tides": True,
        "rotational_flattening": False,
        "general_relativity": False,
        "disk": False,
        "wind": False,
        "evolution": False,
    })
    universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects)

    star_mass = 0.915 # Solar masses
    star_radius_factor = 0.2134
    star_radius = star_radius_factor * posidonius.constants.R_SUN
    star_radius_of_gyration = 4.47e-01 # Brown dwarf
    star_position = posidonius.Axes(0., 0., 0.)
    star_velocity = posidonius.Axes(0., 0., 0.)

    # Initialization of stellar spin
    star_rotation_period = 3.3*24 # hours
    #star_rotation_period = 19.0*24 # hours
    star_angular_frequency = posidonius.constants.TWO_PI/(star_rotation_period/24.) # days^-1
    star_spin = posidonius.Axes(0., 0., star_angular_frequency)

    star_tides_parameters = {
        "dissipation_factor_scale": 0.0,
        "dissipation_factor": 0.0,
        "love_number": 0.0,
    }
    star_tides = posidonius.effects.tides.CentralBody(star_tides_parameters)
    #star_tides = posidonius.effects.tides.OrbitingBody(star_tides_parameters)
    #star_tides = posidonius.effects.tides.Disabled()
    #
    star_rotational_flattening = posidonius.effects.rotational_flattening.Disabled()
    #
    star_general_relativity = posidonius.effects.general_relativity.Disabled()
    #
    star_wind = posidonius.effects.wind.Disabled()
    #
    star_disk = posidonius.effects.disk.Disabled()
    #
    star_evolution = posidonius.NonEvolving()
    #
    star = posidonius.Particle(star_mass, star_radius, star_radius_of_gyration, star_position, star_velocity, star_spin)
    star.set_tides(star_tides)
    star.set_rotational_flattening(star_rotational_flattening)
    star.set_general_relativity(star_general_relativity)
    star.set_wind(star_wind)
    star.set_disk(star_disk)
    star.set_evolution(star_evolution)
    universe.add_particle(star)


    ############################################################################
    # Radiuses in R_EARTH
    planet_radiuses = (1.53, 2.17)
    # Masses in M_SUN
    planet_masses = (3.0e-6*4.73, 3.0e-6*13.56)
    # Semi-major axis in AU
    planet_a = (0.0188, 0.0455)                         
    # Inclination in degrees
    planet_i = (0.0, 0.0)
    # Mean anomaly in degrees
    planet_l = (0.0, 0.0)
    # Planet eccentricity
    planet_e = (0.25, 0.1546)
    # Planet argument of pericentre
    planet_w = (0.0, 0.0)
    # Uniform viscosity coefficients
    planet_uniform_viscosity_coefficients = (1.0e12, 1.0e12)        

    for r, m, a, i, l, e, p, eta in zip(planet_radiuses, planet_masses, planet_a, planet_i, planet_l, planet_e, planet_w, planet_uniform_viscosity_coefficients):
        planet_mass = m # Solar masses (3.0e-6 solar masses = 1 earth mass)
        planet_radius_factor = r # R Earth
        planet_radius = planet_radius_factor * posidonius.constants.R_EARTH
        planet_radius_of_gyration = 5.75e-01
        planet_uniform_viscosity_coefficient = eta

        #////////// Specify initial position and velocity for a stable orbit
        #////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
        i = i * posidonius.constants.DEG2RAD;                      # inclination (degrees)
        n = 0. * posidonius.constants.DEG2RAD;                      # longitude of the ascending node (degrees)
        l = l * posidonius.constants.DEG2RAD;           # mean anomaly (degrees)
        p = (p + n) * posidonius.constants.DEG2RAD;                 # Convert to longitude of perihelion !! (CHECK THIS for I not equal to 0)
        q = a * (1.0 - e);                     # perihelion distance
        planet_position, planet_velocity = posidonius.calculate_cartesian_coordinates(planet_mass, q, e, i, p, n, l, masses=[star_mass], positions=[star_position], velocities=[star_velocity])

        #////// Initialization of planetary spin
        planet_obliquity = 0.0 # rad
        # Pseudo-synchronization period
        planet_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(a, e, star_mass, planet_mass) # days
        planet_angular_frequency = 2.0*posidonius.constants.TWO_PI/(planet_pseudo_synchronization_period) # days^-1
        planet_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
        planet_inclination = planet_keplerian_orbital_elements[3]
        planet_spin = posidonius.calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity)

       # k2pdelta = 2.465278e-3 # Terrestrial planets (no gas)
       # planet_tides_parameters = {
       #     "dissipation_factor_scale": 1.0,
       #     "dissipation_factor": 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet_radius, 5)),
       #     "love_number": 0.299,
       # }

        planet_tides_parameters_creep = {
        "dissipation_factor_scale": 0.0,
        "dissipation_factor": 0.0,
        "love_number": 0.0,
        "uniform_viscosity_coefficient": (planet_uniform_viscosity_coefficient * posidonius.constants.AU * posidonius.constants.DAY / posidonius.constants.M_SUN),
        # value in kg s^-1 m^-1 => M_SUN day^-1 AU^-1 for the viscosity
        }
        #planet_tides = posidonius.effects.tides.CentralBody(planet_tides_parameters)
        #planet_tides = posidonius.effects.tides.OrbitingBody(planet_tides_parameters)
        planet_tides = posidonius.effects.tides.CreepCoplanarOrbitingBody(planet_tides_parameters_creep)
        #planet_tides = posidonius.effects.tides.Disabled()
        #
        planet_rotational_flattening = posidonius.effects.rotational_flattening.Disabled()
        #
        planet_general_relativity = posidonius.effects.general_relativity.Disabled()
        #
        planet_wind = posidonius.effects.wind.Disabled()
        #
        planet_disk = posidonius.effects.disk.Disabled()
        #
        planet_evolution = posidonius.NonEvolving()
        #
        planet = posidonius.Particle(planet_mass, planet_radius, planet_radius_of_gyration, planet_position, planet_velocity, planet_spin)
        planet.set_tides(planet_tides)
        planet.set_rotational_flattening(planet_rotational_flattening)
        planet.set_general_relativity(planet_general_relativity)
        planet.set_wind(planet_wind)
        planet.set_disk(planet_disk)
        planet.set_evolution(planet_evolution)
        universe.add_particle(planet)


    ############################################################################

    whfast_alternative_coordinates="DemocraticHeliocentric"
    #whfast_alternative_coordinates="WHDS"
    #whfast_alternative_coordinates="Jacobi"
    universe.write(filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)
    #universe.write(filename, integrator="IAS15")
    #universe.write(filename, integrator="LeapFrog")


