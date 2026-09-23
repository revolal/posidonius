import posidonius
import numpy as np
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('output_filename', action='store', help='Filename where the initial snapshot will be stored (e.g., universe_integrator.json)')

    args = parser.parse_args()
    filename = args.output_filename
    #filename = posidonius.constants.BASE_DIR+"target/example.json"

    initial_time = 4.5e6*365.25 # time [days] where simulation starts
    time_step = 0.01 # days
    #time_limit   = 4*time_step # days
    time_limit   =  365.25 * 1.0e10 # days
    historic_snapshot_period = 365.25 #* 10. #*time_step # days
    recovery_snapshot_period = 10.*historic_snapshot_period # days
    consider_effects = posidonius.ConsiderEffects({
        "tides": True,
        "rotational_flattening": False,
        "general_relativity": False,
        "disk": False,
        "wind": False,
        "evolution": False,
    })
    universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects)

    # star_mass = 1.0 #0.08 # Solar masses
    # star_radius_factor = 0.845649342247916
    # star_radius = star_radius_factor * posidonius.constants.R_SUN
    # star_radius_of_gyration = 4.47e-01 # Brown dwarf
    star_mass = 1.#0.089 # Solar masses
    star_radius_factor = 1.#0.117
    star_radius = star_radius_factor * posidonius.constants.R_SUN
    star_radius_of_gyration = 4.47e-01 # Brown dwarf
    star_position = posidonius.Axes(0., 0., 0.)
    star_velocity = posidonius.Axes(0., 0., 0.)

    star_position = posidonius.Axes(0., 0., 0.)
    star_velocity = posidonius.Axes(0., 0., 0.)

    # # Initialization of stellar spin
    # star_rotation_period = 70.0 # hours
    # star_angular_frequency = posidonius.constants.TWO_PI/(star_rotation_period/24.) # days^-1
    # star_spin = posidonius.Axes(0., 0., star_angular_frequency)
    # Initialization of stellar spin
    star_rotation_period = 3.3*24 # hours
    #star_rotation_period = 19.0*24 # hours
    star_angular_frequency = posidonius.constants.TWO_PI/(star_rotation_period/24.) # days^-1
    star_spin = posidonius.Axes(0., 0., star_angular_frequency)

    # star_tides_parameters = {
    #     "dissipation_factor_scale": 1.0,
    #     "dissipation_factor": 2.006*3.845764e4,
    #     "love_number": 0.307,
    # }
    star_tides_parameters = {
        "dissipation_factor_scale": 0.0,
        "dissipation_factor": 0.0,
        "love_number": 0.0,
    }
    star_tides_model = posidonius.effects.tides.ConstantTimeLag(star_tides_parameters)
    star_tides = posidonius.effects.tides.CentralBody(star_tides_model)
    #star_tides = posidonius.effects.tides.OrbitingBody(star_tides_model)
    # star_tides = posidonius.effects.tides.Disabled()
    #
    # star_rotational_flattening_parameters = {"love_number": star_tides_parameters["love_number"]}
    # star_rotational_flattening_parameters =  {"uniform_viscosity_coefficient": 1e14}
    # star_rotational_flattening_model = posidonius.effects.rotational_flattening.OblateSpheroid(star_rotational_flattening_parameters)
    # star_rotational_flattening_model = posidonius.effects.rotational_flattening.CreepCoplanar(star_rotational_flattening_parameters)
    # star_rotational_flattening = posidonius.effects.rotational_flattening.CentralBody(star_rotational_flattening_model)
    #star_rotational_flattening = posidonius.effects.rotational_flattening.OrbitingBody(star_rotational_flattening_model)
    star_rotational_flattening = posidonius.effects.rotational_flattening.Disabled()
    #
    # star_general_relativity = posidonius.effects.general_relativity.CentralBody("Kidder1995")
    #star_general_relativity = posidonius.effects.general_relativity.CentralBody("Anderson1975")
    #star_general_relativity = posidonius.effects.general_relativity.CentralBody("Newhall1983")
    #star_general_relativity = posidonius.effects.general_relativity.OrbitingBody()
    star_general_relativity = posidonius.effects.general_relativity.Disabled()
    #
    #star_wind = posidonius.effects.wind.Interaction({
        ## Solar wind parametrisation (Bouvier 1997)
        #"k_factor": 4.0e-18, # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
        #"rotation_saturation": 1.7592918860102842, # 14. * TWO_PI/25.0, in units of the spin of the Sun today
    #})
    star_wind = posidonius.effects.wind.Disabled()
    #
    #disk_surface_density_normalization_gcm = 1000. # g.cm^-2
    #disk_surface_density_normalization_SI = disk_surface_density_normalization_gcm * 1.0e-3 * 1.0e4 # kg.m^-2
    #disk_properties = {
        #'inner_edge_distance': 0.01,  # AU
        #'outer_edge_distance': 100.0, # AU
        #'lifetime': 1.0e5 * 365.25e0, # days
        #'alpha': 1.0e-2,
        #'surface_density_normalization': disk_surface_density_normalization_SI * (1.0/posidonius.constants.M_SUN) * posidonius.constants.AU**2, # Msun.AU^-2
        #'mean_molecular_weight': 2.4,
    #}
    #star_disk = posidonius.effects.disk.CentralBody(disk_properties)
    #star_disk = posidonius.effects.disk.OrbitingBody()
    star_disk = posidonius.effects.disk.Disabled()
    #
    #star_evolution = posidonius.GalletBolmont2017(star_mass) # mass = 0.30 .. 1.40
    #star_evolution = posidonius.BolmontMathis2016(star_mass) # mass = 0.40 .. 1.40
    #star_evolution = posidonius.Baraffe2015(star_mass) # mass = 0.01 .. 1.40
    #star_evolution = posidonius.Leconte2011(star_mass) # mass = 0.01 .. 0.08
    #star_evolution = posidonius.Baraffe1998(star_mass) # Sun (mass = 1.0) or M-Dwarf (mass = 0.1)
    #star_evolution = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
    #star_evolution = posidonius.LeconteChabrier2013(True) # Jupiter with dissipation of dynamical tides
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
    Mearth = 5.9736E+24 #Earth Mass (kg)
    Rearth = 6378.1E+3    #Earth radius (m)
    Mvenus = 4.8685E+24
    Rvenus = 6052E+3
    Mjup = 1.898E+27 # Jupiter mass (kg)
    Rjup = 6.991E+07 # Jupiter radius (m)

    # planet_mass_factor = 1.0
    # planet_mass = planet_mass_factor * posidonius.constants.M_EARTH # Solar masses (3.0e-6 solar masses = 1 earth mass)

    # # Earth-like => mass-radius relationship from Fortney 2007
    # planet_radius_factor = posidonius.tools.mass_radius_relation(planet_mass_factor, planet_mass_type='factor', planet_percent_rock=0.70)
    # planet_radius = planet_radius_factor * posidonius.constants.R_EARTH
    # planet_radius_of_gyration = 5.75e-01 # Earth type planet

    planet_mass = Mvenus / posidonius.constants.M_SUN # Solar masses (3.0e-6 solar masses = 1 earth mass)
    planet_radius = Rvenus / posidonius.constants.AU

    #////////// Specify initial position and velocity for a stable orbit
    #////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    a = 0.072300062228760203                # semi-major axis (in AU)
    e = 0.15;                               # eccentricity
    i = 0. * posidonius.constants.DEG2RAD;  # inclination (degrees)
    p = 0. * posidonius.constants.DEG2RAD;  # argument of pericentre (degrees)
    n = 0. * posidonius.constants.DEG2RAD;  # longitude of the ascending node (degrees)
    l = 0. * posidonius.constants.DEG2RAD;  # mean anomaly (degrees)
    p = (p + n);                            # Convert to longitude of perihelion !!
    q = a * (1.0 - e);                      # perihelion distance
    planet_position, planet_velocity = posidonius.calculate_cartesian_coordinates(planet_mass, q, e, i, p, n, l, masses=[star_mass], positions=[star_position], velocities=[star_velocity])

    #////// Initialization of planetary spin
    planet_obliquity = 0. * posidonius.constants.DEG2RAD # 0.2 rad
    orbital_period = 2.*np.pi*(np.sqrt( (a*posidonius.constants.AU)**3. / ( posidonius.constants.G_SI*( star_mass*posidonius.constants.M_SUN + planet_mass*posidonius.constants.M_SUN ) ) ) )/86400.
    planet_rotation_period = orbital_period/(1.65) #24. # hours orbital_period/(1.0 + 1e-5) #24. # hours
    planet_angular_frequency = posidonius.constants.TWO_PI/(planet_rotation_period) # days^-1
    # Pseudo-synchronization period
    #planet_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
    #planet_semi_major_axis = planet_keplerian_orbital_elements[0]
    #planet_eccentricity = planet_keplerian_orbital_elements[2]
    #planet_semi_major_axis = a
    #planet_eccentricity = e
    #planet_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(planet_semi_major_axis, planet_eccentricity, star_mass, planet_mass)
    #planet_angular_frequency = posidonius.constants.TWO_PI/(planet_pseudo_synchronization_period/24.) # days^-1
    planet_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
    planet_inclination = planet_keplerian_orbital_elements[3]
    planet_spin = posidonius.calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity)

    # --- CTL's parameters --- #
    k2pdelta = 2.465278e-3 # Terrestrial planets (no gas)
    planet_tides_parameters = {
        "dissipation_factor_scale": 1.0,
        "dissipation_factor": 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet_radius, 5)),
        "love_number": 0.305,
    }

    # --- Kaula's developpement parameters --- #
    planet_data_venus = np.loadtxt('./cases/love_number_data/Results_V1_Armann2_atm_freq_Imk2_posidonius.txt',comments='#')
    # planet_data_h = np.loadtxt('Results_Trappist_1h_13.5_150_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_venus, planet_radius_venus, planet_gyration_radius_squared_venus = planet_data_venus[0,:]

    w_lm_venus   = planet_data_venus[1:,0]
    ImK2_venus   = planet_data_venus[1:,1]
    ReK2_venus   = planet_data_venus[1:,2]
    size_venus   = np.size(w_lm_venus)

    w_lmpq = w_lm_venus
    ImK2 = ImK2_venus
    ReK2 = ReK2_venus
    size = size_venus
    gyration_radius_2 = planet_gyration_radius_squared_venus
    planet_radius_of_gyration = np.sqrt(gyration_radius_2)

    Tab_size = 32
    Tab_w_lmpq = np.zeros((Tab_size,Tab_size))
    Tab_ImK2 = np.zeros((Tab_size,Tab_size))
    Tab_ReK2 = np.zeros((Tab_size,Tab_size))
    k = 0
    for i in range(0,int(size/Tab_size)):
        Tab_w_lmpq[i,:] = w_lmpq[i*Tab_size : (1+i)*Tab_size]
        Tab_ReK2[i,:] = ReK2[i*Tab_size : (1+i)*Tab_size]
        Tab_ImK2[i,:] = ImK2[i*Tab_size : (1+i)*Tab_size]
        k+=1
    for i in range (k*Tab_size, size):
        Tab_w_lmpq[k, i-k*Tab_size] = w_lmpq[i]
        Tab_ReK2[k, i-k*Tab_size] = ReK2[i]
        Tab_ImK2[k, i-k*Tab_size] = ImK2[i]

    freq_lmpq = []
    realK2 = []
    imagK2 = []
    freq_lmpq.append(Tab_w_lmpq)
    realK2.append(Tab_ReK2)
    imagK2.append(Tab_ImK2)
    planet_kaula_tidal_parameters_love_numbers = {
        "love_number_excitation_frequency": Tab_w_lmpq.tolist(), #syntaxe tab [[32]32]
        "imaginary_part_love_number": Tab_ImK2.tolist(),
        "real_part_love_number": Tab_ReK2.tolist(),
        "num_datapoints": size,
    }

    # --- Choose the tidal model to use :
    planet_tides_model = posidonius.effects.tides.ConstantTimeLag(planet_tides_parameters)
    # planet_tides_model = posidonius.effects.tides.Kaula(planet_kaula_tidal_parameters_love_numbers)

    # --- Choose the type of particle (central body or orbiting body) :
    #planet_tides = posidonius.effects.tides.CentralBody(planet_tides_model)
    planet_tides = posidonius.effects.tides.OrbitingBody(planet_tides_model)
    #planet_tides = posidonius.effects.tides.Disabled()
    #
    # planet_rotational_flattening_parameters = {"love_number": planet_tides_parameters["love_number"]}
    # planet_rotational_flattening_parameters = {"uniform_viscosity_coefficient": 1e14}
    # planet_rotational_flattening_model = posidonius.effects.rotational_flattening.OblateSpheroid(planet_rotational_flattening_parameters)
    # planet_rotational_flattening_model = posidonius.effects.rotational_flattening.CreepCoplanar(planet_rotational_flattening_parameters)
    #planet_rotational_flattening = posidonius.effects.rotational_flattening.CentralBody(planet_rotational_flattening_model)
    # planet_rotational_flattening = posidonius.effects.rotational_flattening.OrbitingBody(planet_rotational_flattening_model)
    planet_rotational_flattening = posidonius.effects.rotational_flattening.Disabled()
    #
    #planet_general_relativity = posidonius.effects.general_relativity.CentralBody("Kidder1995")
    #planet_general_relativity = posidonius.effects.general_relativity.CentralBody("Anderson1975")
    #planet_general_relativity = posidonius.effects.general_relativity.CentralBody("Newhall1983")
    # planet_general_relativity = posidonius.effects.general_relativity.OrbitingBody()
    planet_general_relativity = posidonius.effects.general_relativity.Disabled()
    #
    #planet_wind = posidonius.effects.wind.Interaction({
        ## Solar wind parametrisation (Bouvier 1997)
        #"k_factor": 4.0e-18, # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
        #"rotation_saturation": 1.7592918860102842, # 14. * TWO_PI/25.0, in units of the spin of the Sun today
    #})
    planet_wind = posidonius.effects.wind.Disabled()
    #
    #disk_surface_density_normalization_gcm = 1000. # g.cm^-2
    #disk_surface_density_normalization_SI = disk_surface_density_normalization_gcm * 1.0e-3 * 1.0e4 # kg.m^-2
    #disk_properties = {
        #'inner_edge_distance': 0.01,  # AU
        #'outer_edge_distance': 100.0, # AU
        #'lifetime': 1.0e5 * 365.25e0, # days
        #'alpha': 1.0e-2,
        #'surface_density_normalization': disk_surface_density_normalization_SI * (1.0/posidonius.constants.M_SUN) * posidonius.constants.AU**2, # Msun.AU^-2
        #'mean_molecular_weight': 2.4,
    #}
    #planet_disk = posidonius.effects.disk.CentralBody(disk_properties)
    #planet_disk = posidonius.effects.disk.OrbitingBody()
    planet_disk = posidonius.effects.disk.Disabled()
    #
    #planet_evolution = posidonius.GalletBolmont2017(planet_mass) # mass = 0.30 .. 1.40
    #planet_evolution = posidonius.BolmontMathis2016(planet_mass) # mass = 0.40 .. 1.40
    #planet_evolution = posidonius.Baraffe2015(planet_mass) # mass = 0.01 .. 1.40
    #planet_evolution = posidonius.Leconte2011(planet_mass) # mass = 0.01 .. 0.08
    #planet_evolution = posidonius.Baraffe1998(planet_mass) # Sun (mass = 1.0) or M-Dwarf (mass = 0.1)
    #planet_evolution = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
    #planet_evolution = posidonius.LeconteChabrier2013(True) # Jupiter with dissipation of dynamical tides
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

    whfast_alternative_coordinates="DemocraticHeliocentric"
    #whfast_alternative_coordinates="WHDS"
    #whfast_alternative_coordinates="Jacobi"
    #universe.write(filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)
    universe.write(filename, integrator="IAS15")
    #universe.write(filename, integrator="LeapFrog")


