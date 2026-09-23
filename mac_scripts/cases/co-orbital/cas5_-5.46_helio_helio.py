import posidonius
import numpy as np
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('output_filename', action='store', help='Filename where the initial snapshot will be stored (e.g., universe_integrator.json)')

    args = parser.parse_args()
    filename = args.output_filename

    initial_time = 4.5e6*365.25 # time [days] where simulation starts
    time_step = 0.01 # days
    #time_limit   = 4*time_step # days
    time_limit   =  365.25 * 1.0e6 # days
    historic_snapshot_period = 365.25 * 10. #*time_step # days
    recovery_snapshot_period = 100.*historic_snapshot_period # days
    consider_effects = posidonius.ConsiderEffects({
        "tides": True,
        "rotational_flattening": False,
        "general_relativity": False,
        "disk": False,
        "wind": False,
        "evolution": False,
    })
    universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects)

    star_mass = 1.#0.089 # Solar masses
    star_radius_factor = 1.#0.117
    star_radius = star_radius_factor * posidonius.constants.R_SUN
    star_radius_of_gyration = 4.47e-01 # Brown dwarf
    star_position = posidonius.Axes(0., 0., 0.)
    star_velocity = posidonius.Axes(0., 0., 0.)

    star_position = posidonius.Axes(0., 0., 0.)
    star_velocity = posidonius.Axes(0., 0., 0.)

    # Initialization of stellar spin
    star_rotation_period = 3.3*24 # hours
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
    # --- Love number planet data
    #-----------------------------------------------------------------------------------------------
    planet_data_b = np.loadtxt('Results_Trappist_1b_20.5_100_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_b, planet_radius_b, planet_gyration_radius_squared_b = planet_data_b[0,:] 
    w_lm_b   = planet_data_b[1:,0]
    ImK2_b   = planet_data_b[1:,1]
    ReK2_b   = planet_data_b[1:,2]
    size_b   = np.size(w_lm_b)

    #-----------------------------------------------------------------------------------------------
    planet_data_c = np.loadtxt('Results_Trappist_1c_07.5_100_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_c, planet_radius_c, planet_gyration_radius_squared_c = planet_data_c[0,:] 
    w_lm_c   = planet_data_c[1:,0]
    ImK2_c   = planet_data_c[1:,1]
    ReK2_c   = planet_data_c[1:,2]
    size_c   = np.size(w_lm_c)

    #-----------------------------------------------------------------------------------------------
    planet_data_d = np.loadtxt('Results_Trappist_1d_19.5_100_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_d, planet_radius_d, planet_gyration_radius_squared_d = planet_data_d[0,:] 
    w_lm_d   = planet_data_d[1:,0]
    ImK2_d   = planet_data_d[1:,1]
    ReK2_d   = planet_data_d[1:,2]
    size_d   = np.size(w_lm_d)

    planet_w_lm = (w_lm_b, w_lm_c, w_lm_d )
    planet_Imk2 = (ImK2_b, ImK2_c, ImK2_d)
    planet_Rek2 = (ReK2_b, ReK2_c, ReK2_d)
    planet_nm_data = (size_b, size_c, size_d)
    planet_gyration_radius_squared = (planet_gyration_radius_squared_b, planet_gyration_radius_squared_c, planet_gyration_radius_squared_d)

    #-----------------------------------------------------------------------------------------------
    # --- FILL HERE ::
    # In the form :
    # data = (planet1, planet2, planet3)

    # Radiuses in m
    planet_radiuses = (0., 0., 0.)

    # Masses in kg
    m0=1.9891e30
    planet_masses = (1.8181818181818182e-05, 1.8181818181818182e-04, 1.0e-4)

    # Semi-major axis in AU
    planet_a = (0.0198874, 0.0198895, 0.0323435)

    # Eccentricity
    planet_e = (0.00183541, 0.00209126, 0.000457054)
    
    # Inclination in degrees
    planet_i = (0., 0., 0.)

    # Mean anomaly in degrees
    planet_mean_ano = (118.975+0.1, -0.0336723, -166.551)

    # argument of pericentre in degree
    planet_arg_peri = (-59.0032, 0.025452, 166.605)

    # longitude of ascending node in degree
    planet_long_asc = (0., 0., 0.)

    # ---
    #-----------------------------------------------------------------------------------------------


    for radius, mass, sma, ecc, inclin, mean_ano, arg_peri, long_asc , w_lmpq, ImK2, ReK2, size, gyration_radius_2 in zip \
        (planet_radiuses, planet_masses, planet_a, planet_e, planet_i, planet_mean_ano, planet_arg_peri, planet_long_asc, planet_w_lm, planet_Imk2, planet_Rek2, planet_nm_data, planet_gyration_radius_squared):
        
        planet_radius = radius / posidonius.constants.AU
        planet_mass = mass / posidonius.constants.M_SUN # Solar masses (3.0e-6 solar masses = 1 earth mass)
        
        # planet_radius_of_gyration = 5.75e-01
        planet_radius_of_gyration = np.sqrt(gyration_radius_2)

        # --- Specify initial position and velocity for a stable orbit
        # Keplerian orbital elements, in the `asteroidal' format of Mercury code
        a = sma;                                        # semi-major axis (in AU)
        e = ecc                                         # eccentricity
        i = inclin * posidonius.constants.DEG2RAD;      # inclination (degrees)
        l = mean_ano * posidonius.constants.DEG2RAD;    # mean anomaly (degrees)
        p = arg_peri * posidonius.constants.DEG2RAD;    # argument of pericentre (degrees)
        n = long_asc * posidonius.constants.DEG2RAD;    # longitude of the ascending node (degrees)
        
        p = (p + n);                            # Convert to longitude of perihelion !!
        q = a * (1.0 - e);                      # perihelion distance
        planet_position, planet_velocity = posidonius.calculate_cartesian_coordinates(planet_mass, q, e, i, p, n, l, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
        # print("PLANET POSITION VELOCITY ",planet_position.get(), planet_velocity.get())

        # --- Initialization of planetary spin
        # spin obliquity
        planet_obliquity = 0. # in rad
        # --- rotation period --- #

        # --- synchronization period
        orbital_period = 2.*np.pi*(np.sqrt( (a*posidonius.constants.AU)**3. / ( posidonius.constants.G_SI*( star_mass*posidonius.constants.M_SUN + planet_mass*posidonius.constants.M_SUN ) ) ) )/86400.
        planet_rotation_period = orbital_period/(1.0) 
        planet_angular_frequency = posidonius.constants.TWO_PI/(planet_rotation_period) # days^-1

        # --- Pseudo-synchronization period
        # planet_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(a, e, star_mass, planet_mass) # days
        # planet_angular_frequency = posidonius.constants.TWO_PI/(planet_pseudo_synchronization_period) # days^-1

        # --- in days
        # day_lenght = 5./12. # 10hrs rotation period
        # planet_angular_frequency = posidonius.constants.TWO_PI/day_lenght # days^-1

        planet_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
        planet_inclination = planet_keplerian_orbital_elements[3]
        planet_spin = posidonius.calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity)
        # ---

        # --- CTL's parameters --- #
        # k2pdelta =  2.465278e-3 #6.7497e-6 #2.465278e-3 # Terrestrial planets (no gas)
        # planet_tides_parameters = {
        #     "dissipation_factor_scale": 1.0,
        #     "dissipation_factor": 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet_radius, 5)),
        #     "love_number": 0.299,
        # }

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
        # planet_tides_model = posidonius.effects.tides.ConstantTimeLag(planet_tides_parameters)
        planet_tides_model = posidonius.effects.tides.Kaula(planet_kaula_tidal_parameters_love_numbers)

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
