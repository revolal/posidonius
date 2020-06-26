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
    initial_time = 1.0e6*365.25 # time [days] where simulation starts
    time_step = 0.01 # days
    #time_step = 0.05 # days
    #time_limit   = 4*time_step # days
    #time_limit   = 365.25 * 1.0e8 # days
    time_limit   = 1.#365.25 * 1.0e3 # days
    historic_snapshot_period = 0.01#*365.25 # days
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

    star_mass = 0.089 # Solar masses
    star_radius_factor = 0.117
    star_radius = star_radius_factor * posidonius.constants.R_SUN
    star_radius_of_gyration = 4.47e-01 # Brown dwarf
    star_position = posidonius.Axes(0., 0., 0.)
    star_velocity = posidonius.Axes(0., 0., 0.)

    # Initialization of stellar spin
    star_rotation_period = 3.3*24 # hours
    #star_rotation_period = 19.0*24 # hours
    star_angular_frequency = posidonius.constants.TWO_PI/(star_rotation_period/24.) # days^-1
    star_spin = posidonius.Axes(0., 0., star_angular_frequency)

    # star_tides_parameters = {
    #     "dissipation_factor_scale": 0.01,
    #     "dissipation_factor": 2.006*3.845764e4,
    #     "love_number": 0.307,
    # }
    star_tides_parameters = {
        "dissipation_factor_scale": 0.0,
        "dissipation_factor": 0.0,
        "love_number": 0.0,
    }
    star_tides = posidonius.effects.tides.CentralBody(star_tides_parameters)
    #star_tides = posidonius.effects.tides.OrbitingBody(star_tides_parameters)
    #star_tides = posidonius.effects.tides.Disabled()
    #
    star_rotational_flattening_parameters = {"love_number": star_tides_parameters["love_number"] }
    #star_rotational_flattening = posidonius.effects.rotational_flattening.CentralBody(star_rotational_flattening_parameters)
    #star_rotational_flattening = posidonius.effects.rotational_flattening.OrbitingBody(star_rotational_flattening_parameters)
    star_rotational_flattening = posidonius.effects.rotational_flattening.Disabled()
    #
    #star_general_relativity = posidonius.effects.general_relativity.CentralBody("Kidder1995")
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

    #==================================================================================================================
    print("Recover K2 data")
    

    #===============================================================================================
    planet_data_b = np.loadtxt('Results_Trappist_1b_20.5_100_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_b, planet_radius_b, planet_gyration_radius_squared_b = planet_data_b[0,:] 

    w_lm_b   = planet_data_b[1:,0]
    ImK2_b   = planet_data_b[1:,1]
    ReK2_b   = planet_data_b[1:,2]
    size_b   = np.size(w_lm_b)

    #===============================================================================================
    planet_data_c = np.loadtxt('Results_Trappist_1c_07.5_100_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_c, planet_radius_c, planet_gyration_radius_squared_c = planet_data_c[0,:] 

    w_lm_c   = planet_data_b[1:,0]
    ImK2_c   = planet_data_b[1:,1]
    ReK2_c   = planet_data_b[1:,2]
    size_c   = np.size(w_lm_b)

    #===============================================================================================
    planet_data_d = np.loadtxt('Results_Trappist_1d_19.5_100_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_d, planet_radius_d, planet_gyration_radius_squared_d = planet_data_d[0,:] 

    w_lm_d   = planet_data_b[1:,0]
    ImK2_d   = planet_data_b[1:,1]
    ReK2_d   = planet_data_b[1:,2]
    size_d   = np.size(w_lm_b)

    #===============================================================================================
    planet_data_e = np.loadtxt('Results_Trappist_1e_05.0_150_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_e, planet_radius_e, planet_gyration_radius_squared_e = planet_data_e[0,:] 

    w_lm_e   = planet_data_b[1:,0]
    ImK2_e   = planet_data_b[1:,1]
    ReK2_e   = planet_data_b[1:,2]
    size_e   = np.size(w_lm_b)

    #===============================================================================================
    planet_data_f = np.loadtxt('Results_Trappist_1f_14.5_150_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_f, planet_radius_f, planet_gyration_radius_squared_f = planet_data_f[0,:] 

    w_lm_f   = planet_data_b[1:,0]
    ImK2_f   = planet_data_b[1:,1]
    ReK2_f   = planet_data_b[1:,2]
    size_f   = np.size(w_lm_b)

    #===============================================================================================
    planet_data_g = np.loadtxt('Results_Trappist_1g_22.5_150_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_g, planet_radius_g, planet_gyration_radius_squared_g = planet_data_g[0,:] 

    w_lm_g   = planet_data_b[1:,0]
    ImK2_g   = planet_data_b[1:,1]
    ReK2_g   = planet_data_b[1:,2]
    size_g   = np.size(w_lm_b)

    #===============================================================================================
    planet_data_h = np.loadtxt('Results_Trappist_1h_13.5_150_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_h, planet_radius_h, planet_gyration_radius_squared_h = planet_data_h[0,:] 

    w_lm_h   = planet_data_b[1:,0]
    ImK2_h   = planet_data_b[1:,1]
    ReK2_h   = planet_data_b[1:,2]
    size_h   = np.size(w_lm_b)


    ############################################################################
    # Radiuses in R_EARTH
    # planet_radiuses = (1.086e0, 1.056e0, 0.772e0, 0.918e0, 1.045e0, 1.127e0, 0.755e0)
    planet_radiuses = (planet_radius_b, planet_radius_c, planet_radius_d, planet_radius_e, planet_radius_f, planet_radius_g, planet_radius_h)
    # Masses in M_SUN
    # planet_masses = (2.5843686e-06, 4.1957984e-06, 1.2465778e-06, 1.8850689e-06, 2.0674949e-06, 4.0741811e-06, 1.2465778e-06)
    planet_masses = (planet_mass_b, planet_mass_c, planet_mass_d, planet_mass_e, planet_mass_f, planet_mass_g, planet_mass_h)
    # Semi-major axis in AU
    # planet_a = (0.01111, 0.01521, 0.02144, 0.02817, 0.0371, 0.0451, 0.0596)
    planet_a = (0.01154775, 0.01581512, 0.02228038, 0.02928285, 0.03853361, 0.04687692, 0.06193488)
    # Inclination in degrees
    # planet_i = (0.3500000, 0.330000, 0.250000, 0.140000, 0.320000, 0.290000, 0.130000)
    planet_i = (0., 0., 0., 0., 0., 0., 0.)
    # Mean anomaly in degrees
    # planet_l = (323.732652895, 96.4925777097, 111.770368348, 165.724187804, 254.117367005, 161.020362506, 134.724813585)
    planet_l = (0., 0., 0., 0., 0., 0., 0.)

    planet_e = (0., 0., 0., 0., 0., 0., 0.)

    planet_w_lm = (w_lm_b, w_lm_c, w_lm_d, w_lm_e, w_lm_f, w_lm_g, w_lm_h)

    planet_Imk2 = (ImK2_b, ImK2_c, ImK2_d, ImK2_e, ImK2_f, ImK2_g, ImK2_h)

    planet_Rek2 = (ReK2_b, ReK2_c, ReK2_d, ReK2_e, ReK2_f, ReK2_g, ReK2_h)

    planet_nm_data = (size_b, size_c, size_d, size_e, size_f, size_g, size_h)

    planet_gyration_radius_squared = (planet_gyration_radius_squared_b, planet_gyration_radius_squared_c, planet_gyration_radius_squared_d, planet_gyration_radius_squared_e, planet_gyration_radius_squared_f, planet_gyration_radius_squared_g, planet_gyration_radius_squared_h)

    for r, m, a, i, l, e, w_lmpq, ImK2, ReK2, size, gyration_radius_2 in zip(planet_radiuses, planet_masses, planet_a, planet_i, planet_l, planet_e, planet_w_lm, planet_Imk2, planet_Rek2, planet_nm_data, planet_gyration_radius_squared):
        planet_mass = m / posidonius.constants.M_SUN # Solar masses (3.0e-6 solar masses = 1 earth mass)
        planet_radius_factor = r # R Earth
        planet_radius = planet_radius_factor * posidonius.constants.R_EARTH / 6.3781e6
        # planet_radius_of_gyration = 5.75e-01
        planet_radius_of_gyration = np.sqrt(gyration_radius_2)

        #////////// Specify initial position and velocity for a stable orbit
        #////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
        # a = 0.01111;                             # semi-major axis (in AU)
        # e = #0000010;                              # eccentricity
        i = i * posidonius.constants.DEG2RAD;                      # inclination (degrees)
        p = 0. * posidonius.constants.DEG2RAD;                                # argument of pericentre (degrees)
        n = 0. * posidonius.constants.DEG2RAD;                      # longitude of the ascending node (degrees)
        l = l * posidonius.constants.DEG2RAD;           # mean anomaly (degrees)
        p = (p + n);                 # Convert to longitude of perihelion !!
        q = a * (1.0 - e);                     # perihelion distance
        planet_position, planet_velocity = posidonius.calculate_cartesian_coordinates(planet_mass, q, e, i, p, n, l, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
        # print("PLANET POSITION VELOCITY ",planet_position.get(), planet_velocity.get())

        #////// Initialization of planetary spin
        planet_obliquity = 0. #1.0e-4 # rad
        # Pseudo-synchronization period
        planet_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(a, e, star_mass, planet_mass) # days
        # planet_angular_frequency = posidonius.constants.TWO_PI/(planet_pseudo_synchronization_period) # days^-1
        planet_angular_frequency = posidonius.constants.TWO_PI/(5./12.) # days^-1
        planet_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
        planet_inclination = planet_keplerian_orbital_elements[3]
        planet_spin = posidonius.calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity)

        k2pdelta = 2.465278e-3 # Terrestrial planets (no gas)
        planet_tides_parameters = {
            "dissipation_factor_scale": 1.0,
            "dissipation_factor": 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet_radius, 5)),
            "love_number": 0.299,
        }

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


        planet_tides_parameters_love_number_kaula = {
            "dissipation_factor_scale": 0.0, #1.0,
            "dissipation_factor": 0.0, #2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet_radius, 5)),
            "love_number": 0.0, #0.299,

            "love_number_excitation_frequency": Tab_w_lmpq.tolist(), #syntaxe tab [[32]32]
            "imaginary_part_love_number": Tab_ImK2.tolist(),
            "real_part_love_number": Tab_ReK2.tolist(),
            "num_datapoints": size,
        }


        #planet_tides = posidonius.effects.tides.CentralBody(planet_tides_parameters)
        #planet_tides = posidonius.effects.tides.OrbitingBody(planet_tides_parameters)
        # planet_tides = posidonius.effects.tides.ConstTimeLagOrbitingBody(planet_tides_parameters)
        planet_tides = posidonius.effects.tides.KaulaCoplanarOrbitingBody(planet_tides_parameters_love_number_kaula)
        #planet_tides = posidonius.effects.tides.Disabled()
        #
        planet_rotational_flattening_parameters = {"love_number": 0.9532}
        #planet_rotational_flattening = posidonius.effects.rotational_flattening.CentralBody(planet_rotational_flattening_parameters)
        # planet_rotational_flattening = posidonius.effects.rotational_flattening.OrbitingBody(planet_rotational_flattening_parameters)
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


    ############################################################################

    whfast_alternative_coordinates="DemocraticHeliocentric"
    #whfast_alternative_coordinates="WHDS"
    #whfast_alternative_coordinates="Jacobi"
    universe.write(filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)
    #universe.write(filename, integrator="IAS15")
    #universe.write(filename, integrator="LeapFrog")


