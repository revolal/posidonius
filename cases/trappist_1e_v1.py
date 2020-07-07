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
    initial_time = 1.2e6*365.25 # time [days] where simulation starts
    time_step = 0.08 # 0.08 days
    #time_step = 0.05 # days
    #time_limit   = 4*time_step # days
    #time_limit   = 365.25 * 1.0e8 # days
    time_limit   =  365.25 * 1e6  #365.25 * 1e10 #365.25 * 1.0e6 # days
    historic_snapshot_period = 50.*365.25 #100.*365.25 # days
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
    # THE SUN https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
    star_mass = 0.089 # Solar masses
    # star_mass = 1.0
    star_radius_factor = 0.117
    # star_radius_factor = 1.0
    star_radius = star_radius_factor * posidonius.constants.R_SUN
    star_radius_of_gyration = 4.47e-01 # Brown dwarf
    # star_radius_of_gyration = 0.070 #The sun
    star_position = posidonius.Axes(0., 0., 0.)
    star_velocity = posidonius.Axes(0., 0., 0.)

    # Initialization of stellar spin
    star_rotation_period = 3.3*24 # hours HERE For T1-e
    # star_rotation_period = 609.12 # The sun in hours
    #star_rotation_period = 19.0*24 # hours
    star_angular_frequency = posidonius.constants.TWO_PI/(star_rotation_period/24.) # days^-1
    star_spin = posidonius.Axes(0., 0., star_angular_frequency)

    star_tides_parameters = {
        "dissipation_factor_scale": 0.0, #0.01,
        "dissipation_factor": 0.0, #2.006*3.845764e4,
        "love_number": 0.0, #0.307,
    }
    
    star_tides = posidonius.effects.tides.CentralBody(star_tides_parameters)
    #star_tides = posidonius.effects.tides.OrbitingBody(star_tides_parameters)
    # star_tides = posidonius.effects.tides.Disabled()
    #

    #star_rotational_flattening_parameters = {"love_number": star_tides_parameters["love_number"] }
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

    #############################################################################
    Fichiertxt = 'Results_Trappist_1e_00.0_138_freq_Imk2_posidonius.txt'
    # Fichiertxt = 'Results_Trappist_1e_05.0_150_freq_Imk2_posidonius.txt'
    # planet_data = np.loadtxt('Results_Trappist_1e_00.0_138_freq_Imk2_posidonius.txt',comments='#')
    # planet_data = np.loadtxt('Results_Trappist_1e_05.0_150_freq_Imk2_posidonius.txt',comments='#')

    planet_data = np.loadtxt(Fichiertxt,comments='#')
    print("Fichier ", Fichiertxt)
    
    planet_masses, planet_radiuses, planet_gyration_radius = planet_data[0,:]  #in kg, m and radius of gyration squared 
 
    w_lmpq = planet_data [1:,0]
    ImK2 = planet_data [1:,1]
    ReK2 = planet_data [1:,2]
    size = np.size(w_lmpq)

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

    ############################################################################
    # #Trappist 1e only:
    # #Radiuses in R_EARTH
    # planet_radiuses = 0.918e0
    # #Masses in M_SUN
    # planet_masses = 1.8850689e-06
    # #Semi-major axis in AU
    # planet_a = 0.02817
    # #Inclination in degrees
    # planet_i = 0.140000
    # #Mean anomaly in degrees
    # planet_l = 165.724187804


    # ############################################################################
    # https://nssdc.gsfc.nasa.gov/planetary/factsheet/
    #MERCURY
    # planet_radiuses = 4879. / 6378. #in R_EARTH
    # planet_i = 0.0
    # planet_inclination = 0.0
    # planet_obliquity = 0.0
    # e = 0.205
    # planet_a = 0.387 #AU
    # planet_masses = 0.330e24/posidonius.constants.M_SUN # Mearth in  Msol
    # planet_angular_frequency = posidonius.constants.TWO_PI/(10.) # 1. Days
    # planet_l = 165.724187804
    ##############################################################################
    # planet_mass = planet_masses /posidonius.constants.M_SUN # kg in Solar masses (3.0e-6 solar masses = 1 earth mass)
    # planet_radius_factor = planet_radiuses/63781137. # m in R Earth
    # planet_radius = planet_radius_factor * posidonius.constants.R_EARTH # R_earth in AU
    # planet_radius_of_gyration = planet_gyration_radius # 5.75e-01
    # #Semi-major axis in AU
    # planet_a = 0.02817
    # #Inclination in degrees
    # planet_i = 0.140000
    # #Mean anomaly in degrees
    # planet_l = 165.724187804
    # planet_obliquity = 0.0
    # a = planet_a
    # i = planet_i
    # l = planet_l
    #################################################################################
    planet_l = 0. #165.724187804
    planet_radiuses = 0.91e0 #0.4e0 ==========================================================================================PlanetRaduis
    planet_i = 0.0
    planet_inclination = 0.0
    planet_obliquity = 0.0
    e = 0. #0.0051 #0.07 #========================================================================================================ECC
    planet_a = 0.02928285# 0.02817 #0.015 #0.02817 #AU ==========================================================================SemiMaj 
    planet_masses = 0.772 * 5.9736E+24/1.9818e30 #1.898e27/1.989e30 #1.8850689e-06 #2.98e-6 # Mearth in  Msol ===================PLanetMass

    planet_angular_frequency = posidonius.constants.TWO_PI/(5./12.) #TEST 1. Days#(5/12) #10 hours periods in day ===============PlanetSpin
    # planet_angular_frequency = posidonius.constants.TWO_PI / 6.144379741742022 #TEST Pseudosyncrho in rad/Days================================================PlanetSpin
    print("Verif angular frequency",posidonius.constants.TWO_PI/planet_angular_frequency )

    planet_mass = planet_masses # m # Solar masses (3.0e-6 solar masses = 1 earth mass)
    planet_radius_factor = planet_radiuses # r # R Earth
    planet_radius = planet_radius_factor * posidonius.constants.R_EARTH
    planet_radius_of_gyration = np.sqrt( planet_gyration_radius ) # 5.75e-01
    a = planet_a
    i = planet_i
    l = planet_l
    # #Verification
    print("Planet mass in Msol", planet_masses)
    print("Planet radius in Rearth", planet_radiuses)
    # periode_orbital = posidonius.constants.TWO_PI * pow( planet_a*posidonius.constants.AU, 3./2.) / np.sqrt(posidonius.constants.G_SI*( star_mass )*posidonius.constants.M_SUN ) 
    # print("The Orbital period calculated", periode_orbital/posidonius.constants.DAY)
    # masse_planete = (pow(posidonius.constants.TWO_PI,2)*pow(planet_a*posidonius.constants.AU,3))/(posidonius.constants.G_SI* pow(6.09*posidonius.constants.DAY,2)) - star_mass*posidonius.constants.M_SUN
    # print("masse planet", masse_planete)
    # ############################################################################
    # # Radiuses in R_EARTH
    # planet_radiuses = (1.086e0, 1.056e0, 0.772e0, 0.918e0, 1.045e0, 1.127e0, 0.755e0)
    # # Masses in M_SUN
    # planet_masses = (2.5843686e-06, 4.1957984e-06, 1.2465778e-06, 1.8850689e-06, 2.0674949e-06, 4.0741811e-06, 1.2465778e-06)
    # # Semi-major axis in AU
    # planet_a = (0.01111, 0.01521, 0.02144, 0.02817, 0.0371, 0.0451, 0.0596)
    # # Inclination in degrees
    # planet_i = (0.3500000, 0.330000, 0.250000, 0.140000, 0.320000, 0.290000, 0.130000)
    # # Mean anomaly in degrees
    # planet_l = (323.732652895, 96.4925777097, 111.770368348, 165.724187804, 254.117367005, 161.020362506, 134.724813585)
    ###############################################################################
    #for r, m, a, i, l in zip(planet_radiuses, planet_masses, planet_a, planet_i, planet_l):
    # #HERE
    # planet_mass = planet_masses # m # Solar masses (3.0e-6 solar masses = 1 earth mass)
    # planet_radius_factor = planet_radiuses # r # R Earth
    # planet_radius = planet_radius_factor * posidonius.constants.R_EARTH
    # planet_radius_of_gyration = planet_gyration_radius # 5.75e-01
    # a = planet_a
    # i = planet_i
    # l = planet_l


    #////////// Specify initial position and velocity for a stable orbit
    #////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    #a = 0.01111;                             # semi-major axis (in AU)
    # e = 0.0000010;                              # eccentricity
    i = i * posidonius.constants.DEG2RAD;                      # inclination (degrees)
    p = 0. * posidonius.constants.DEG2RAD;                                # argument of pericentre (degrees)
    n = 0. * posidonius.constants.DEG2RAD;                      # longitude of the ascending node (degrees)
    l = l * posidonius.constants.DEG2RAD;           # mean anomaly (degrees)
    p = (p + n);                 # Convert to longitude of perihelion !!
    q = a * (1.0 - e);                     # perihelion distance
    planet_position, planet_velocity = posidonius.calculate_cartesian_coordinates(planet_mass, q, e, i, p, n, l, masses=[star_mass], positions=[star_position], velocities=[star_velocity])

    #////// Initialization of planetary spin
    # planet_obliquity = 1.0e-4 # rad
    # Pseudo-synchronization period
    # planet_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(a, e, star_mass, planet_mass) # days
    # planet_angular_frequency = posidonius.constants.TWO_PI/(planet_pseudo_synchronization_period) # days^-1
    planet_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
    # planet_inclination = planet_keplerian_orbital_elements[3]
    # planet_spin = posidonius.calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity)

    planet_spin = posidonius.calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity)

    k2pdelta = 2.465278e-3 # Terrestrial planets (no gas)

    planet_tides_parameters = {
        "dissipation_factor_scale": 1.0,
        "dissipation_factor": 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet_radius, 5)),
        "love_number": 0.299,
        # "love_number_excitation_frequency": 0.99,
        # "imaginary_part_love_number": 0.98,
        # "real_part_love_number": 0.97,
        # "num_datapoints": 0.96,
    }

    planet_tides_parameters_love_number_kaula = {
        "dissipation_factor_scale": 0.0, #1.0,
        "dissipation_factor": 0.0, #2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet_radius, 5)),
        "love_number": 0.0, #0.299,
        # "love_number_excitation_frequency": freq_lmpq,
        # "imaginary_part_love_number": imagK2,
        # "real_part_love_number": realK2,
        # "num_datapoints": size,
        "love_number_excitation_frequency": Tab_w_lmpq.tolist(), #syntaxe tab [[32]32]
        "imaginary_part_love_number": Tab_ImK2.tolist(),
        "real_part_love_number": Tab_ReK2.tolist(),
        "num_datapoints": size,
        # "love_number_excitation_frequency": 0.99,
        # "imaginary_part_love_number": 0.98,
        # "real_part_love_number": 0.97,
        # "num_datapoints": 0.96,
    }

    # planet_tides = posidonius.effects.tides.CentralBody(planet_tides_parameters)
    # planet_tides = posidonius.effects.tides.OrbitingBody(planet_tides_parameters)
    planet_tides = posidonius.effects.tides.KaulaCoplanarOrbitingBody(planet_tides_parameters_love_number_kaula)    
    #planet_tides = posidonius.effects.tides.Disabled()

    #planet_tides = posidonius.effects.tides.ConstTimeLagCentralBody(planet_tides_parameters)
    #planet_tides = posidonius.effects.tides.ConstTimeLagOrbitingBody(planet_tides_parameters)


    planet_rotational_flattening_parameters = {"love_number": 0.9532}
    #planet_rotational_flattening = posidonius.effects.rotational_flattening.CentralBody(planet_rotational_flattening_parameters)
    #planet_rotational_flattening = posidonius.effects.rotational_flattening.OrbitingBody(planet_rotational_flattening_parameters)
    planet_rotational_flattening = posidonius.effects.rotational_flattening.Disabled()
    #
    #planet_general_relativity = posidonius.effects.general_relativity.CentralBody("Kidder1995")
    #planet_general_relativity = posidonius.effects.general_relativity.CentralBody("Anderson1975")
    #planet_general_relativity = posidonius.effects.general_relativity.CentralBody("Newhall1983")
    #planet_general_relativity = posidonius.effects.general_relativity.OrbitingBody()
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
    # universe.write(filename, integrator="IAS15")
    #universe.write(filename, integrator="LeapFrog")


