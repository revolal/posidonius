import posidonius
import numpy as np
import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('output_filename', action='store', help='Filename where the initial snapshot will be stored (e.g., universe_integrator.json)')

    args = parser.parse_args()
    filename = args.output_filename
    #filename = posidonius.constants.BASE_DIR+"target/case7.json"

    #initial_time = 4.5e6*365.25 # time [days] where simulation starts
    initial_time = 1.0*365.25 # time [days] where simulation starts
    time_step = 0.02 # days
    #time_step = 0.07 # days
    #time_limit   = 4*time_step # days
    #time_limit   = 365.25 * 1.0e8 # days
    time_limit   = 2.*365.25 #* 1.0e2 # days
    historic_snapshot_period = 1. #0.01*365.25 # days
    recovery_snapshot_period = historic_snapshot_period # days
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
    star_tides_model = posidonius.effects.tides.ConstantTimeLag(star_tides_parameters)
    star_tides = posidonius.effects.tides.CentralBody(star_tides_model)
    #star_tides = posidonius.effects.tides.OrbitingBody(star_tides_model)
    # star_tides = posidonius.effects.tides.Disabled()
    #
    # star_rotational_flattening_parameters = {"love_number": star_tides_parameters["love_number"] }
    #star_rotational_flattening_model = posidonius.effects.rotational_flattening.OblateSpheroid(star_rotational_flattening_parameters)
    #star_rotational_flattening = posidonius.effects.rotational_flattening.CentralBody(star_rotational_flattening_model)
    #star_rotational_flattening = posidonius.effects.rotational_flattening.OrbitingBody(star_rotational_flattening_model)
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

    # === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === ===
    # dir_path = './cases/Trappist_1_BurnMan/'
    dir_path = './cases/love_number_data/k2_trappist-1_earth/'
    res = []
    planet_param = []
    planet_masses = []
    planet_radiuses = []
    planet_gyration_radius_squared = []
    planet_w_lm = []
    planet_Imk2 = []
    planet_Rek2 = []
    planet_nm_data = []

    for path in os.listdir(dir_path):
        if os.path.isfile(os.path.join(dir_path, path)):
            res.append(path)
    res.sort()
    for path in res:
        print(path)
        planet_data = np.loadtxt(dir_path+path,comments='#')
        mass, radius, gyration_radius_squared = planet_data[0,:] 
        planet_masses.append(mass)
        planet_radiuses.append(radius)
        planet_gyration_radius_squared.append(gyration_radius_squared)
        planet_param.append(planet_data[0,:])
        planet_w_lm.append(planet_data[1:,0])
        planet_Imk2.append(planet_data[1:,1])
        planet_Rek2.append(planet_data[1:,2])
        planet_nm_data.append(np.size(planet_data[1:,0]))
    # res.sort()
    print(res)


    ############################################################################

     # --- Planets periods in  days
    planet_period =(1.5112868036098193, 2.4227038555269917, 4.0504368087076257, 6.1029264019452052, 9.2103976604134150, 12.3562321332659266, 18.7787920964476669) # in days From AS Libert
    # --- Planets masses
    planet_masses = (4.1297145833e-06, 3.9564821698e-06, 1.1720208671e-06, 2.1148831802e-06, 3.1362928581e-06, 3.9669131016e-06, 9.8160839680e-07) # in Msun From AS Libert
    tmp = [ mm *posidonius.constants.M_SUN  for mm in planet_masses]
    print(f'From AS Libert :: \t {tmp}') # Avoid the effect of the masses on the dynamisc --> only take the k2 from the files
    # --- Semi-major axis in AU
    planet_a = ( 0.0115418669425055, 0.0158094975622961, 0.0222700389997286, 0.0292694953252985, 0.0385106150958128, 0.0468445028220842, 0.0619224236258543 ) # in AU From AS Libert
    # --- Inclination in degrees
    planet_i = (0., 0., 0., 0., 0., 0., 14.)
    # --- Eccentricity
    planet_e = ( 0.0066945347902135, 0.0019585600861175, 0.0072217598426537, 0.0052151102407526, 0.0095027023237257, 0.0032634092345139, 0.0044037179966920) # From AS Libert 
    # --- Mean longitude 
    planet_lamdba = ( 6.2807871243880449, 3.2752103444471459, 0.2445003452092198, 1.4578342582705632, 1.0048074576831532, 1.4551190036825181, 5.0763862710529217)  # in rad From AS Libert 
    # --- Argument of pericentre
    # planet_p = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # --- longitude of ascending node
    planet_n = (0., 0., 0., 0., 0., 0., 20.)
    # --- longitude of pericentre = argument of pericenter + longitude ascending node
    planet_p = ( 3.9138033639162217, 3.9568810173922708, 3.6462719004641233, 0.9017217708596807, 3.1105825606949371, 5.8092293524825127, 3.2359374491128419) # in rad From AS Libert 
    # planet_p = (0., 0., 0., 0., 0., 0., 0.)
    # --- Mean anomaly in degrees
    planet_mean_ano =  [ mean_lon - long_peri for (mean_lon, long_peri) in zip( planet_lamdba, planet_p) ] # compute from Mean longitude and longitude of pericentre
    print(f' -- Recomputed Mean anomaly {planet_mean_ano}')

    planets = ['T1-b','T1-c', 'T1-d', 'T1-e', 'T1-f', 'T1-g', 'T1-h']

    key=6

    (planet, r, m, a, i, l, e, p, n, w_lmpq, ImK2, ReK2, size, gyration_radius_2) = (planets[key], planet_radiuses[key], planet_masses[key], planet_a[key], planet_i[key], planet_mean_ano[key], planet_e[key], planet_p[key], planet_n[key], planet_w_lm[key], planet_Imk2[key], planet_Rek2[key], planet_nm_data[key], planet_gyration_radius_squared[key])
    planet_mass = m #/ posidonius.constants.M_SUN # Solar masses (3.0e-6 solar masses = 1 earth mass)
    # planet_radius_factor = r # R Earth
    planet_radius = r / posidonius.constants.AU # planet_radius_factor * posidonius.constants.R_EARTH
    planet_radius_of_gyration = np.sqrt(gyration_radius_2) #5.75e-01

    # e = 0.16

    # --- Specify initial position and velocity for a stable orbit
    # -- Keplerian orbital elements, in the `asteroidal' format of Mercury code
    i = i * posidonius.constants.DEG2RAD;   # inclination (degrees)
    # p = p * posidonius.constants.DEG2RAD;   # argument of pericentre (degrees)
    n = n * posidonius.constants.DEG2RAD;  # longitude of the ascending node (degrees)
    l = l * posidonius.constants.DEG2RAD;   # mean anomaly (degrees)
    # p = (p + n);# Convert to longitude of perihelion !!
    q = a * (1.0 - e);# perihelion distance
    planet_position, planet_velocity = posidonius.calculate_cartesian_coordinates(planet_mass, q, e, i, p, n, l, masses=[star_mass], positions=[star_position], velocities=[star_velocity])

    # --- Initialization of planetary spin
    planet_obliquity = 0.0 # 1.0e-4 # rad
    # --- Initialize spin
    orbital_period = 2.*np.pi*(np.sqrt( (a*posidonius.constants.AU)**3. / ( posidonius.constants.G_SI*( star_mass*posidonius.constants.M_SUN + planet_mass*posidonius.constants.M_SUN ) ) ) )/86400.
    planet_rotation_period = orbital_period/(1.0) #24. # hours orbital_period/(1.0 + 1e-5) #24. # hours
    planet_angular_frequency = posidonius.constants.TWO_PI/(planet_rotation_period) # days^-1
    # --- Pseudo-synchronization period
    # planet_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(a, e, star_mass, planet_mass) # days
    # planet_angular_frequency = posidonius.constants.TWO_PI/(planet_pseudo_synchronization_period) # days^-1
    # ---
    # planet_keplerian_orbital_elements = posidonius.tools.calculate_keplerian_orbital_elements_bis(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
    planet_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
    planet_inclination = planet_keplerian_orbital_elements[3]
    planet_spin = posidonius.calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity)

    print("\n|--- Check "+planet+' : ')
    print("| \t Mass {} Radius {} Rg^2 {} ".format(m,r,gyration_radius_2))
    print("| \t-- Initial orbital elements : ")
    print("| \t sma {}\n| \t ecc {}\n| \t incl {}\n| \t long_asc {}\n| \t long_peri {}\n| \t mean_ano {}\n| \t mean_long {} ".format(a,e,i,n%(2.*np.pi),p%(2.*np.pi),l%(2.*np.pi),(l+p)%(2.*np.pi)))
    (sma, _q, ecc, incl, long_peri, long_asc, mean_ano, _true_ano) = planet_keplerian_orbital_elements
    print("| \t-- Computed orbital elements : ")
    print("| \t sma {}\n| \t ecc {}\n| \t incl {}\n| \t long_asc {}\n| \t long_peri {}\n| \t mean_ano {}\n| \t mean_long {} ".format(sma,ecc,incl, long_asc%(2.*np.pi), long_peri%(2.*np.pi), mean_ano%(2.*np.pi), (long_peri+mean_ano)%(2.*np.pi) ))

    # k2pdelta = 2.465278e-3 # Terrestrial planets (no gas)
    # planet_tides_parameters = {
    #     "dissipation_factor_scale": 1.0,
    #     "dissipation_factor": 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet_radius, 5)),
    #     "love_number": 0.299,
    # }

    # --- Kaula parameters
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

    k2pdelta = 2.465278e-3 # Terrestrial planets (no gas)
    planet_tides_parameters = {
        "dissipation_factor_scale": 1.0,
        "dissipation_factor": 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet_radius, 5)),
        "love_number": 0.299,
    }

    # --- Choose the tidal model to use :
    planet_tides_model = posidonius.effects.tides.Kaula(planet_kaula_tidal_parameters_love_numbers)
    # planet_tides_model = posidonius.effects.tides.ConstantTimeLag(planet_tides_parameters)

    # --- Choose the type of particle (central body or orbiting body) :
    #planet_tides = posidonius.effects.tides.CentralBody(planet_tides_model)
    planet_tides = posidonius.effects.tides.OrbitingBody(planet_tides_model)
    # planet_tides = posidonius.effects.tides.Disabled()
    # ---
    #planet_rotational_flattening_parameters = {"love_number": 0.9532}
    #planet_rotational_flattening_model = posidonius.effects.rotational_flattening.OblateSpheroid(planet_rotational_flattening_parameters)
    #planet_rotational_flattening = posidonius.effects.rotational_flattening.CentralBody(planet_rotational_flattening_model)
    # planet_rotational_flattening = posidonius.effects.rotational_flattening.OrbitingBody(planet_rotational_flattening_model)
    planet_rotational_flattening = posidonius.effects.rotational_flattening.Disabled()
    # ---
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
    # universe.write(filename, integrator="LeapFrog")


