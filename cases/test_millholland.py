import posidonius
import numpy as np
import argparse
import os, sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('output_filename', action='store', help='Filename where the initial snapshot will be stored (e.g., universe_integrator.json)')

    args = parser.parse_args()
    filename = args.output_filename
    #filename = posidonius.constants.BASE_DIR+"target/case7.json"

    #initial_time = 4.5e6*365.25 # time [days] where simulation starts
    initial_time = 1.0*365.25 # time [days] where simulation starts
    time_step = 0.3 # days
    #time_limit   = 4*time_step # days
    time_limit = 3e6 *365.25 # 1e4 *365.25 #* 1.0e10 # days
    historic_snapshot_period = 5e0*365.25  # days
    recovery_snapshot_period = 1e4*365.25  # days
    consider_effects = posidonius.ConsiderEffects({
        "tides": True,
        "rotational_flattening": False, 
        "general_relativity": False,
        "disk": True,
        "wind": False,
        "evolution": False,
    })
    universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects)

    star_mass = 0.0898 # Solar masses Agol+21
    star_radius_factor = 0.119 #0.117 # Agol+21
    star_radius = star_radius_factor * posidonius.constants.R_SUN
    star_radius_of_gyration = 4.47e-01 # Brown dwarf
    star_position = posidonius.Axes(0., 0., 0.)
    star_velocity = posidonius.Axes(0., 0., 0.)

    # Initialization of stellar spin
    # star_rotation_period = 3.3*24 # hours
    #star_rotation_period = 19.0*24 # hours
    # star_angular_frequency = posidonius.constants.TWO_PI/(star_rotation_period/24.) # days^-1
    star_rotation_period = 7.5 # days
    star_angular_frequency = posidonius.constants.TWO_PI/(star_rotation_period) # days^-1
    star_spin = posidonius.Axes(0., 0., star_angular_frequency)

    star_tides_parameters = {
        "dissipation_factor_scale": 0.01,
        "dissipation_factor": 2.006*3.845764e4,
        "love_number": 0.307,
    }
    # star_tides_parameters = {
    #     "dissipation_factor_scale": 0.0,
    #     "dissipation_factor": 0.0,
    #     "love_number": 0.0,
    # }
    star_tides_model = posidonius.effects.tides.ConstantTimeLag(star_tides_parameters)
    star_tides = posidonius.effects.tides.CentralBody(star_tides_model)
    # star_tides = posidonius.effects.tides.OrbitingBody(star_tides_model)
    # star_tides = posidonius.effects.tides.Disabled()
    #
    star_rotational_flattening_parameters = {"love_number": star_tides_parameters["love_number"] }
    # star_rotational_flattening_model = posidonius.effects.rotational_flattening.OblateSpheroid(star_rotational_flattening_parameters)
    # star_rotational_flattening = posidonius.effects.rotational_flattening.CentralBody(star_rotational_flattening_model)
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
    disk_surface_density_normalization_gcm = 10. # g.cm^-2
    disk_surface_density_normalization_SI = disk_surface_density_normalization_gcm * 1.0e-3 * 1.0e4 # kg.m^-2
    disk_properties = {
        'inner_edge_distance': 0.03,  # AU 0.0292694953252985, 0.0385106150958128
        'outer_edge_distance': 0.04, # AU
        'lifetime': 1e6 * 365.25e0, # days
        'alpha': 1.0e-2,
        'surface_density_normalization': disk_surface_density_normalization_SI * (1.0/posidonius.constants.M_SUN) * posidonius.constants.AU**2, # Msun.AU^-2
        'mean_molecular_weight': 2.4,
    }
    star_disk = posidonius.effects.disk.CentralBody(disk_properties)
    #star_disk = posidonius.effects.disk.OrbitingBody()
    # star_disk = posidonius.effects.disk.Disabled()
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

    # === PLANET PARAMETERS ===


    # --- Trappist1 Earth planets profiles
    print("  > Earth frequ selection")
    # dir_path = './input/love_numbers/trappist_love_numbers/Result_Trappist1_Earth/'
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
    res.sort()

    
    ############################################################################
    print(f'From Burnman :: \t {planet_masses} ')

    # --- Planets periods in  days
    # planet_period = (1.510826, 2.421937, 4.049219, 6.101013, 9.207540, 12.352446, 18.772866) #From Agol+21 Table 2 in days
    # planet_period =(1.5112868036098193, 2.4227038555269917, 4.0504368087076257, 6.1029264019452052, 9.2103976604134150, 12.3562321332659266, 18.7787920964476669) # in days From AS Libert
    planet_period =( 6.1029264019452052, 9.2103976604134150) # in days From AS Libert
    # planet_period =( 6.1029264019452052, 9.154389602917808) # 2nd planet tuned for 3:2

    # --- Planets masses and radiuses
    # planet_masses = (4.596, 4.374, 1.297, 2.313, 3.475, 4.418, 1.088) # From Agol+21 Table 2 in e-5*Msun
    # planet_masses = (4.1297145833e-06, 3.9564821698e-06, 1.1720208671e-06, 2.1148831802e-06, 3.1362928581e-06, 3.9669131016e-06, 9.8160839680e-07) # in Msun From AS Libert
    planet_masses = ( 2.1148831802e-06, 3.1362928581e-06 ) # in Msun From AS Libert

    tmp = [ mm *posidonius.constants.M_SUN  for mm in planet_masses]
    print(f'From AS Libert :: \t {tmp}') # Avoid the effect of the masses on the dynamisc --> only take the ≠ on the k2
    print(f'planet radiuses {planet_radiuses}')

    # --- Semi-major axis in AU
    # planet_a = (0.01111, 0.01521, 0.02144, 0.02817, 0.0371, 0.0451, 0.0596)
    # planet_a = (0.01110318e0, 0.01520668e0, 0.02142513e0, 0.02815839e0, 0.03705241e0, 0.04508048e0, 0.05955922e0) #From Bolmont 2020 TTV
    # planet_a = (0.01154e0, 0.01580e0, 0.02227e0, 0.02925e0, 0.03849e0, 0.04683e0, 0.06189e0) #From Agol et al. 2021 Table 6 
    # planet_a = (1.154e-2, 1.580e-2, 2.227e-2, 2.925e-2, 3.849e-2, 4.683e-2, 6.189e-2) #From Agol+21 Table 6 in AU
    # planet_a = ( 0.0115418669425055, 0.0158094975622961, 0.0222700389997286, 0.0292694953252985, 0.0385106150958128, 0.0468445028220842, 0.0619224236258543 ) # in AU From AS Libert
    a1 = 0.0292694953252985
    a2 = a1*1.05*np.power(1.5,2./3.)
    planet_a = ( a1, a2 ) # in AU From AS Libert
    
    # --- Inclination in degrees
    # planet_i = (0.3500000, 0.330000, 0.250000, 0.140000, 0.320000, 0.290000, 0.130000)
    # planet_i = (0.59e0, 0.5e0, 0.30e0, 0.40e0, 0.080e0, 0.37e0, 0.20e0) #From Bolmont 2020 TTV
    planet_i = (0.1, 5.0)
    
    # --- Eccentricity
    # planet_e = (0.00622, 0.00654, 0.0837, 0.0051, 0.01007, 0.00208, 0.00567)
    # planet_e = (3.878e-3, 0.162e-3, 6e-3, 6.24e-3, 8.6e-3, 3.9e-3, 3.73e-3)
    # planet_e = (0.006795, 0.001841, 0.007366, 0.005078, 0.009578, 0.004734, 0.004448) # recomputed from Maximum-Likelihood in Fig~22 of Agol+21
    # planet_e = ( 0.0066945347902135, 0.0019585600861175, 0.0072217598426537, 0.0052151102407526, 0.0095027023237257, 0.0032634092345139, 0.0044037179966920) # From AS Libert 
    planet_e = ( 0.0052151102407526, 0.0095027023237257) # From AS Libert 
    
    # planet_e = (0., 0., 0., 0., 0., 0., 0.) # circular version

    # --- Mean longitude 
    # planet_lamdba = ( 6.2807871243880449, 3.2752103444471459, 0.2445003452092198, 1.4578342582705632, 1.0048074576831532, 1.4551190036825181, 5.0763862710529217)  # in rad From AS Libert 
    planet_lamdba = (1.4578342582705632, 1.0048074576831532)  # in rad From AS Libert 
    
    # --- Argument of pericentre
    # planet_p = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # planet_p = (2.384988, 2.268317, 2.640451, -0.885834, -3.124688, 0.516811, 3.080891) # recomputed from Maximum-Likelihood in Fig~22 of Agol+21 # in rad
    # planet_p = (-89.76412520427041, -89.99217800352784, -89.36936434335934, -89.70701557248958, 62.574227545625696, 89.609586549535, -86.32066658907719)
    
    # --- longitude of ascending node
    # planet_n = (0., 0., 0., 0., 0., 0., 0.)
    planet_long_asc_node = ( 0., 184. )
    
    # --- longitude of pericentre = argument of pericenter + longitude ascending node
    # planet_long_peri = ( 3.9138033639162217, 3.9568810173922708, 3.6462719004641233, 0.9017217708596807, 3.1105825606949371, 5.8092293524825127, 3.2359374491128419) # in rad From AS Libert 
    planet_long_peri = ( 0.9017217708596807, 3.1105825606949371 ) # in rad From AS Libert 
    
    # planet_long_peri = (0., 0., 0., 0., 0., 0., 0.)
    
    # --- Mean anomaly in degrees
    # planet_l = (323.732652895, 96.4925777097, 111.770368348, 165.724187804, 254.117367005, 161.020362506, 134.724813585)
    # planet_l = (90.0e0, 51.5815880e0, 84.5759410e0, 305.4552470e0, 283.5599420e0, 233.7735200e0, 1.3139080e0) #From Bolmont+20 Table~A.2 # in degree
    planet_l =  [ lamdba - g for (lamdba, g) in zip( planet_lamdba, planet_long_peri) ] # compute from Mean longitude and longitude of pericentre
    print(f' -- Recomputed Mean anomaly {planet_l}')

    planets = ['b', 'c', 'd', 'e', 'f', 'g', 'h']
    for planet, r, m, t, a, i, l, e, long_peri, long_asc_node, w_lmpq, ImK2, ReK2, size, gyration_radius_2 in zip(
        planets, planet_radiuses, planet_masses, planet_period, planet_a, planet_i, planet_l, planet_e, planet_long_peri, planet_long_asc_node, planet_w_lm, planet_Imk2, planet_Rek2, planet_nm_data, planet_gyration_radius_squared):
        print(f"\n ===================== Planet {planet} ===================== ")
        print("Nem masses from Agol+21 :", m *posidonius.constants.M_SUN)
        m = m *posidonius.constants.M_SUN
        planet_mass = m /posidonius.constants.M_SUN # Solar masses (3.0e-6 solar masses = 1 earth mass)
        print(f'Mass Trappist \t {m} Kg \t {planet_mass} Msun ' )
        # planet_radius_factor = r # R Earth
        planet_radius = r / posidonius.constants.AU # planet_radius_factor * posidonius.constants.R_EARTH
        planet_radius_of_gyration = np.sqrt(gyration_radius_2) #5.75e-01

        mu = posidonius.constants.G_SI*( star_mass*posidonius.constants.M_SUN + planet_mass*posidonius.constants.M_SUN )
        tmp = np.power(t*86400. /(2.0*np.pi),2)*mu

        sma = np.power(tmp, 1./3.)
        tmp = sma /posidonius.constants.AU
        print("\n--- Check recomputed semi-major axis --- ")
        print(f' From IC {a} | {tmp} | err = {(1-a/tmp)*100.} pc ')
        # a = tmp # sma 

        # --- Specify initial position and velocity for a stable orbit
        # -- Keplerian orbital elements, in the `asteroidal' format of Mercury code
        i = i * posidonius.constants.DEG2RAD;   # inclination (degrees to rad)
        # long_peri = long_peri already in rad * posidonius.constants.DEG2RAD;   # argument of pericentre (degrees)
        long_asc_node = long_asc_node * posidonius.constants.DEG2RAD;  # longitude of the ascending node (degrees)
        # l = l * posidonius.constants.DEG2RAD;   # mean anomaly (degrees)
        long_peri = (long_peri + long_asc_node) # Convert to longitude of perihelion !!
        q = a * (1.0 - e);# perihelion distance
        planet_position, planet_velocity = posidonius.calculate_cartesian_coordinates(planet_mass, q, e, i, long_peri, long_asc_node, l, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
        print("\n--- Check computed cartesian parameters --- ")
        print(f' \t position {planet_position.x()} | {planet_position.y()} | {planet_position.z()} | ')
        print(f' \t angles {np.arctan2(planet_position.y(), planet_position.x())} | ')
        # --- Initialization of planetary spin
        planet_obliquity = 0.0 *posidonius.constants.DEG2RAD; # 1.0e-4 # rad
        # --- Initialize spin
        # period_tmp = 2.*np.pi*(np.sqrt( (1.0*posidonius.constants.AU)**3. / ( posidonius.constants.G_SI*( 1.0*1988500e24 +1.0*5.9722e24 ) ) ) )/86400.
        # print(period_tmp)
        # period_tmp = 2.*np.pi*(np.sqrt( (0.24*posidonius.constants.AU)**3. / ( posidonius.constants.G_SI*( 1.0*1988500e24 + 5.0*5.9722e24 ) ) ) )/86400.
        # print(period_tmp)
        # period_tmp = 2.*np.pi*(np.sqrt( (0.17*posidonius.constants.AU)**3. / ( posidonius.constants.G_SI*( 1.0*1988500e24 + 5.0*5.9722e24 ) ) ) )/86400.
        # print(period_tmp)
        orbital_period = t # 2.*np.pi*(np.sqrt( (a*posidonius.constants.AU)**3. / ( posidonius.constants.G_SI*( star_mass*posidonius.constants.M_SUN + planet_mass*posidonius.constants.M_SUN ) ) ) )/86400.
        planet_rotation_period = orbital_period/(10) # in Days 
        planet_angular_frequency = posidonius.constants.TWO_PI/(planet_rotation_period) # days^-1
        # --- Pseudo-synchronization period
        # planet_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(a, e, star_mass, planet_mass) # days
        # planet_angular_frequency = posidonius.constants.TWO_PI/(planet_pseudo_synchronization_period) # days^-1
        # ---
        # a - sma, q - perihelion distance, e -ecc, i -inclination, p -longitude perihelion, n -longitude asc node, l -mean ano, f -true ano
        planet_keplerian_orbital_elements = posidonius.tools.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
        # planet_inclination = planet_keplerian_orbital_elements[3]
        (planet_sma, planet_peri_distance, planet_ecc, planet_inclination, planet_long_peri, planet_long_asc_node, planet_mean_anomaly, planet_true_ano) = planet_keplerian_orbital_elements
        # planet_spin = posidonius.calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity)
        planet_spin = posidonius.calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity, planet_long_asc_node)

        # print("\n--- Check recomputed orbital parameters --- ")
        # print(f"Mass {m} MSun  {m/5.9722e+24} , Radius {r} Gyration {gyration_radius_2} Nm data {size}  for {32*32}")
        # print(f"sma {planet_sma} | ecc {planet_ecc} | incl {planet_inclination} | long peri {planet_long_peri} | long asc node {planet_long_asc_node} | mean ano {planet_mean_anomaly} | true ano {planet_true_ano} ")
        # print(f"Mean anomaly verif | {planet_mean_anomaly} ?= {l} | {2*np.pi -l} | {2*np.pi +l}" )
        # # print(f'Check Orbital period {orbital_period} with CI {t}')

        print("\n|--- Check "+planet+' : ')
        print("| \t Mass {} Radius {} Rg^2 {} ".format(m,r,gyration_radius_2))
        print("| \t-- Initial orbital elements : ")
        print("| \t sma {}\n| \t ecc {}\n| \t incl {}\n| \t long_asc {}\n| \t long_peri {}\n| \t mean_ano {}\n| \t mean_long {} ".format(a,e,i,long_asc_node%(2.*np.pi),long_peri%(2.*np.pi),l%(2.*np.pi),(l+long_peri)%(2.*np.pi)))
        (sma, _q, ecc, incl, long_peri, long_asc, mean_ano, _true_ano) = planet_keplerian_orbital_elements
        print("| \t-- Computed orbital elements : ")
        print("| \t sma {}\n| \t ecc {}\n| \t incl {}\n| \t long_asc {}\n| \t long_peri {}\n| \t mean_ano {}\n| \t mean_long {} ".format(sma,ecc,incl, long_asc%(2.*np.pi), long_peri%(2.*np.pi), mean_ano%(2.*np.pi), (long_peri+mean_ano)%(2.*np.pi) ))

        k2pdelta = 2.465278e-3 # Terrestrial planets (no gas) 0.0013380687002540815
        planet_tides_parameters = {
            "dissipation_factor_scale": 1.0,
            "dissipation_factor": 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet_radius, 5)),
            "love_number": 0.4 ,#0.299,
        }

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
    
        # --- Choose the tidal model to use :
        # planet_tides_model = posidonius.effects.tides.Kaula(planet_kaula_tidal_parameters_love_numbers)
        planet_tides_model = posidonius.effects.tides.ConstantTimeLag(planet_tides_parameters)

        # --- Choose the type of particle (central body or orbiting body) :
        #planet_tides = posidonius.effects.tides.CentralBody(planet_tides_model)
        planet_tides = posidonius.effects.tides.OrbitingBody(planet_tides_model)
        # planet_tides = posidonius.effects.tides.Disabled()
        # ---
        ReK2_rota_flat = ReK2[0]
        print(f'\t Rek2 for Rotational flattening {ReK2_rota_flat} ')
        planet_rotational_flattening_parameters = {"love_number": ReK2_rota_flat}
        # planet_rotational_flattening_model = posidonius.effects.rotational_flattening.OblateSpheroid(planet_rotational_flattening_parameters)
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
        # disk_surface_density_normalization_gcm = 1000. # g.cm^-2
        # disk_surface_density_normalization_SI = disk_surface_density_normalization_gcm * 1.0e-3 * 1.0e4 # kg.m^-2
        # disk_properties = {
        #     'inner_edge_distance': 0.01,  # AU
        #     'outer_edge_distance': 10.0, # AU
        #     'lifetime': 1.0e5 * 365.25e0, # days
        #     'alpha': 1.0e-2,
        #     'surface_density_normalization': disk_surface_density_normalization_SI * (1.0/posidonius.constants.M_SUN) * posidonius.constants.AU**2, # Msun.AU^-2
        #     'mean_molecular_weight': 2.4,
        # }
        # planet_disk = posidonius.effects.disk.CentralBody(disk_properties)
        planet_disk = posidonius.effects.disk.OrbitingBody()
        # planet_disk = posidonius.effects.disk.Disabled()
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


