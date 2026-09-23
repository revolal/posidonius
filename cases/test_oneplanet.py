import posidonius
import numpy as np
import argparse
import os, sys
from posidonius.particles.axes import Axes

if __name__ == "__main__":
    def calculate_eccentricity_vector(gm, position, velocity):

        # // --- Input --- //
        x = position.x()
        y = position.y()
        z = position.z()
        u = velocity.x()
        v = velocity.y()
        w = velocity.z()

        # // --- Output --- //
        # eccentricity_vector

        # // --- Local --- //
        # // Angular momentum
        hx = y * w  -  z * v
        hy = z * u  -  x * w
        hz = x * v  -  y * u
        # // v vectorial h
        v_vect_h_x = v * hz - w * hy
        v_vect_h_y = w * hx - u * hz
        v_vect_h_z = u * hy - v * hx
        # // distance
        r = np.sqrt(np.power(x,2) + np.power(y,2) + np.power(z,2) )
        # // eccentricity componant
        eccentricity_vector = []
        eccentricity_vector.append( (v_vect_h_x / gm) - (x / r))
        eccentricity_vector.append( (v_vect_h_y / gm) - (y / r))
        eccentricity_vector.append( (v_vect_h_z / gm) - (z / r))

        return eccentricity_vector

    parser = argparse.ArgumentParser()
    parser.add_argument('output_filename', action='store', help='Filename where the initial snapshot will be stored (e.g., universe_integrator.json)')

    args = parser.parse_args()
    filename = args.output_filename
    #filename = posidonius.constants.BASE_DIR+"target/case7.json"

    initial_time = 1.0*365.25 # time [days] where simulation starts
    time_step = 0.3 # 10.0 # days
    time_limit = 1.5e8 *365.25 # 50.0 * 365.25 # 18.772866 *100.# 2.*365.25 #* 1.0e2 # days
    historic_snapshot_period = 1.0e3 #100.*time_step # 1.*365.25 # days
    recovery_snapshot_period =  1.0e6 #10.*historic_snapshot_period # days
    consider_effects = posidonius.ConsiderEffects({
        "tides": True,
        "rotational_flattening": False,
        "general_relativity": False,
        "disk": False,
        "wind": False,
        "evolution": False,
    })
    universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects)

    star_mass = 1.0 # 0.0898 # Solar masses Agol+21
    star_radius_factor = 0.119 #0.117 # Agol+21
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
    # star_tides = posidonius.effects.tides.OrbitingBody(star_tides_model)
    # star_tides = posidonius.effects.tides.Disabled()
    #
    star_rotational_flattening_parameters = {"love_number": star_tides_parameters["love_number"] }
    # star_rotational_flattening_model = posidonius.effects.rotational_flattening.OblateSpheroid(star_rotational_flattening_parameters)
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
    disk_surface_density_normalization_gcm = 1000. # g.cm^-2
    disk_surface_density_normalization_SI = disk_surface_density_normalization_gcm * 1.0e-3 * 1.0e4 # kg.m^-2
    disk_properties = {
        'inner_edge_distance': 0.02,  # AU
        'outer_edge_distance': 0.04, # AU
        'lifetime': 100 * 365.25e0, # days
        'alpha': 1.0e-2,
        'surface_density_normalization': disk_surface_density_normalization_SI * (1.0/posidonius.constants.M_SUN) * posidonius.constants.AU**2, # Msun.AU^-2
        'mean_molecular_weight': 2.4,
    }
    # star_disk = posidonius.effects.disk.CentralBody(disk_properties)
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

    # === PLANET PARAMETERS ===
    # dir_path = './input/love_numbers/corrected_profiles/andrade/earth/'
    # dir_path = './input/love_numbers/corrected_profiles/sunbercooper/earth/'
    # dir_path = '../Revol_letter/imk2_t1_earth/'
    # dir_path = '../Revol_letter/imk2_t1_earth_visccorr/'
    dir_path = './input/love_numbers/TRAPPIST-1_Earth-like/' # Love numbers files Earth-like core composition
    # dir_path = "/home/arevol/Desktop/Dynaclim/posidonius/input/love_numbers/data_love_number_Homogen_for_posidonius/" # Love numbers files Earth-like core composition

    file_dir_path = []
    for path in os.listdir(dir_path):
        print(path)
        if os.path.isfile(os.path.join(dir_path, path)):
            if path.endswith('.txt'):
                file_dir_path.append(path)
    file_dir_path.sort()
    print(file_dir_path)

    planet_param = []
    # planet_masses = []
    # planet_radiuses = []
    # planet_gyration_radius_squared = []
    planet_w_lm = []
    planet_Imk2 = []
    planet_Rek2 = []
    planet_nm_data = []

    for path in file_dir_path:
        print(dir_path+path)
        planet_data = np.loadtxt(dir_path+path,comments='#',skiprows=1,)#delimiter=',')
        # mass, radius, gyration_radius_squared = planet_data[0,:]
        print("infile ",planet_data[0,:])
        # planet_masses.append(mass)
        # planet_radiuses.append(radius)
        # planet_gyration_radius_squared.append(gyration_radius_squared)
        planet_param.append(planet_data[0,:])
        planet_w_lm.append(planet_data[1:,0])
        planet_Imk2.append(planet_data[1:,1])
        planet_Rek2.append(planet_data[1:,2])
        planet_nm_data.append(np.size(planet_data[1:,0]))
    planet_param = planet_param[0]
    planet_w_lm = planet_w_lm[0]
    planet_Imk2 = planet_Imk2[0]
    planet_Rek2 = planet_Rek2[0]
    planet_nm_data = planet_nm_data[0]
    


    ############################################################################
    Mearth = 5.9736E+24 #Earth Mass (kg)
    Rearth = 6378.1E+3    #Earth radius (m)
    Mvenus = 4.8685E+24
    Rvenus = 6052E+3
    Mjup = 1.898E+27 # Jupiter mass (kg)
    Rjup = 6.991E+07 # Jupiter radius (m)

    # --- Planets periods in  days
    planet_period = 1.5112868036098193

    # --- Planets masses and radiuses
    planet_masses = Mvenus /posidonius.constants.M_SUN # in M_SUN 
    planet_radiuses = Rvenus
    planet_gyration_radius = 0.5559429772345046
    planet_gyration_radius_squared = planet_gyration_radius *planet_gyration_radius 

    # --- Semi-major axis in AU
    planet_a = 0.0723

    # --- Inclination in degrees
    planet_i = 0.

    # --- Eccentricity
    planet_e = 0.1
    # planet_e =0.

    # --- mean_longitude
    planet_mean_long = 0.0 # 6.2807871243880449

    # --- Argument of pericentre
    # planet_p = 0.0

    # --- longitude of ascending node
    planet_n = 0.

    # --- longitude of pericentre = argument of pericenter + longitude ascending node
    planet_p = 0.0 # 3.9138033639162217
    # planet_p = (0., 0., 0., 0., 0., 0., 0.)

    # --- Mean anomaly in degrees
    # planet_l = (323.732652895, 96.4925777097, 111.770368348, 165.724187804, 254.117367005, 161.020362506, 134.724813585)
    # planet_l = (90.0e0, 51.5815880e0, 84.5759410e0, 305.4552470e0, 283.5599420e0, 233.7735200e0, 1.3139080e0) #From Bolmont+20 Table~A.2 # in degree
    planet_l = (planet_mean_long - planet_p )%(2.*np.pi)  # compute from Mean longitude and longitude of pericentre
    print(f' -- Recomputed Mean anomaly {planet_l}')

    planets = ['Venus']
    (planet, r, m, t, a, i, l, e, p, n, w_lmpq, ImK2, ReK2, size, gyration_radius_2) = (
    planets, planet_radiuses, planet_masses, planet_period, planet_a, planet_i, planet_l, planet_e, planet_p, planet_n, planet_w_lm, planet_Imk2, planet_Rek2, planet_nm_data, planet_gyration_radius_squared)
    print(f"\n ===================== Planet {planet} ===================== ")
    m = m *posidonius.constants.M_SUN
    planet_mass = m /posidonius.constants.M_SUN # Solar masses (3.0e-6 solar masses = 1 earth mass)
    print(f'Mass \t {m} Kg \t {planet_mass} Msun ' )
    # planet_radius_factor = r # R Earth
    planet_radius = r / posidonius.constants.AU # planet_radius_factor * posidonius.constants.R_EARTH
    planet_radius_of_gyration = np.sqrt(gyration_radius_2) #5.75e-01

    mu = posidonius.constants.G_SI*( star_mass*posidonius.constants.M_SUN + planet_mass*posidonius.constants.M_SUN )
    t = np.sqrt(np.power(a * posidonius.constants.AU,3) /mu  )*2.0*np.pi / 86400.
    # tmp = np.power(t*86400. /(2.0*np.pi),2)*mu
    # sma = np.power(tmp, 1./3.)
    # sma = sma /posidonius.constants.AU

    print("\n--- Check orbital elements --- ")
    print(f'\t Sma {a} | ecc {e} | incl {i} | long_peri {p} | long_asc_node {n} | mean_ano {l} | mean_long {(l+p)%(2.*np.pi)}  ')
    print(f'\t Orbital period {t} | Orbital frequency {2.0*np.pi/(t*24.*3600.)}')
    # a = sma # sma
    # print(" - Check recomputed sma --- ")
    # print(f' From IC {a} | Cacl {sma} | err = {(1-a/tmp)*100.} pc ')

    # --- Specify initial position and velocity for a stable orbit
    # -- Keplerian orbital elements, in the `asteroidal' format of Mercury code
    i = 0. * posidonius.constants.DEG2RAD;   # inclination (degrees)
    # p = p already in rad * posidonius.constants.DEG2RAD;   # argument of pericentre (degrees)
    # n = n * posidonius.constants.DEG2RAD;  # longitude of the ascending node (degrees)
    # l = l * posidonius.constants.DEG2RAD;   # mean anomaly (degrees)
    p = (p + n) # Convert to longitude of perihelion !!
    q = a * (1.0 - e);# perihelion distance
    planet_position, planet_velocity = posidonius.calculate_cartesian_coordinates(planet_mass, q, e, i, p, n, l, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
    print("\n--- Check computed cartesian parameters --- ")
    print(f' \t position {planet_position.x()} | {planet_position.y()} | {planet_position.z()} | ')
    angle_position_xy = np.arctan2(planet_position.y(), planet_position.x())
    print(f' \t angles position planet {angle_position_xy} | {(angle_position_xy)%(2.*np.pi)} | ')
    # --- Initialization of planetary spin
    planet_obliquity = 0.0 # 1.0e-4 # rad
    # --- Initialize spin
    orbital_period = t # 2.*np.pi*(np.sqrt( (a*posidonius.constants.AU)**3. / ( posidonius.constants.G_SI*( star_mass*posidonius.constants.M_SUN + planet_mass*posidonius.constants.M_SUN ) ) ) )/86400.
    planet_rotation_period = orbital_period/(2.6) # in Days
    planet_angular_frequency = posidonius.constants.TWO_PI/(planet_rotation_period) # days^-1
    # --- Pseudo-synchronization period
    # planet_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(a, e, star_mass, planet_mass) # days
    # planet_angular_frequency = posidonius.constants.TWO_PI/(planet_pseudo_synchronization_period) # days^-1
    # ---
    # a - sma, q - perihelion distance, e -ecc, i -inclination, p -longitude perihelion, n -longitude asc node, l -mean ano, f -true ano
    planet_keplerian_orbital_elements = posidonius.tools.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
    # planet_inclination = planet_keplerian_orbital_elements[3]
    (planet_sma, planet_peri_distance, planet_ecc, planet_inclination, planet_long_peri, planet_long_asc_node, planet_mean_anomaly, planet_true_ano) = planet_keplerian_orbital_elements
    planet_inclination = planet_keplerian_orbital_elements[3]
    planet_long_asc_node = planet_keplerian_orbital_elements[5]
    planet_spin = posidonius.calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity, planet_long_asc_node)
    planet_true_longitude  = (planet_true_ano + planet_long_peri)%(2.*np.pi)

    # ---
    print('> Comtpute True anomaly ')
    gm = posidonius.constants.G *(star_mass+planet_mass)
    vec_position = planet_position
    vec_velocity = planet_velocity
    (ex, ey, ez) = calculate_eccentricity_vector(gm, vec_position, vec_velocity)

    check_long_per = np.arctan2(ey,ex)
    print(f" \t Check long peri  {check_long_per} | {check_long_per%(2.*np.pi)}")
    # ---

    print("\n--- Check recomputed orbital parameters --- ")
    print(f"Mass {m} MSun  {m/5.9722e+24} MEarth , Radius {r} Gyration {gyration_radius_2} Nm data {size}  for {32*32}")
    # print(f'\t  Sma {planet_sma} | ecc {e} | incl {i} | long_peri {p} | mean_ano {l} | mean_long {l+p} ')
    print(f"\t Sma {planet_sma} | ecc {planet_ecc} | incl {planet_inclination} | long peri {planet_long_peri} | mean_ano {planet_mean_anomaly}  | long_asc_node {planet_long_asc_node} |  mean_long {planet_long_peri+planet_mean_anomaly} | true_ano {planet_true_ano} | true_long {planet_true_longitude} ")
    # print(f'Check Orbital period {orbital_period} with CI {t}')

    print("\n--- Compare recomputed orbital parameters --- ")
    print(f"\t Sma {a} | {planet_sma} ")
    print(f"\t ecc {e} | {planet_ecc} ")
    print(f"\t incl {i} | {planet_inclination} ")
    print(f"\t long peri {p} | {planet_long_peri%(2.*np.pi)}  ")
    print(f"\t mean_ano {l} | {planet_mean_anomaly}  ")
    print(f"\t long_asc_node {n} | {planet_long_asc_node} ")
    print(f"\t mean_long {(l+p)%(2.*np.pi)} |  {(planet_long_peri+planet_mean_anomaly)%(2.*np.pi)} ")
    print("\n--- --- ")

    # k2pdelta = 2.465278e-3 # Terrestrial planets (no gas) 0.0013380687002540815
    # planet_tides_parameters = {
    #     "dissipation_factor_scale": 1.0,
    #     "dissipation_factor": 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet_radius, 5)),
    #     "love_number": 0.4 ,#0.299,
    # }

    expected_size = 1024
    if size > expected_size:
        raise Exception("This should not happen")
    planet_kaula_tidal_parameters_love_numbers = {
        "love_numbers": {
            "spectrum_excitation_frequency": w_lmpq.tolist() + [0] * (expected_size - len(w_lmpq)),
            "spectrum_real_part": ReK2.tolist() + [0] * (expected_size - len(ReK2)),
            "spectrum_imaginary_part": ImK2.tolist() + [0] * (expected_size - len(ImK2)),
        },
    }

    # --- Choose the tidal model to use :
    planet_tides_model = posidonius.effects.tides.Kaula(planet_kaula_tidal_parameters_love_numbers)
    # planet_tides_model = posidonius.effects.tides.ConstantTimeLag(planet_tides_parameters)

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
    # planet_disk = posidonius.effects.disk.OrbitingBody()
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


