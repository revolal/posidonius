import posidonius
import numpy as np
import argparse
import os
from posidonius.particles.axes import Axes
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('output_filename', action='store', help='Filename where the initial snapshot will be stored (e.g., universe_integrator.json)')

    args = parser.parse_args()
    filename = args.output_filename

    initial_time = 1.0 * 365.25 # time [days] where simulation starts
    time_step = 0.5 # days
    time_limit = 1.0e6 * 356.25 # DAY
    historic_snapshot_period = 1.0e2 * 365.25 # DAY
    recovery_snapshot_period = 1.0e3 * 365.25  # DAY
    consider_effects = posidonius.ConsiderEffects({
        "tides": False,
        "rotational_flattening": False,
        "general_relativity": False,
        "disk": True,
        "wind": False,
        "evolution": False,
    })
    universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects)

    star_mass = 0.77 # M_SUN 
    star_radius_factor = 0.77 # R_SUN
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
    # star_tides = posidonius.effects.tides.CentralBody(star_tides_model)
    # star_tides = posidonius.effects.tides.OrbitingBody(star_tides_model)
    star_tides = posidonius.effects.tides.Disabled()
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
    disk_surface_density_normalization_gcm = 800. # g.cm⁻²
    disk_surface_density_normalization_SI = disk_surface_density_normalization_gcm * 1.0e-3 * 1.0e4 # kg.m⁻²
    disk_properties = {
        'inner_edge_distance': 0.15,    # AU
        'outer_edge_distance': 2,       # AU
        'lifetime': 2.5e4 * 365.25e0,     # DAY
        'alpha': 2.5e-3,
        'surface_density_normalization': disk_surface_density_normalization_SI * (1.0/posidonius.constants.M_SUN) * posidonius.constants.AU**2, # M_SUN.AU⁻²
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


    # === PLANET TIDAL PARAMETERS ======================
    dir_path = './input/love_numbers/Rheo_homogeneous/'
    file_dir_path = []
    for path in os.listdir(dir_path):
        print(path)
        if os.path.isfile(os.path.join(dir_path, path)):
            if path.endswith('.csv'):
                file_dir_path.append(path)
    file_dir_path.sort()
    print(file_dir_path)

    planets_w_lm = []
    planets_Imk2 = []
    planets_Rek2 = []
    planets_nm_data = []
    for path in file_dir_path:
        print(dir_path+path)
        df = pd.read_csv(dir_path+path, skiprows=1)
        w = df.iloc[:, 0].values 
        i = df.iloc[:, 1].values 
        r = df.iloc[:, 2].values 
        planets_w_lm.append(w)
        planets_Imk2.append(i)
        planets_Rek2.append(r)
        planets_nm_data.append(np.size(w))
    

    # === PLANET ORBITAL PARAMETERS
    # print(f'From data files :: \t {planet_masses} ')

    # --- Planets periods
    planets_period =(
        17.16073,
        34.525528,
        ) # in DAY

    planets_sma =(
        0.125, #0.175, #0.1192,
        0.20, # 0.262, #0.1900
        ) # in AU

    planets_masses = (
        19,
        178,
        ) # in M_EARTH
    planets_masses = [p / posidonius.constants.M2EARTH for p in planets_masses ] # in M_SUN
    planets_radiuses = (
        8.0,
        10.1,
        ) # in R_EARTH
    planets_radiuses = [ r / posidonius.constants.AU for r in planets_radiuses ] # in AU

    planets_gyration_radius = (
        0.5559429772345046, 0.5550938677788629,
        ) # [TO CHANGE]
    planets_gyration_radius_squared = [rg*rg for rg in planets_gyration_radius]


    # --- Inclination in degrees
    planets_i = (
        0.0,
        0.0,
    )

    # --- Eccentricity
    planets_e = (
        0.0, #0.16,
        0.0, #0.0946,
    )

    # --- mean_longitude
    # planet_mean_long = ( 6.2807871243880449, 3.2752103444471459, 0.2445003452092198, 1.4578342582705632, 1.0048074576831532, 1.4551190036825181, 5.0763862710529217)  # in rad From AS Libert
    planets_mean_long = (
        0.0,
        0.0,
    )

    # --- longitude of ascending node
    planets_n = ( 0., 0., 0., 0., 0., 0., 0.)

    # --- longitude of pericentre = argument of pericenter + longitude ascending node
    planets_p = (
        0.0, #291.8,
        0.0, #190.0,
    )
    planets_p = [ p * np.pi / 180. for p in planets_p ]

    # --- Mean anomaly 
    planets_l =  [ (mean_long - long_peri)%(2.*np.pi) for (mean_long, long_peri) in zip(planets_mean_long, planets_p) ] # compute from Mean longitude and longitude of pericentre
    print(f' -- Recomputed Mean anomaly {planets_l}')


    planets = ['b', 'c',]
    for planet, planet_radius, planet_mass, planet_period, planet_sma, planet_i, planet_l, planet_e, planet_p, planet_n, w_lmpq, ImK2, ReK2, size, gyration_radius_2 in zip(
        planets, planets_radiuses, planets_masses, planets_period, planets_sma, planets_i, planets_l, planets_e, planets_p, planets_n, planets_w_lm, planets_Imk2, planets_Rek2, planets_nm_data, planets_gyration_radius_squared):
        print(f"\n ===================== Planet {planet} ===================== ")
        print(f'Mass \t {planet_mass *posidonius.constants.M_SUN } kg \t {planet_mass / posidonius.constants.M_EARTH } Mearth \t {planet_mass} Msun ' )
        planet_radius_of_gyration = np.sqrt(gyration_radius_2) #5.75e-01

        mu = posidonius.constants.G_SI*( star_mass*posidonius.constants.M_SUN + planet_mass*posidonius.constants.M_SUN )
        tmp = np.power(planet_period *86400. /(2.0*np.pi),2)*mu

        sma = np.power(tmp, 1./3.) / posidonius.constants.AU
        sma = planet_sma # sma /posidonius.constants.AU
        print("\n--- Check orbital elements --- ")
        print(f'\t Sma {sma} | ecc {planet_e} | incl {planet_i} | long_peri {planet_p} | long_asc_node {planet_n} | mean_ano {planet_l} | mean_long {(planet_l+planet_p)%(2.*np.pi)}  ')
        print(f'\t Orbital period {planet_period} | Orbital frequency {2.0*np.pi/(planet_period*24.*3600.)}')

        # --- Specify initial position and velocity for a stable orbit
        # -- Keplerian orbital elements, in the `asteroidal' format of Mercury code
        # i = 0. * posidonius.constants.DEG2RAD;   # inclination (degrees)
        # p = p already in rad * posidonius.constants.DEG2RAD;   # argument of pericentre (degrees)
        # n = n * posidonius.constants.DEG2RAD;  # longitude of the ascending node (degrees)
        # l = l * posidonius.constants.DEG2RAD;   # mean anomaly (degrees)
        # planet_p = (planet_p + planet_n) # Convert to longitude of perihelion !!
        planet_q = sma * (1.0 - planet_e);# perihelion distance
        planet_position, planet_velocity = posidonius.calculate_cartesian_coordinates(planet_mass, planet_q, planet_e, planet_i, planet_p, planet_n, planet_l, masses=[star_mass], positions=[star_position], velocities=[star_velocity])

        # --- Initialization of planetary spin
        planet_obliquity = 0.0 # 1.0e-4 # rad
        # --- Initialize spin
        # orbital_period = 2.*np.pi*(np.sqrt( (sma*posidonius.constants.AU)**3. / ( posidonius.constants.G_SI*( star_mass*posidonius.constants.M_SUN + planet_mass*posidonius.constants.M_SUN ) ) ) )/86400.
        planet_rotation_period = planet_period # orbital_period/(1.0) # in Days
        # print("> Planet rotation period :", planet_rotation_period, orbital_period/planet_rotation_period )
        planet_angular_frequency = posidonius.constants.TWO_PI/(planet_rotation_period) # days^-1

        # ---
        # a - sma, q - perihelion distance, e -ecc, i -inclination, p -longitude perihelion, n -longitude asc node, l -mean ano, f -true ano
        planet_keplerian_orbital_elements = posidonius.tools.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
        # planet_inclination = planet_keplerian_orbital_elements[3]
        (planet_sma, planet_peri_distance, planet_ecc, planet_inclination, planet_long_peri, planet_long_asc_node, planet_mean_anomaly, planet_true_ano) = planet_keplerian_orbital_elements
        planet_inclination = planet_keplerian_orbital_elements[3]
        planet_long_asc_node = planet_keplerian_orbital_elements[5]
        planet_spin = posidonius.calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity, planet_long_asc_node)
        planet_true_longitude  = (planet_true_ano + planet_long_peri)%(2.*np.pi)


        print("\n--- Check recomputed orbital parameters --- ")
        # print(f'\t  Sma {planet_sma} | ecc {e} | incl {i} | long_peri {p} | mean_ano {l} | mean_long {l+p} ')
        print(f"\t Sma {planet_sma} | ecc {planet_ecc} | incl {planet_inclination} | long peri {planet_long_peri} | mean_ano {planet_mean_anomaly}  | long_asc_node {planet_long_asc_node} |  mean_long {planet_long_peri+planet_mean_anomaly} | true_ano {planet_true_ano} | true_long {planet_true_longitude} ")
        # print(f'Check Orbital period {orbital_period} with CI {t}')

        # Set Tidal parameters
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
        # planet_tides = posidonius.effects.tides.OrbitingBody(planet_tides_model)
        planet_tides = posidonius.effects.tides.Disabled()
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


