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
    time_step = 0.01 # days
    time_limit = 1500. # 1.* 365.25 # 18.772866 *100.# 2.*365.25 #* 1.0e2 # days
    historic_snapshot_period = 1.0 # 0.1*365.25 #100.*time_step # 1.*365.25 # days
    recovery_snapshot_period =  10. * 365.25  #10.*historic_snapshot_period # days
    consider_effects = posidonius.ConsiderEffects({
        "tides": False,
        "rotational_flattening": False,
        "general_relativity": False,
        "disk": False,
        "wind": False,
        "evolution": False,
    })
    universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects)

    center_of_mass_position = (1.4173038216929595e-06, 0.0, 3.223406102928513e-06)
    center_of_mass_velocity = (2.183475087802898e-06, 0.0, -1.610681439768188e-06)

    star_mass = 0.0898 # Solar masses Agol+21
    star_radius_factor = 0.119 #0.117 # Agol+21
    star_radius = star_radius_factor * posidonius.constants.R_SUN
    star_radius_of_gyration = 4.47e-01 # Brown dwarf
    # star_position = posidonius.Axes(0., 0., 0.)
    # star_velocity = posidonius.Axes(0., 0., 0.)
    # star_position = posidonius.Axes(-1.4173038216892925e-06-center_of_mass_position[0], 0.0-center_of_mass_position[1], -3.22340610293452e-06-center_of_mass_position[2])
    # star_velocity = posidonius.Axes(-2.18347508780533e-06-center_of_mass_velocity[0], 0.0-center_of_mass_velocity[1], 1.6106814397688342e-06-center_of_mass_velocity[2])
    # star_position = posidonius.Axes(-0.006228200996557742-center_of_mass_position[0], 0.0-center_of_mass_position[1], -0.014164938240686467-center_of_mass_position[2])
    # star_velocity = posidonius.Axes(-0.009595064593538568-center_of_mass_velocity[0], 0.0-center_of_mass_velocity[1], 0.0070779797491257405-center_of_mass_velocity[2])
    star_position = posidonius.Axes(-1.4173038216892925e-06, 0.0, -3.22340610293452e-06)
    star_velocity = posidonius.Axes(-2.18347508780533e-06, 0.0, 1.6106814397688342e-06)

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
    file_dir_path = []
    for path in os.listdir(dir_path):
        print(path)
        if os.path.isfile(os.path.join(dir_path, path)):
            if path.endswith('.txt'):
                file_dir_path.append(path)
    file_dir_path.sort()
    print(file_dir_path)

    planet_param = []
    planet_masses = []
    planet_radiuses = []
    planet_gyration_radius_squared = []
    planet_w_lm = []
    planet_Imk2 = []
    planet_Rek2 = []
    planet_nm_data = []

    for path in file_dir_path:
        print(dir_path+path)
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


    ############################################################################
    print(f'From Burnman :: \t {planet_masses} ')

    # === From Teyssandier IC in Msun
    planet_masses = (
        4.129714583336683e-06,
        3.9564821698106564e-06, 
        1.172020867077638e-06, 
        2.11488318016312e-06, 
        3.136292858082315e-06, 
        3.966913101584511e-06, 
        9.81608396760074e-07) # in Msun
    planet_radiuses = (
        1.9383135539813295e-09, 
        1.905313592040787e-09,
        1.3686300004814403e-09, 
        1.597892893963103e-09,
        1.8149979067298288e-09, 
        1.960892475309069e-09,
        1.311314277111025e-09 
        ) # in AU
    planet_gyration_radius = [0.5559429772345046, 0.5550938677788629, 0.5730245432973096, 0.5646295658126751, 0.566372994597031, 0.5688976599660345, 0.5823930468590962]
    # planet_gyration_radius_squared = [rg*rg for rg in planet_gyration_radius]

    # # Inertial position-velocities in x-z plane from Rebound-Teyssandier IC
    planet_positions = (
        (0.011595646564314931,  0.0, 7.692286051067042e-05  ),
        (-0.01565080907036833,  0.0, -0.0020674747364663017 ), 
        (0.021736842692296483,  0.0, 0.005505383144540553   ),
        (0.00312307880454434,   0.0, 0.028968857594787387   ),
        (0.021279615835211273,  0.0, 0.03232075140847502    ),
        (0.0051280978729623845, 0.0, 0.04661585882380491    ),
        (0.022561454967478176,  0.0, -0.057744464283335996  ),
    ) # in AU
    # planet_positions = [ [ p-c for p,c in zip(position, center_of_mass_position)] for position in planet_positions]
    planet_velocities = (
        (0.00010529805824325997, 0.0, -0.04775357838716252),
        (-0.005422352581315724,  0.0, 0.04070376647409209),
        (0.008363205071675964,   0.0, -0.03326888883914474),
        (0.030081673878636698,   0.0, -0.0033276842928732803),
        (0.02194981180982474,    0.0, -0.014196911633469764),
        (0.02364154251490609,    0.0, -0.002673872398002983),
        (-0.01930648365922732,   0.0, -0.007449083614085322),
    ) # in AU/days
    # planet_velocities = [ [ v-c for v,c in zip(velocity, center_of_mass_velocity)] for velocity in planet_velocities]
    # ## Test new inertial coordinate (corrected from center of mass in Rebound)
    # planet_positions = (
    #     (0.005368862871578878,       0.0, -0.014084791974072863),
    #     (-0.021877592763104383,      0.0, -0.016229189571049834),
    #     (0.015510058999560431,       0.0, -0.00865633169004298),
    #     (-0.003103704888191713,      0.0, 0.014807142760203854),
    #     (0.015052832142475221,       0.0, 0.018159036573891486),
    #     (-0.0010986858197736685,     0.0, 0.032454143989221376),
    #     (0.016334671274742123,       0.0, -0.07190617911791952)
    # ) # in AU
    # planet_positions = [ [ p-c for p,c in zip(position, center_of_mass_position)] for position in planet_positions]
    # planet_velocities = (
    #     (-0.009487583060207503 ,     0.0, -0.04067720931947655),
    #     (-0.015015233699766487 ,     0.0, 0.04778013554177806),
    #     (-0.0012296760467747983,     0.0, -0.026192519771458773),
    #     (0.020488792760185934  ,     0.0, 0.003748684774812691),
    #     (0.01235693069137398   ,     0.0, -0.007120542565783793),
    #     (0.014048661396455329  ,     0.0, 0.004402496669682988),
    #     (-0.028899364777678085 ,     0.0, -0.0003727145463993508)
    # ) # in AU/days
    # planet_velocities = [ [ v-c for v,c in zip(velocity, center_of_mass_velocity)] for velocity in planet_velocities]


    # # Heliocentric position-velocities in x-z plane from Rebound-Teyssandier IC
    # planet_positions = (
    #     (0.011597063868136621,  0.0, 8.014626661360493e-05),
    #     (-0.01564939176654664,  0.0, -0.0020642513303633673),
    #     (0.021738259996118173,  0.0, 0.0055086065506434875),
    #     (0.0031244961083660293, 0.0, 0.028972081000890323),
    #     (0.021281033139032964,  0.0, 0.03232397481457795),
    #     (0.005129515176784074,  0.0, 0.04661908222990784),
    #     (0.022562872271299866,  0.0, -0.05774124087723306),
    # ) # in AU

    # planet_velocities = (
    #     (0.0001074815333310653, 0.0, -0.047755189068602294),
    #     (-0.005420169106227919, 0.0, 0.040702155792652316),
    #     (0.00836538854676377,   0.0, -0.033270499520584515),
    #     (0.030083857353724503,  0.0, -0.003329294974313049),
    #     (0.021951995284912547,  0.0, -0.014198522314909532),
    #     (0.023643725989993897,  0.0, -0.0026754830794427516),
    #     (-0.019304300184139516, 0.0, -0.007450694295525091),
    # ) # in AU/days

    # Heliocentric position-velocities in x-y plane from Rebound-Teyssandier IC
    # planet_positions = (
    #     (0.011597063868136621,  8.014626661360493e-05 , 0.0),
    #     (-0.01564939176654664,  -0.0020642513303633673, 0.0),
    #     (0.021738259996118173,  0.0055086065506434875 , 0.0),
    #     (0.0031244961083660293, 0.028972081000890323  , 0.0),
    #     (0.021281033139032964,  0.03232397481457795   , 0.0),
    #     (0.005129515176784074,  0.04661908222990784   , 0.0),
    #     (0.022562872271299866,  -0.05774124087723306  , 0.0),
    # ) # in AU
    # planet_velocities = (
    #     (0.0001074815333310653, -0.047755189068602294 , 0.0),
    #     (-0.005420169106227919, 0.040702155792652316  , 0.0),
    #     (0.00836538854676377,   -0.033270499520584515 , 0.0),
    #     (0.030083857353724503,  -0.003329294974313049 , 0.0),
    #     (0.021951995284912547,  -0.014198522314909532 , 0.0),
    #     (0.023643725989993897,  -0.0026754830794427516, 0.0),
    #     (-0.019304300184139516, -0.007450694295525091 , 0.0),
    # ) # in AU/days

    planet_period =(1.5112868036098193, 2.4227038555269917, 4.0504368087076257, 6.1029264019452052, 9.2103976604134150, 12.3562321332659266, 18.7787920964476669) # in days From AS Libert


    planets = ['b', 'c', 'd', 'e', 'f', 'g', 'h']
    for planet, r, m, t, xyz, vxvyvz, w_lmpq, ImK2, ReK2, size, gyration_radius in zip(
        planets, planet_radiuses, planet_masses, planet_period, planet_positions, planet_velocities, planet_w_lm, planet_Imk2, planet_Rek2, planet_nm_data, planet_gyration_radius):
        print(f"\n ===================== Planet {planet} ===================== ")
        print("Nem masses from Teyssandier + :", m *posidonius.constants.M_SUN)
        planet_mass = m # Solar masses (3.0e-6 solar masses = 1 earth mass)
        print(f'Mass Trappist \t {m *posidonius.constants.M_SUN} Kg \t {planet_mass} Msun ' )
        planet_radius = r
        planet_radius_of_gyration = gyration_radius #5.75e-01

        print("\n--- Check c cartesian parameters --- ")
        print(f' \t position {xyz[0]} | {xyz[1]} | {xyz[2]} | ')
        # --- Initialization of planetary spin
        planet_obliquity = 0.0 # rad
        # --- Initialize spin
        orbital_period = t 
        planet_rotation_period = orbital_period/(1.0) # in Days
        planet_angular_frequency = posidonius.constants.TWO_PI/(planet_rotation_period) # days^-1µ
        planet_spin = posidonius.calculate_spin(planet_angular_frequency, 0.0, 0.0, 0.0)

        planet_position = posidonius.Axes(xyz[0], xyz[1], xyz[2])
        planet_velocity = posidonius.Axes(vxvyvz[0], vxvyvz[1], vxvyvz[2])

        # === Check cartesian input 
        planet_keplerian_orbital_elements = posidonius.tools.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
        a, q, eccentricity, i, p, n, l, f = planet_keplerian_orbital_elements
        print(f"Check : sma : {a} | ecc : {eccentricity} | inclination : {i} ")


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

    # whfast_alternative_coordinates="DemocraticHeliocentric"
    #whfast_alternative_coordinates="WHDS"
    whfast_alternative_coordinates="Jacobi"
    universe.write(filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)
    # universe.write(filename, integrator="IAS15")
    #universe.write(filename, integrator="LeapFrog")


