import posidonius
from posidonius import constants
import posidonius.tools
import numpy as np
import argparse
import math

if __name__ == "__main__":
    print("\n -------------------------------------- ")
    print(" --- Trappist 1e .py initialization --- \n")
    parser = argparse.ArgumentParser()
    parser.add_argument('output_filename', action='store', help='Filename where the initial snapshot will be stored (e.g., universe_integrator.json)')

    args = parser.parse_args()
    filename = args.output_filename
    #filename = posidonius.constants.BASE_DIR+"target/case7.json"

    #initial_time = 4.5e6*365.25 # time [days] where simulation starts
    initial_time = 1.0e6*365.25 # time [days] where simulation starts
    time_step = 0.01 # days
    time_step = 0.000008 # days
    #time_limit   = 4*time_step # days
    #time_limit   = 365.25 * 1.0e8 # days
    time_limit   = 365.25 * 1.0e10 # days
    historic_snapshot_period = 10.*365.25 # days
    recovery_snapshot_period = historic_snapshot_period#10.*historic_snapshot_period # days
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
    # star_tides = posidonius.effects.tides.OrbitingBody(star_tides_parameters)
    # star_tides = posidonius.effects.tides.Disabled()
    #
    # star_rotational_flattening_parameters = {"love_number": star_tides_parameters["love_number"] }
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
    print(" --- Recover K2 data --- ")
    

    #===============================================================================================
    # planet_data_b = np.loadtxt('Results_Trappist_1b_20.5_100_freq_Imk2_posidonius.txt',comments='#')
    planet_data_b = np.loadtxt('Results_Creep_Trappist_1b_20.5_100_freq_Imk2_posidonius.txt',comments='#')
    # planet_data_b = np.loadtxt('Results_Creep-2_Trappist_1b_20.5_100_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_b, planet_radius_b, planet_gyration_radius_squared_b = planet_data_b[0,:] 

    w_lm_b   = planet_data_b[1:,0]
    ImK2_b   = planet_data_b[1:,1]
    ReK2_b   = planet_data_b[1:,2]
    size_b   = np.size(w_lm_b)

    #===============================================================================================
    planet_data_c = np.loadtxt('Results_Trappist_1c_07.5_100_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_c, planet_radius_c, planet_gyration_radius_squared_c = planet_data_c[0,:] 

    w_lm_c   = planet_data_c[1:,0]
    ImK2_c   = planet_data_c[1:,1]
    ReK2_c   = planet_data_c[1:,2]
    size_c   = np.size(w_lm_c)

    #===============================================================================================
    planet_data_d = np.loadtxt('Results_Trappist_1d_19.5_100_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_d, planet_radius_d, planet_gyration_radius_squared_d = planet_data_d[0,:] 

    w_lm_d   = planet_data_d[1:,0]
    ImK2_d   = planet_data_d[1:,1]
    ReK2_d   = planet_data_d[1:,2]
    size_d   = np.size(w_lm_d)

    #===============================================================================================
    # planet_data_e = np.loadtxt('Results_Trappist_1e_00.0_138_freq_Imk2_posidonius.txt',comments='#')
    planet_data_e = np.loadtxt('Results_Trappist_1e_05.0_150_freq_Imk2_posidonius.txt',comments='#') 
    planet_mass_e, planet_radius_e, planet_gyration_radius_squared_e = planet_data_e[0,:] 

    w_lm_e   = planet_data_e[1:,0]
    ImK2_e   = planet_data_e[1:,1]
    ReK2_e   = planet_data_e[1:,2]
    size_e   = np.size(w_lm_e)

    #===============================================================================================
    planet_data_f = np.loadtxt('Results_Trappist_1f_11.0_100_freq_Imk2_posidonius.txt',comments='#')
    # planet_data_f = np.loadtxt('Results_Trappist_1f_14.5_150_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_f, planet_radius_f, planet_gyration_radius_squared_f = planet_data_f[0,:] 

    w_lm_f   = planet_data_f[1:,0]
    ImK2_f   = planet_data_f[1:,1]
    ReK2_f   = planet_data_f[1:,2]
    size_f   = np.size(w_lm_f)

    #===============================================================================================
    planet_data_g = np.loadtxt('Results_Trappist_1g_19.0_100_freq_Imk2_posidonius.txt',comments='#')
    # planet_data_g = np.loadtxt('Results_Trappist_1g_22.5_150_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_g, planet_radius_g, planet_gyration_radius_squared_g = planet_data_g[0,:] 

    w_lm_g   = planet_data_g[1:,0]
    ImK2_g   = planet_data_g[1:,1]
    ReK2_g   = planet_data_g[1:,2]
    size_g   = np.size(w_lm_g)

    #===============================================================================================
    planet_data_h = np.loadtxt('Results_Trappist_1h_10.0_100_freq_Imk2_posidonius.txt',comments='#')
    # planet_data_h = np.loadtxt('Results_Trappist_1h_13.5_150_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_h, planet_radius_h, planet_gyration_radius_squared_h = planet_data_h[0,:] 

    w_lm_h   = planet_data_h[1:,0]
    ImK2_h   = planet_data_h[1:,1]
    ReK2_h   = planet_data_h[1:,2]
    size_h   = np.size(w_lm_h)

    #===============================================================================================
    planet_data_venus = np.loadtxt('Results_V1_Armann2_atm_freq_Imk2_posidonius.txt',comments='#')
    # planet_data_h = np.loadtxt('Results_Trappist_1h_13.5_150_freq_Imk2_posidonius.txt',comments='#')
    planet_mass_venus, planet_radius_venus, planet_gyration_radius_squared_venus = planet_data_venus[0,:] 

    w_lm_venus   = planet_data_venus[1:,0]
    ImK2_venus   = planet_data_venus[1:,1]
    ReK2_venus   = planet_data_venus[1:,2]
    size_venus   = np.size(w_lm_venus)

    # ############################################################################
    # # Radiuses in R_EARTH
    # # planet_radiuses = (1.086e0, 1.056e0, 0.772e0, 0.918e0, 1.045e0, 1.127e0, 0.755e0)
    # planet_radiuses = (planet_radius_b, planet_radius_c, planet_radius_d, planet_radius_e, planet_radius_f, planet_radius_g, planet_radius_h)
    # # Masses in M_SUN
    # # planet_masses = (2.5843686e-06, 4.1957984e-06, 1.2465778e-06, 1.8850689e-06, 2.0674949e-06, 4.0741811e-06, 1.2465778e-06)
    # planet_masses = (planet_mass_b, planet_mass_c, planet_mass_d, planet_mass_e, planet_mass_f, planet_mass_g, planet_mass_h)
    # # Semi-major axis in AU
    # # planet_a = (0.01111, 0.01521, 0.02144, 0.02817, 0.0371, 0.0451, 0.0596)
    # planet_a = (0.01154775, 0.01581512, 0.02228038, 0.02928285, 0.03853361, 0.04687692, 0.06193488)
    # # Inclination in degrees
    # # planet_i = (0.3500000, 0.330000, 0.250000, 0.140000, 0.320000, 0.290000, 0.130000)
    # planet_i = (0., 0., 0., 0., 0., 0., 0.)
    # # Mean anomaly in degrees
    # # planet_l = (323.732652895, 96.4925777097, 111.770368348, 165.724187804, 254.117367005, 161.020362506, 134.724813585)
    # planet_l = (0., 0., 0., 0., 0., 0., 0.)

    # planet_e = (0., 0., 0., 0., 0., 0., 0.)
    # # planet_e = (0.00622, 0.00654, 0.0837, 0.0051, 0.01007, 0.00208, 0.00567)
    # planet_w_lm = (w_lm_b, w_lm_c, w_lm_d, w_lm_e, w_lm_f, w_lm_g, w_lm_h)
    # planet_Imk2 = (ImK2_b, ImK2_c, ImK2_d, ImK2_e, ImK2_f, ImK2_g, ImK2_h)
    # planet_Rek2 = (ReK2_b, ReK2_c, ReK2_d, ReK2_e, ReK2_f, ReK2_g, ReK2_h)
    # planet_nm_data = (size_b, size_c, size_d, size_e, size_f, size_g, size_h)
    # planet_gyration_radius_squared = (planet_gyration_radius_squared_b, planet_gyration_radius_squared_c, planet_gyration_radius_squared_d, planet_gyration_radius_squared_e, planet_gyration_radius_squared_f, planet_gyration_radius_squared_g, planet_gyration_radius_squared_h)
    # ############################################################################
    # Radiuses in R_EARTH
    # planet_radiuses = (1.086e0, 1.056e0, 0.772e0, 0.918e0, 1.045e0, 1.127e0, 0.755e0)
    planet_radiuses = (planet_radius_e)
    # Masses in M_SUN
    # planet_masses = (2.5843686e-06, 4.1957984e-06, 1.2465778e-06, 1.8850689e-06, 2.0674949e-06, 4.0741811e-06, 1.2465778e-06)
    planet_masses = (planet_mass_e)
    # Semi-major axis in AU
    # planet_a = (0.01111, 0.01521, 0.02144, 0.02817, 0.0371, 0.0451, 0.0596)
    planet_a = (0.02928285)
    # Inclination in degrees
    # planet_i = (0.3500000, 0.330000, 0.250000, 0.140000, 0.320000, 0.290000, 0.130000)
    planet_i = (0.)
    # Mean anomaly in degrees
    # planet_l = (323.732652895, 96.4925777097, 111.770368348, 165.724187804, 254.117367005, 161.020362506, 134.724813585)
    planet_l = (0.)

    # planet_e = (0.00622, 0.00654, 0.0837, 0.0051, 0.01007, 0.00208, 0.00567)
    planet_e = (0.0051,)

    planet_w_lm = (w_lm_e)
    planet_Imk2 = (ImK2_e)
    planet_Rek2 = (ReK2_e)
    planet_nm_data = (size_e)
    planet_gyration_radius_squared = (planet_gyration_radius_squared_e)

    # for r, m, a, i, l, e, w_lmpq, ImK2, ReK2, size, gyration_radius_2 in zip(planet_radiuses, planet_masses, planet_a, planet_i, planet_l, planet_e, planet_w_lm, planet_Imk2, planet_Rek2, planet_nm_data, planet_gyration_radius_squared):
    Mearth = 5.9736E+24 #Earth Mass (kg)
    Rearth = 6378.1E+3    #Earth radius (m)
    Mvenus = 4.8685E+24
    Rvenus = 6052E+3
    Mjup = 1.898E+27 # Jupiter mass (kg)
    Rjup = 6.991E+07 # Jupiter radius (m)
    #===========================T1e
    # r = planet_radius_e
    # m = planet_mass_e
    # a = 0.02928285
    # i = 0.
    # l = 0.
    # e = 0.0051#planet_e[0] #0.0#0.07
    # w_lmpq = w_lm_e
    # ImK2 = ImK2_e
    # ReK2 = ReK2_e
    # size = size_e
    # gyration_radius_2 = planet_gyration_radius_squared_e
    #==========================T1b
    # r = planet_radius_b
    # m = planet_mass_b
    # a = 0.01111
    # i = 0.
    # l = 0.
    # e = 0.0#1 #0.00622
    # w_lmpq = w_lm_b
    # ImK2 = ImK2_b
    # ReK2 = ReK2_b
    # size = size_b
    # gyration_radius_2 = planet_gyration_radius_squared_b
    # #=========================

    #==========================Venus
    r = Rvenus # planet_radius_venus
    m = Mvenus # planet_mass_venus
    # a = 0.0723 # AU
    # e = 0.05#0.802 #0.00622
    w_lmpq = w_lm_venus
    ImK2 = ImK2_venus
    ReK2 = ReK2_venus
    size = size_venus
    gyration_radius_2 = planet_gyration_radius_squared_venus
    # ------------
    planet_mass = m / posidonius.constants.M_SUN # Solar masses (3.0e-6 solar masses = 1 earth mass)
    planet_radius_factor = r # R Earth
    # planet_radius = planet_radius_factor * posidonius.constants.R_EARTH / 6.3781e6
    planet_radius = 0.950*6.3781E+6 / posidonius.constants.AU
    # planet_radius_of_gyration = 5.75e-01
    planet_radius_of_gyration = np.sqrt(gyration_radius_2)


    # --------------------------------------------------------
    # --- Specify initial position and velocity for a stable orbit
    # --- Keplerian orbital elements, in the `asteroidal' format of Mercury code
    a = 0.0723                                    # semi-major axis (in AU)
    e = 0.901                                      # eccentricity
    i = 300. * posidonius.constants.DEG2RAD;        # inclination (degrees to radiant)
    p = 20. * posidonius.constants.DEG2RAD;     # argument of pericentre (degrees to radiant) !! is != p is the longitude of perihelion in tools
    n = 30. * posidonius.constants.DEG2RAD;      # longitude of the ascending node (degrees to radiant)
    l = 2. * posidonius.constants.DEG2RAD;        # mean anomaly (degrees to radiant)
    p = (p + n);                                 # Convert argument pericentre to longitude of perihelion !!
    q = a * (1.0 - e);                           # perihelion distance
    # 
    # 
    print("Verifications: \n SemiMajor axis {} AU, eccentricity {} radius {} meter \n".format(a,e,r))
    print(" \t planet mass {} radius {} \n\t star mass {} ".format(planet_mass,planet_radius,star_mass))
    #
    #
    # --- Calculate planet position and velocity from the keplerian orbital elements
    print("\n Input orbital elements \n\t a={} q={} e={} \n\t i={}rad ={}deg \n\t p={}rad ={}deg (or ={}deg)  \n\t n={}rad ={}deg  \n\t l={}rad ={}deg \n" \
        .format(a,q,e,i,i/posidonius.constants.DEG2RAD,p,p/posidonius.constants.DEG2RAD,p/posidonius.constants.DEG2RAD-360,n,n/posidonius.constants.DEG2RAD,l,l/posidonius.constants.DEG2RAD))
    # print("\n Input orbital elements \n\t a={} q={} e={} \n\t i={}rad ={}deg \n\t p={}rad ={}deg={}deg  \n\t n={}rad ={}deg  \n\t p-n={}rad ={}deg or {}  \n\t l={}rad ={}deg \n" \
    #     .format(a,q,e,i,i/posidonius.constants.DEG2RAD,p,p/posidonius.constants.DEG2RAD,p/posidonius.constants.DEG2RAD-360,n,n/posidonius.constants.DEG2RAD,p-n,(p-n)/posidonius.constants.DEG2RAD,360-(p-n)/posidonius.constants.DEG2RAD,l,l/posidonius.constants.DEG2RAD))
    
    # --- Calculate the position and velocity vectors
    print(" -------------------------- ")
    planet_position, planet_velocity = posidonius.tools.calculate_cartesian_coordinates(planet_mass, q, e, i, p, n, l, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
    print("\n Planet Position {} Velocity {} " \
        .format(planet_position.get(), planet_velocity.get()))

    # --- Calculate Old Keplerian orbital elements
    out_a, out_q, out_e, out_i, out_p, out_n, out_l = posidonius.tools.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
    print("\n Old Module orbital elements \n\t a={} q={} e={}  \n\t i={}rad ={}deg \n\t p={}rad ={}deg  \n\t n={}rad ={}deg \n\t l={}rad ={}deg" \
        .format(out_a, out_q, out_e, out_i, out_i/posidonius.constants.DEG2RAD, out_p, out_p/posidonius.constants.DEG2RAD, out_n,out_n/posidonius.constants.DEG2RAD, out_l, out_l/posidonius.constants.DEG2RAD))
    
    # --- Local calculation
    Local_calculation = False
    if Local_calculation:
        print(" \n --- My calculation")
        x = planet_position.x()
        y = planet_position.y()
        z = planet_position.z()
        u = planet_velocity.x()
        v = planet_velocity.y()
        w = planet_velocity.z()
        # --- Eccentricity vector
        ex,ey,ez = posidonius.tools.calculate_eccentricity_vector(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
        print("\tEccentricity vector {} {} {} ".format(ex,ey,ez))
        print("\tecc vector norm", np.sqrt(np.power(ex,2) +np.power(ey,2) +np.power(ez,2) ))
        e_scal_r = ex*x + ey*y + ez*z
        # --- Angular momentum
        hx = y * w  -  z * v
        hy = z * u  -  x * w
        hz = x * v  -  y * u
        print("\tAngular momentum {} {} {}".format(hx,hy,hz))
        # --- Some local variable
        h2 = np.power(hx, 2) + np.power(hy, 2) + np.power(hz, 2)
        r_cal = np.sqrt( np.power(x,2) +np.power(y,2) +np.power(z,2) )
        v_cal_2 = ( np.power(u,2) +np.power(v,2) +np.power(w,2) )
        r_scal_v = x*u +y*v +z*w
        reduced_mass = posidonius.constants.G*(planet_mass+star_mass)
        
        # --- Semi major axis
        a_cal = 1./(2./r_cal -v_cal_2/reduced_mass)
        print("\ta Calculated", a_cal)
        
        # --- eccentricity
        e_sqare = 1. -h2/(reduced_mass*a_cal)
        if e_sqare < np.power(3.0e-8,2):
            e_cal = 0. #np.sqrt(1. -h2/(reduced_mass*a_cal))
        else:
            e_cal = np.sqrt(1. -h2/(reduced_mass*a_cal))
        print("\te Calculated", e_cal)
        q_cal = a_cal*(1-e_cal)

        # --- inclination
        # i_cal = np.arccos(hz/np.sqrt(h2))
        # print("\tiCalculated", i_cal, i_cal/posidonius.constants.DEG2RAD)
        # sign = math.copysign(1,hx_rot)
        # i_cal = sign*np.arccos(hz/np.sqrt(h2))
        # # print("SIGNE ", math.copysign(1,hx_rot), "-1 ", np.sign(-1), "0", np.sign(0), "1", np.sign(1))
        # print("\tiReCalculated", i_cal, i_cal/posidonius.constants.DEG2RAD)
        i_cal = np.arccos(hz/np.sqrt(h2))
        if hy>0.:
            i_cal = 2*np.pi - i_cal
        print("\ti ReCalculated", i_cal, i_cal/posidonius.constants.DEG2RAD)

        # --- The longitude of ascending node 
        if i_cal==0: # need an inclination to be correctly determined
            long_asc_node_cal = 0.
        else:
            long_asc_node_cal = np.arctan(hx/-hy)
            # if hy<0.: 
            #     long_asc_node_cal = np.pi + long_asc_node_cal
            # if (hy>0. and hx>0.):
            #     long_asc_node_cal = 2*np.pi + long_asc_node_cal
            print("\tlong asc node cal", long_asc_node_cal, long_asc_node_cal/posidonius.constants.DEG2RAD)


        # --- The longitude of ascending node 
        # long_asc_node_cal = -np.arctan2(hx, hy)
        hx_rot = hx*np.cos(long_asc_node_cal) + hy*np.sin(long_asc_node_cal)
        hy_rot = - hx*np.sin(long_asc_node_cal) + hy*np.cos(long_asc_node_cal)
        hz_rot = hz
        print("\t Angular momentum Rotate{} {} {}".format(hx_rot,hy_rot,hz_rot))

        # --- argument of pericenter
        # Through the Matrix R3 rotation passage of the longitude of ascending node
        ex_rot = ex*np.cos(long_asc_node_cal) + ey*np.sin(long_asc_node_cal)
        ey_rot = - ex*np.sin(long_asc_node_cal) + ey*np.cos(long_asc_node_cal)
        ez_rot = ez
        print("\teccentricity vector tho rotation matrix {} {} {} ".format(ex_rot, ey_rot, ez_rot))
        if e_cal==0.:
            arg_pericentre = 0.
        else:    
            arg_pericentre = np.arctan2(ey_rot,ex_rot)
        print("\targument pericenter ", arg_pericentre, arg_pericentre/posidonius.constants.DEG2RAD)

        # --- The longitude of pericentre 
        long_pericentre = long_asc_node_cal + arg_pericentre
        # long_pericentre = np.arctan2(ey,ex)
        print("\tlongitude pericenter ", long_pericentre, long_pericentre/posidonius.constants.DEG2RAD)
        # long_pericentre = np.arctan(ex/ey)
        # print("\tlongitude pericenter ", long_pericentre, long_pericentre/posidonius.constants.DEG2RAD)

        # --- The true anomaly
        small_ecentricity = e_cal < 1.0e-7
        if small_ecentricity:
            x_rot = x*np.cos(long_asc_node_cal) + y*np.sin(long_asc_node_cal)
            y_rot = - x*np.sin(long_asc_node_cal) + y*np.cos(long_asc_node_cal)
            z_rot = z
            cosf_cal = x_rot/r_cal
            f = np.arccos(cosf_cal)
            print("\tTrue long in small ecc {} rad and {} deg".format(f, f/posidonius.constants.DEG2RAD))
        else:
            cosf_cal = e_scal_r / (e_cal*r_cal)
            # print(" x {} y {} z {} =  r {} ".format(x,y,z,r_cal))
            # print("e scal r {} e_cal {}".format( e_scal_r, e_cal))
            # print(" cosf {} ".format(cosf_cal))
            if r_scal_v<0:
                print("\tr scal v is neg")
                f = 2*np.pi -np.arccos(cosf_cal)
                # print(np.arccos(cosf_cal))
            else :
                f = np.arccos(cosf_cal)
            print("\tTrue ano cal {} rad and {} deg".format( f, f/posidonius.constants.DEG2RAD))
            cosf_cal = (1/e_cal)*( (a_cal*(1-e_cal*e_cal)/r_cal) -1.)
            f = np.arccos(cosf_cal) 
            print("\tTrue ano cal {} rad and {} deg".format( f, f/posidonius.constants.DEG2RAD))
        
        # --- mean anomaly 
        if small_ecentricity: # if small ecentricity :: mean_anomaly = true_longitude = true_anomaly + long_pericentre
            mean_ano_cal = f #+ long_pericentre
        else:
            cos_E = (1/e_cal)*(1.-(r_cal/a_cal))
            mean_ano_cal = np.arccos(cos_E) - e_cal*np.sin(np.arccos(cos_E) )
        print("\tMean anomaly cal", mean_ano_cal, mean_ano_cal/posidonius.constants.DEG2RAD)
        
        
        # --- Output
        p_cal = arg_pericentre + long_asc_node_cal
        n_cal = long_asc_node_cal
        l_cal = mean_ano_cal
        print("\n    --- Calculated orbital elements \n\t a={} q={} e={}  \n\t i={}rad ={}deg \n\t p={}rad ={}deg  \n\t n={}rad ={}deg \n\t l={}rad ={}deg" \
            .format(a_cal, q_cal, e_cal, i_cal, i_cal/posidonius.constants.DEG2RAD, p_cal, p_cal/posidonius.constants.DEG2RAD, n_cal, n_cal/posidonius.constants.DEG2RAD, l_cal, l_cal/posidonius.constants.DEG2RAD))

    # --- My module of calculate keplerian orbital elements
    out_a_cal, out_q_cal, out_e_cal, out_i_cal, out_p_cal, out_n_cal, out_l_cal = posidonius.tools.calculate_keplerian_orbital_elements_2(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
    print("\n New Module orbital elements \n\t a={} q={} e={}  \n\t i={}rad ={}deg \n\t p={}rad ={}deg  \n\t n={}rad ={}deg \n\t l={}rad ={}deg" \
        .format(out_a_cal, out_q_cal, out_e_cal, out_i_cal, out_i_cal/posidonius.constants.DEG2RAD, out_p_cal, out_p_cal/posidonius.constants.DEG2RAD, out_n_cal, out_n_cal/posidonius.constants.DEG2RAD, out_l_cal, out_l_cal/posidonius.constants.DEG2RAD))
    
    # --- 
    print(" -------------------------- ")
    out_planet_position, out_planet_velocity = posidonius.tools.calculate_cartesian_coordinates(planet_mass, out_q, out_e, out_i, out_p, out_n, out_l, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
    print("\n1_ Planet Position {} Velocity {} " \
        .format(out_planet_position.get(), out_planet_velocity.get()))

    print(" -------------------------- ")
    out_planet_position, out_planet_velocity = posidonius.tools.calculate_cartesian_coordinates(planet_mass, out_q_cal, out_e_cal, out_i_cal, out_p_cal, out_n_cal, out_l_cal, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
    print("\n2_ Planet Position {} Velocity {} " \
        .format(out_planet_position.get(), out_planet_velocity.get()))

    # # --- My module of calculate keplerian orbital elements
    # print(" -------------------------- ")
    # out_a_cal, out_q_cal, out_e_cal, out_i_cal, out_p_cal, out_n_cal, out_l_cal = posidonius.tools.calculate_keplerian_orbital_elements_2(planet_mass, out_planet_position, out_planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
    # print("\n New Module orbital elements \n\t a={} q={} e={}  \n\t i={}rad ={}deg \n\t p={}rad ={}deg  \n\t n={}rad ={}deg \n\t l={}rad ={}deg" \
    #     .format(out_a_cal, out_q_cal, out_e_cal, out_i_cal, out_i_cal/posidonius.constants.DEG2RAD, out_p_cal, out_p_cal/posidonius.constants.DEG2RAD, out_n_cal, out_n_cal/posidonius.constants.DEG2RAD, out_l_cal, out_l_cal/posidonius.constants.DEG2RAD))

    # print(" -------------------------- ")
    # out_planet_position, out_planet_velocity = posidonius.tools.calculate_cartesian_coordinates(planet_mass, out_q_cal, out_e_cal, out_i_cal, out_p_cal, out_n_cal, out_l_cal, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
    # print("\n2_ Planet Position {} Velocity {} " \
    #     .format(out_planet_position.get(), out_planet_velocity.get()))

    print(" -------------------------- ")
    
    # --- Initialization of planetary spin
    planet_obliquity = 0. #1.0e-4 # rad
    
    # --- Pseudo-synchronization period
    planet_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(a, e, star_mass, planet_mass) # days
    # planet_angular_frequency = posidonius.constants.TWO_PI/(planet_pseudo_synchronization_period) # days^-1
    # planet_angular_frequency = posidonius.constants.TWO_PI/(15./24.) # days^-1
    # planet_angular_frequency = posidonius.constants.TWO_PI/planet_keplerian_orbital_elements[7]# days^-1
    orbital_period  = 2.*np.pi*(np.sqrt( (a*posidonius.constants.AU)**3. / ( posidonius.constants.G_SI*( star_mass*posidonius.constants.M_SUN + planet_mass*posidonius.constants.M_SUN ) ) ) )/86400.
    # print("Orbtial period is ",orbital_period)
    rotation_period = orbital_period/1.5
    planet_angular_frequency = posidonius.constants.TWO_PI/(rotation_period)


    planet_inclination = i
    planet_spin = posidonius.calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity)

    k2pdelta = 2.465278e-3 # Terrestrial planets (no gas)
    planet_tides_parameters = {
        "dissipation_factor_scale": 1.0,
        "dissipation_factor": 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet_radius, 5)),
        "love_number": 1.0 * 0.299 #0.307#0.299,
    }

    Tab_size = 32
    Tab_w_lmpq = np.zeros((Tab_size,Tab_size))
    Tab_ImK2 = np.zeros((Tab_size,Tab_size))
    Tab_ReK2 = np.zeros((Tab_size,Tab_size))
    # print('\n Size table ', Tab_size, ' nombre de case ', Tab_size*Tab_size)
    # print('Nombre de points', size)
    # print('stocke donc ', int(size/Tab_size) , 'lignes de 32 cases, plus  ', size-32*int(size/Tab_size), ' cases soit ',32*int(size/Tab_size)+(size-32*int(size/Tab_size)) )
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

    print("\t")
    # planet_tides = posidonius.effects.tides.CentralBody(planet_tides_parameters)
    # planet_tides = posidonius.effects.tides.OrbitingBody(planet_tides_parameters)
    planet_tides = posidonius.effects.tides.KaulaCoplanarOrbitingBody(planet_tides_parameters_love_number_kaula)
    # planet_tides = posidonius.effects.tides.Disabled()
    #
    # planet_rotational_flattening_parameters = {"love_number": 0.9532}
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
    # universe.write(filename, integrator="IAS15")
    #universe.write(filename, integrator="LeapFrog")


