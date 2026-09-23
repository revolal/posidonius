import os
import numpy as np
import pandas as pd
import matplotlib as mpl
# mpl.use('Agg')
mpl.rcParams['agg.path.chunksize'] = 10000
import matplotlib.pyplot as plt
import argparse
import posidonius
import json
from statistics import mean, median
import matplotlib.ticker as ticker
from matplotlib.ticker import (MultipleLocator,
                               FormatStrFormatter,
                               AutoMinorLocator)
# from sklearn.linear_model import LinearRegression
# from sklearn import datasets, linear_model
from scipy.integrate import cumulative_trapezoid
from pylaplace import LaplaceCoefficient
from numpy import diff

def calculate_eccentricity_vector(gm, x, y, z, u, v, w):

    # // --- Input --- //
    # x = position[0]
    # y = position[1]
    # z = position[2]
    # u = velocity[0]
    # v = velocity[1]
    # w = velocity[2]
    
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

def calculate_tidal_excitation_frequency_mode_sigma_2mpq(m, p, q, spin_angle, mean_motion):
    return (2.0 -2.0*p +q)*mean_motion -m*spin_angle

def calculate_eccentricity_function_g_20q(eccentricity):
    eccentricity_function = [0.0] * 15

    ecc = eccentricity
    ecc_2 = ecc * ecc
    ecc_3 = ecc_2 * ecc
    ecc_4 = ecc_3 * ecc
    ecc_5 = ecc_4 * ecc
    ecc_6 = ecc_5 * ecc
    ecc_7 = ecc_6 * ecc

    eccentricity_function[0] = (15625. / 129024.) * ecc_7
    eccentricity_function[1] = (4. / 45.) * ecc_6
    eccentricity_function[2] = (81. / 1280.) * ecc_5 + (81. / 2048.) * ecc_7
    eccentricity_function[3] = (1. / 24.) * ecc_4 + (7. / 240.) * ecc_6
    eccentricity_function[4] = (1. / 48.) * ecc_3 + (11. / 768.) * ecc_5 + (313. / 30720.) * ecc_7
    eccentricity_function[5] = 0.0 * ecc
    eccentricity_function[6] = -0.5 * ecc + (1. / 16.) * ecc_3 - (5. / 384.) * ecc_5 - (143. / 18432.) * ecc_7
    eccentricity_function[7] = 1.0 - (5. / 2.) * ecc_2 + (13. / 16.) * ecc_4 - (35. / 288.) * ecc_6
    eccentricity_function[8] = (7. / 2.) * ecc - (123. / 16.) * ecc_3 + (489. / 128.) * ecc_5 - (1763. / 2048.) * ecc_7
    eccentricity_function[9] = (17. / 2.) * ecc_2 - (115. / 16.) * ecc_4 + (601. / 48.) * ecc_6
    eccentricity_function[10] = (845. / 48.) * ecc_3 - (32525. / 768.) * ecc_5 + (208225. / 6144.) * ecc_7
    eccentricity_function[11] = (533. / 16.) * ecc_4 - (13827. / 160.) * ecc_6
    eccentricity_function[12] = (228347. / 3840.) * ecc_5 - (3071075. / 18432.) * ecc_7
    eccentricity_function[13] = (73369. / 720.) * ecc_6
    eccentricity_function[14] = (12144273. / 71680.) * ecc_7

    return eccentricity_function

def calculate_eccentricity_function_g_21q(eccentricity):
    eccentricity_function = [0.0] * 15

    ecc = eccentricity
    ecc_2 = ecc * ecc
    ecc_3 = ecc_2 * ecc
    ecc_4 = ecc_3 * ecc
    ecc_5 = ecc_4 * ecc
    ecc_6 = ecc_5 * ecc
    ecc_7 = ecc_6 * ecc

    eccentricity_function[0] = (432091. / 30720.) * ecc_7
    eccentricity_function[1] = (3167. / 320.) * ecc_6
    eccentricity_function[2] = (1773. / 256.) * ecc_5 - (4987. / 6144.) * ecc_7
    eccentricity_function[3] = (77. / 16.) * ecc_4 + (129. / 160.) * ecc_6
    eccentricity_function[4] = (53. / 16.) * ecc_3 + (393. / 256.) * ecc_5 + (24753. / 10240.) * ecc_7
    eccentricity_function[5] = (9. / 4.) * ecc_2 + (7. / 4.) * ecc_4 + (141. / 64.) * ecc_6
    eccentricity_function[6] = (3. / 2.) * ecc + (27. / 16.) * ecc_3 + (261. / 128.) * ecc_5 + (14309. / 6144.) * ecc_7
    eccentricity_function[7] = (1. - ecc_2) ** (-3. / 2.)
    eccentricity_function[8] = eccentricity_function[6]
    eccentricity_function[9] = eccentricity_function[5]
    eccentricity_function[10] = eccentricity_function[4]
    eccentricity_function[11] = eccentricity_function[3]
    eccentricity_function[12] = eccentricity_function[2]
    eccentricity_function[13] = eccentricity_function[1]
    eccentricity_function[14] = eccentricity_function[0]

    return eccentricity_function
1000
def calculate_eccentricity_function_g_22q(eccentricity):
    eccentricity_function = [0.0] * 15

    ecc = eccentricity
    ecc_2 = ecc * ecc
    ecc_3 = ecc_2 * ecc
    ecc_4 = ecc_3 * ecc
    ecc_5 = ecc_4 * ecc
    ecc_6 = ecc_5 * ecc
    ecc_7 = ecc_6 * ecc

    eccentricity_function[14] = (15625. / 129024.) * ecc_7
    eccentricity_function[13] = (4. / 45.) * ecc_6
    eccentricity_function[12] = (81. / 1280.) * ecc_5 + (81. / 2048.) * ecc_7
    eccentricity_function[11] = (1. / 24.) * ecc_4 + (7. / 240.) * ecc_6
    eccentricity_function[10] = (1. / 48.) * ecc_3 + (11. / 768.) * ecc_5 + (313. / 30720.) * ecc_7
    eccentricity_function[9] = 0.0
    eccentricity_function[8] = -0.5 * ecc + (1. / 16.) * ecc_3 - (5. / 384.) * ecc_5 - (143. / 18432.) * ecc_7
    eccentricity_function[7] = 1.0 - (5. / 2.) * ecc_2 + (13. / 16.) * ecc_4 - (35. / 288.) * ecc_6
    eccentricity_function[6] = (7. / 2.) * ecc - (123. / 16.) * ecc_3 + (489. / 128.) * ecc_5 - (1763. / 2048.) * ecc_7
    eccentricity_function[5] = (17. / 2.) * ecc_2 - (115. / 16.) * ecc_4 + (601. / 48.) * ecc_6
    eccentricity_function[4] = (845. / 48.) * ecc_3 - (32525. / 768.) * ecc_5 + (208225. / 6144.) * ecc_7
    eccentricity_function[3] = (533. / 16.) * ecc_4 - (13827. / 160.) * ecc_6
    eccentricity_function[2] = (228347. / 3840.) * ecc_5 - (3071075. / 18432.) * ecc_7
    eccentricity_function[1] = (73369. / 720.) * ecc_6
    eccentricity_function[0] = (12144273. / 71680.) * ecc_7

    return eccentricity_function



if __name__ == "__main__":
    print("\n \t -- START PLOT --")
    parser = argparse.ArgumentParser()
    parser.add_argument('start_case_filename', action='store', help='Filename with the initial conditions of the simulation (e.g., universe_integrator.json)')
    parser.add_argument('historic_snapshot_filename', action='store', help='Filename with the historic snapshots of the simulation (e.g., universe_integrator_history.bin)')

    args = parser.parse_args()
    print("\t Load json file")
    universe_integrator_json = json.load(open(args.start_case_filename, "r"))
    
    print("\t Read filename")
    filename = args.historic_snapshot_filename
    n_particles, data = posidonius.analysis.history.read(filename)
    most_massive_particle_index = universe_integrator_json['universe']['hosts']['index']['most_massive']
    print("\t Transforming positions/velocities to heliocentric coordinates using the most masssive particle at index '{}'...".format(most_massive_particle_index))
    star_data, planets_data, planets_keys = posidonius.analysis.history.classify(n_particles, data, reference_particle_index=most_massive_particle_index, discard_first_years=False, first_years=50., discard_last_years=False, last_years=300.)
    # star_data, planets_data, planets_keys = posidonius.analysis.history.classify(n_particles, data, reference_particle_index=most_massive_particle_index, discard_first_hundred_years=False)
    # star_data = star_data[star_data['current_time'] >= 10. ]
    for key in planets_keys:
        # planets_data[key] = planets_data[key][planets_data[key]['current_time'] >= 10. ]
        print('Time : ', planets_data[key]['current_time'][0] ,  planets_data[key]['current_time'][ len(planets_data[key]['current_time'])-1 ], ' | Ndata : ', len(planets_data[key]['current_time']))
    star_mass = star_data['mass'][0]


    ### Select one every two data points
    #one_every_two = np.arange(len(star_data)) % 2 == 1
    #star_data = star_data[one_every_two]
    #for key in planets_keys:
        #planets_data[key] = planets_data[key][one_every_two]

    #-------------------------------------------------------------------------------
    # Main
    #-------------------------------------------------------------------------------
    print("Computing keplerian and supplementary values...")

    ################################################################################
    ## Star
    ################################################################################
    star_norm_spin = np.sqrt(np.power(star_data['spin_x'], 2) + np.power(star_data['spin_y'], 2) + np.power(star_data['spin_z'], 2))
    star_rotation_period = 2*np.pi / star_norm_spin


    ################################################################################
    plot_1      = False # Plot AM conservation
    plot_2      = False # Plot data
    plot_2b     = True # Plot data reduced
    plot_20     = False # Plot data
    plot_21     = False # Plot positions
    plot_22     = False # Plot precession periods
    plot_23     = False # Print data 
    plot_3e     = False # plot ecc data
    plot_3ebis  = True # plot ecc data
    plot_3s     = False # plot spin data
    plot_30     = False # plot spin / mean motion
    plot_31     = False # plot sub-stellar point 
    plot_31b    = False # plot sub-stellar point for planet b
    plot_31all  = False # plot sub-stellar point for each planet
    plot_32     = False # plot sub-stellar point 2 calculs
    plot_33     = False # Check angles
    plot_33_1p  = False # Check true anomaly
    plot_33_meanano = False
    plot_34     = False # Plot Spin each planet 
    plot_4      = False # compare spin mean motion
    plot_41     = False # compare mean motion 
    plot_42     = False # compare mean_anomaly true anomaly
    plot_43     = False # compare eccentricities
    plot_5      = False # compare orbital period with observed orbital period from Agol+21
    plot_6      = False # compare sma 
    plot_7      = False # compare period ratio pi+1/pi with observed 
    plot_70     = False # compare period ratio offset with observed - END value
    plot_71     = False # compare period ratio offset with observed - Evolution each planets
    plot_72     = False # compare period ratio offset with observed - Evolution 
    plot_73     = False # Fig-1 Charalambous et al. 2023 | Departure from 2-body MMRs
    plot_8      = False # planet period ratio chain
    plot_9      = False # Resonance angle 3 bodies
    plot_10     = False # Resonance angle 2 bodies
    plot_101    = False # Resonance angle 2 planets 2 bodies
    plot_102    = True # Look for resonance angle 2 planets by pairs 
    plot_11     = False # Resonance modes lmpq frequency angle two bodies
    plot_12     = False # Print Hansen coeff
    history_output = True

    ################################################################################
    ## Planets
    ################################################################################
    total_planets_angular_momentum_x = 0.
    total_planets_angular_momentum_y = 0.
    total_planets_angular_momentum_z = 0.
    total_planets_angular_momentum   = 0.
    total_planets_orbital_angular_momentum = 0.
    planets_computed_data = {}
    print('## CONTROL CARTESIAN COORDINATES:')
    # print("Time 0")
    # for key in planets_keys:
    #     planet_data = planets_data[key]
    #     print(planet_data['current_time'][0])
    # print("Star")
    # print("Position : ", star_data["position_x"],star_data["position_z"],star_data["position_z"] )
    # print("Velocity : ", star_data["velocity_x"],star_data["velocity_z"],star_data["velocity_z"] )
    # print("Planets")
    # print("positions")
    # for key in planets_keys:
    #     planet_data = planets_data[key]
    #     print(planet_data['position_x'][0], planet_data['position_y'][0], planet_data['position_z'][0])
    # print("velocity")
    # for key in planets_keys:
    #     planet_data = planets_data[key]
    #     print(planet_data['velocity_x'][0], planet_data['velocity_y'][0], planet_data['velocity_z'][0])
    print("> Check final orbtial parameters")
    for key in planets_keys:
        planet_data = planets_data[key]
        list = ('semi-major_axis', 'eccentricity', 'inclination', 'longitude_of_perihelion', 'longitude_of_ascending_node', 'mean_anomaly', 'true_anomaly')
        for l in list :
            print(l, planet_data[l][-1])

    for key in planets_keys:
        print("\n| ===== Compute for planet ", key, " ===== | ")
        planet_data = planets_data[key]
        planet_computed_data = {}

        planet_norm_spin = np.sqrt(np.power(planet_data['spin_x'], 2) + np.power(planet_data['spin_y'], 2) + np.power(planet_data['spin_z'], 2))
        planet_rotation_period = 2*np.pi / planet_norm_spin
        planet_computed_data['planet_rotation_period'] = planet_rotation_period

        

        ## Planet obliquity
        # Calculation of orbital angular momentum (without mass and in AU^2/day)
        orbital_angular_momentum_x = planet_data['position_y'] * planet_data['velocity_z'] - planet_data['position_z'] * planet_data['velocity_y'];
        orbital_angular_momentum_y = planet_data['position_z'] * planet_data['velocity_x'] - planet_data['position_x'] * planet_data['velocity_z'];
        orbital_angular_momentum_z = planet_data['position_x'] * planet_data['velocity_y'] - planet_data['position_y'] * planet_data['velocity_x'];
        orbital_angular_momentum = np.sqrt(np.power(orbital_angular_momentum_x, 2) + np.power(orbital_angular_momentum_y, 2) + np.power(orbital_angular_momentum_z, 2));
        relative_orbital_angular_momentum_x = orbital_angular_momentum_x/orbital_angular_momentum
        relative_orbital_angular_momentum_y = orbital_angular_momentum_y/orbital_angular_momentum
        relative_orbital_angular_momentum_z = orbital_angular_momentum_z/orbital_angular_momentum
        numerator = relative_orbital_angular_momentum_x * planet_data['spin_x'] + \
                            relative_orbital_angular_momentum_y * planet_data['spin_y'] + \
                            relative_orbital_angular_momentum_z * planet_data['spin_z']
        denominator= np.sqrt(np.power(relative_orbital_angular_momentum_x, 2) + \
                                np.power(relative_orbital_angular_momentum_y, 2) + \
                                np.power(relative_orbital_angular_momentum_z, 2)) * \
                                np.sqrt(np.power(planet_data['spin_x'], 2) + np.power(planet_data['spin_y'], 2) + np.power(planet_data['spin_z'], 2))
        planet_obliquity = numerator / denominator
        ofilter = planet_obliquity <= 1.
        planet_obliquity[ofilter] = np.arccos(planet_obliquity[ofilter])*180./np.pi
        planet_obliquity[np.logical_not(ofilter)] = 1.e-6
        planet_computed_data['planet_obliquity'] = planet_obliquity

        ## Planet precession angle
        # https://en.wikipedia.org/wiki/Apsidal_precession
        # Angle defined between an arbitrary direction and the semi-major axis direction
        numerator = 1./np.sin(planet_obliquity) \
                    *(planet_data['position_x']*planet_data['spin_x'] \
                        + planet_data['position_y']*planet_data['spin_y'] \
                        + planet_data['position_z']*planet_data['spin_z'])

        denominator = np.sqrt(np.power(planet_data['position_x'], 2) + np.power(planet_data['position_y'], 2) + np.power(planet_data['position_z'], 2)) \
                        * planet_norm_spin
        planet_precession_angle = numerator / denominator
        ofilter = planet_precession_angle <= 1.
        planet_precession_angle[ofilter] = np.arccos(planet_precession_angle[ofilter])*180./np.pi
        planet_precession_angle[np.logical_not(ofilter)] = 1.e-6
        planet_computed_data['planet_precession_angle'] = planet_precession_angle

        planet_mass = planet_data['mass'][0]
        norm_spin = np.sqrt(np.power(star_data['spin_x'], 2) + np.power(star_data['spin_y'], 2) + np.power(star_data['spin_z'], 2))
        corotation_radius = ((posidonius.constants.G_SI*posidonius.constants.M_SUN*(star_mass+planet_mass))**(1/3.)) * ((norm_spin/posidonius.constants.DAY)**(-2./3.))/posidonius.constants.AU
        planet_computed_data['corotation_radius'] = corotation_radius

        e = planet_data['eccentricity']
        non_zero = e != 0
        alpha = np.ones(len(e))
        alpha[non_zero] = (1.+15./2.*e[non_zero]**2+45./8.*e[non_zero]**4+5./16.*e[non_zero]**6)*1./(1.+3.*e[non_zero]**2+3./8.*e[non_zero]**4)*1./(1.-e[non_zero]**2)**1.5
        pseudo_rot = alpha * np.sqrt(posidonius.constants.G_SI*posidonius.constants.M_SUN*(star_mass+planet_mass))
        #pseudo_synchronization_period  = 2.*np.pi / (pseudo_rot * (planet_data['semi-major_axis']*AU)**(-3./2.) * posidonius.constants.HOUR) # Hours
        pseudo_synchronization_period  = 2.*np.pi / (pseudo_rot * (planet_data['semi-major_axis']*posidonius.constants.AU)**(-3./2.) * posidonius.constants.DAY) # Days
        pseudo_synchronization_period[non_zero] = np.nan
        planet_computed_data['pseudo_synchronization_period'] = pseudo_synchronization_period
        Pseudosynchro_ratio = (1. +(15./2.)*np.power(planet_data['eccentricity'],2.) +(45./8.)*np.power(planet_data['eccentricity'],4.) +(5./16.)*np.power(planet_data['eccentricity'],6.) ) \
            / ( ( 1. +3.*np.power(planet_data['eccentricity'],2.) +(3./8.)*np.power(planet_data['eccentricity'],4.) ) *np.power(1. -np.power(planet_data['eccentricity'],2.),3./2.) )

        if universe_integrator_json['universe']['consider_effects']['tides']:
            ### Calculation of energydot and tidal flux, in W/m2
            # Gravitationl energy lost of the system due to dissipation
            # Masses in kg
            k2pdelta = 2.465278e-3 # Terrestrial planets (no gas)
            gravitational_energy_lost = posidonius.analysis.computation.energydot(planet_data['semi-major_axis']*posidonius.constants.AU, \
                                                    planet_data['eccentricity'], \
                                                    planet_norm_spin / posidonius.constants.DAY, \
                                                    planet_obliquity * np.pi/180.0, \
                                                    posidonius.constants.G_SI, \
                                                    planet_data['mass'] * posidonius.constants.M_SUN, \
                                                    star_mass * posidonius.constants.M_SUN, \
                                                    planet_data['radius'] * posidonius.constants.AU, \
                                                    k2pdelta * posidonius.constants.DAY)
            planet_computed_data['gravitational_energy_lost'] = gravitational_energy_lost

            # The tidal heat flux depends on the eccentricity and on the obliquity of the planet.
            # If the planet has no obliquity, no eccentricity and if its rotation is synchronized, the tidal heat flux is zero.
            mean_tidal_flux = gravitational_energy_lost / (4 * np.pi * np.power(planet_data['radius'] * posidonius.constants.AU, 2))
            tidal_effect = universe_integrator_json['universe']['particles'][int(key)]['tides']['effect']
            if 'OrbitingBody' in tidal_effect:
                dissipation_factor_scale = tidal_effect.get('ConstantTimeLag', {}).get('dissipation_factor_scale', 0.)
            elif 'CentralBody' in tidal_effect:
                dissipation_factor_scale = tidal_effect.get('ConstantTimeLag', {}).get('dissipation_factor_scale', 0.)
            else:
                dissipation_factor_scale = 0.
            mean_tidal_flux *= dissipation_factor_scale
            planet_computed_data['mean_tidal_flux'] = mean_tidal_flux

            denergy_dt = planet_data['denergy_dt'] * 6.90125e37 # conversation from Msun.AU^2.day^-3 to W
            planet_computed_data['denergy_dt'] = denergy_dt

            inst_tidal_flux = denergy_dt / (4 * np.pi * np.power(planet_data['radius'] * posidonius.constants.AU, 2))
            planet_computed_data['inst_tidal_flux'] = inst_tidal_flux
        else:
            zeros = np.zeros(len(planet_computed_data['planet_rotation_period']))
            planet_computed_data['gravitational_energy_lost'] = zeros
            planet_computed_data['mean_tidal_flux'] = zeros
            planet_computed_data['denergy_dt'] = zeros
            planet_computed_data['inst_tidal_flux'] = zeros


        ################################################################################
        ## Star obliquity (with respect to a planet)
        # https://en.wikipedia.org/wiki/Axial_tilt
        # Angle between the spin axis of the star and the planet's orbital plane (in other words, it's the inclination of the planet)
        # or, equivalently, the angle between its equatorial plane and orbital plane
        # At an obliquity of zero, the two axes point in the same direction; i.e., the rotational axis is perpendicular to the orbital plane.
        numerator = relative_orbital_angular_momentum_x * star_data['spin_x'] + \
                            relative_orbital_angular_momentum_y * star_data['spin_y'] + \
                            relative_orbital_angular_momentum_z * star_data['spin_z']
        denominator= np.sqrt(np.power(relative_orbital_angular_momentum_x, 2) + \
                                np.power(relative_orbital_angular_momentum_y, 2) + \
                                np.power(relative_orbital_angular_momentum_z, 2)) * \
                                np.sqrt(np.power(star_data['spin_x'], 2) + np.power(star_data['spin_y'], 2) + np.power(star_data['spin_z'], 2))
        star_obliquity = numerator / denominator
        ofilter = star_obliquity <= 1.
        star_obliquity[ofilter] = np.arccos(star_obliquity[ofilter])*180./np.pi
        star_obliquity[np.logical_not(ofilter)] = 1.e-6
        planet_computed_data['star_obliquity'] = star_obliquity
        ################################################################################
        # print(' >> >>', star_data['mass'])

        ################################################################################
        #G         =  6.6742367e-11            # m^3.kg^-1.s^-2
        #AU        =  1.49598e11               # m
        #Msun      =  1.98892e30               # kg
        planet_orbital_period = 2*np.pi*np.sqrt(np.power(planet_data['semi-major_axis']*posidonius.constants.AU, 3)/(posidonius.constants.G_SI*posidonius.constants.M_SUN*(star_data['mass']+planet_data['mass'])))
        planet_computed_data['orbital_period'] = planet_orbital_period/(posidonius.constants.DAY) # in days  
        planet_computed_data['mean_motion'] = 2.0*np.pi / planet_computed_data['orbital_period'] # in rad / days
        print(f"|\t Orbital period (t=0) : {planet_computed_data['orbital_period'][0]} days ")

        ################################################################################
        # 
        x = planet_data['position_x']
        y = planet_data['position_y']
        z = planet_data['position_z']
        r = np.sqrt( np.power(x,2) +np.power(y,2) +np.power(z,2) )
        vx = planet_data['velocity_x']
        vy = planet_data['velocity_y']
        vz = planet_data['velocity_z']
        eccentricity = planet_data['eccentricity']

        ################################################################################

        planet_computed_data['mean_longitude'] = planet_data['mean_anomaly'] + planet_data['longitude_of_perihelion']

        planet_computed_data['omega'] = np.sqrt( np.power(planet_data['spin_x'],2) +np.power(planet_data['spin_y'],2) +np.power(planet_data['spin_z'],2)  ) # in rad / days

        planet_computed_data['mean_motion'] = 2.0*np.pi / planet_computed_data['orbital_period'] # in rad / days


        ################################################################################

        if plot_22: # Plot precession periods
            print('|\t Compute Laplace coeff')

            k2 = planet_data['love_number']
            planet_mass = planet_data['mass'][0]
            planet_radius = planet_data['radius'][0]
            rg_2 = np.power(planet_data['radius'],2)*planet_data['radius_of_gyration_2']
            C = rg_2 * planet_data['mass'] / np.power(planet_data['radius'],2) / planet_data['mass']
            precess_cste = 0.5*(star_mass/planet_mass) *np.power((planet_radius/planet_data['semi-major_axis']),3) *planet_computed_data['omega']*k2 / C
            print(f'|\t For precess cste = {precess_cste} | check k2 {k2[0]} and C {C[0]} ')
            planet_computed_data['spin_precession_cste'] = precess_cste # in rad/days (from omega)
            planet_computed_data['spin_precession_period'] = 2.0*np.pi / (precess_cste*np.cos(planet_computed_data['planet_obliquity']*np.pi/180.)) # in years

            # time = planet_data['current_time']
            # long_asc_node  = planet_data['longitude_of_ascending_node']
            # g = diff(long_asc_node) / diff(time)
            # print(f"|\t len g {len(g)} len time {len(time)} ")
            planet_computed_data['node_precession_period'] = 0.#2.0*np.pi / g

        if plot_31 or plot_31b or plot_31all or plot_32 or plot_33 or history_output: # Check orbital angles
            # Check Synodic period
            print('> Compute synodic frequency')
            # planet_computed_data['omega'] = np.sqrt( np.power(planet_data['spin_x'],2) +np.power(planet_data['spin_y'],2) +np.power(planet_data['spin_z'],2)  ) # in rad / days
            print('  check spin ', planet_data['spin_z'][0], '(=', planet_computed_data['omega'][0], ')',  2.*np.pi / planet_data['spin_z'][0], ' days | ', ) # spin is rad / days
            # planet_computed_data['mean_motion'] = 2.0*np.pi / planet_computed_data['orbital_period'] # in rad / days
            print('  check orbital period ', planet_computed_data['orbital_period'][0],  ' days ' ) # in days
            print('  check mean motion ', planet_computed_data['mean_motion'][0], 2.*np.pi / planet_computed_data['mean_motion'][0], ' days '  ) # mean_motion is rad / days
            # 
            planet_computed_data['synodic_frequency'] = ( ( planet_computed_data['omega'] / planet_computed_data['mean_motion'] ) -1.0) *planet_computed_data['mean_motion'] # in rad / days
            planet_computed_data['synodic_period'] = (2.*np.pi / planet_computed_data['synodic_frequency'])/365.25  # in days 
            planet_computed_data['synodic_frequency'] = mean( ( planet_computed_data['omega'] / planet_computed_data['mean_motion'] ) -1.0) *planet_computed_data['mean_motion'] # in rad / days
            planet_computed_data['synodic_period_mean'] = (2.*np.pi / planet_computed_data['synodic_frequency'])/365.25  # in days 

            time = planet_data['current_time']
            time_day = [ t*365.25 for t in time]
            sub_stellar_integrated = cumulative_trapezoid(planet_computed_data['synodic_frequency'], x=time_day, initial=0.)
            planet_computed_data['synodic_rotation'] = sub_stellar_integrated # in rad

            # Compute Sub-stellar points
            print('|\t > Compute Sub-stellar point')
        
            # - Stellar position in the sky
            anglet = np.arctan2(y,x)*180./np.pi
            tmp = (anglet)+np.abs(anglet[0])
            print("   Solr Init :",anglet[0],(anglet[0])*180./np.pi, tmp[0]/180.*np.pi, tmp[0] )
            angle_cumul = []
            # angle_cumul.append(tmp[0])
            angle_cumul.append(0.)
            for dd in range(1,len(tmp)):
                if tmp[dd] > tmp[dd-1]:
                    bef = angle_cumul[len(angle_cumul)-1]
                    delta = bef + (tmp[dd] - tmp[dd-1])
                    angle_cumul.append(delta)
                else:
                    d1 = 360.- tmp[dd-1]
                    d2 = tmp[dd]
                    delta = d1 +d2
                    bef = angle_cumul[len(angle_cumul)-1]
                    angle_cumul.append(delta +bef)
            # ---
            # test _spin angle integrated 
            spin = planet_computed_data['omega']
            time_day = [ t*365.25 for t in planet_data['current_time']]
            theta = cumulative_trapezoid(spin, x=time_day, initial=0.)
            theta = [t*180./np.pi for t in theta]
            planet_computed_data['angle_cumul'] = angle_cumul
            diff = [ a-b for (a,b) in zip(angle_cumul,theta)] # diff = [ a-b for (a,b) in zip(angle_cumul,spin_angle)]
            init = diff[0]
            diff = [d - init for d in diff]
            planet_computed_data['substellar_point'] = diff

        if plot_33 or plot_33_1p or plot_32 or history_output: # Check orbital angles | Compute angle sub-stellar
            # --- compute cartesian position angle
            anglet = np.arctan2(y,x)*180./np.pi
            tmp = (anglet)+np.abs(anglet[0])
            print(f"|\t Cartesian position angle : {anglet[0]} deg |\t Computed init: {tmp[0]/180.*np.pi} rad  {tmp[0]} deg " )
            cartesian_position_angle = []
            # angle_cumul.append(tmp[0])
            cartesian_position_angle.append(0.)
            for dd in range(1,len(tmp)):
                if tmp[dd] > tmp[dd-1]:
                    bef = cartesian_position_angle[len(cartesian_position_angle)-1]
                    delta = bef + (tmp[dd] - tmp[dd-1])
                    cartesian_position_angle.append(delta)
                else:
                    d1 = 360.- tmp[dd-1]
                    d2 = tmp[dd]
                    delta = d1 +d2
                    bef = cartesian_position_angle[len(cartesian_position_angle)-1]
                    cartesian_position_angle.append(delta +bef)
            # ---
            planet_computed_data['cartesian_position_angle'] = cartesian_position_angle

            true_anomaly = planet_data['true_anomaly']
            # time = [ t*365.25 for t in planet_data['current_time']]
            spin = planet_computed_data['omega']
            mean_motion = planet_computed_data['mean_motion']
            # Compute true ano cumul
            true0 = true_anomaly[0]
            count = 0.
            table = []
            for true1 in true_anomaly:
                if (np.abs(true1-true0)>np.pi) :
                    count +=1.
                true0 = true1
                table.append(count)
            true_ano_cumul = [(tr +ta*2.*np.pi) for (tr,ta) in zip(true_anomaly, table)]
            planet_computed_data['true_anomaly_cumul'] = true_ano_cumul
            # Compute theta | Integrate Spin rate
            theta = cumulative_trapezoid(spin, x=time_day, initial=0.)
            planet_computed_data['theta'] = theta
            # Compute mean ano cumul
            mean_ano_cumul = cumulative_trapezoid(mean_motion, x=time_day, initial=0.)
            planet_computed_data['mean_ano_cumul'] = mean_ano_cumul
            # Compute true ano - theta
            theta_true_ano = [o-f for o,f in zip(theta, true_ano_cumul)]
            planet_computed_data['theta_true_ano'] = [t-theta_true_ano[0] for t in theta_true_ano]
            # Compute theta true ano (long asc node correction) - theta
            long_peri = planet_data['longitude_of_perihelion']
            long_peri = [(l+4.*np.pi*np.pi)%(2.0*np.pi) for l in long_peri]
            l0 = long_peri[0]
            long_peri_0 = [l-l0 for l in long_peri]
            planet_computed_data['longitude_of_perihelion_rescaled'] = long_peri_0
            true_ano_longperi_corrected = [t - l for t,l in zip(true_ano_cumul,long_peri_0)]
            planet_computed_data['true_ano_longperi_corrected'] = true_ano_longperi_corrected
            theta_true_ano_long_peri = [ f-l for f,l in zip(planet_computed_data['theta_true_ano'], long_peri_0)]
            planet_computed_data['theta_true_ano_long_peri'] = theta_true_ano_long_peri

            theta_mean_ano = [o-n for o,n in zip(theta, mean_ano_cumul)]
            planet_computed_data['theta_mean_ano'] = [t-theta_mean_ano[0] for t in theta_mean_ano]
            theta_cartesian_position_angle = [o-(n*np.pi/180.) for o,n in zip(theta, cartesian_position_angle)]
            planet_computed_data['theta_cartesian_position_angle'] = [t-theta_cartesian_position_angle[0] for t in theta_cartesian_position_angle]

        if plot_33_meanano :
            # semi_major_axis = planet_data['semi-major_axis']
            # eccentricity = planet_data['eccentricity']
            # (new_x, new_y, new_z) = (x, y, z)

            # long_asc_node = planet_data['longitude_of_ascending_node']
            # long_peri = planet_data['longitude_of_perihelion']
            # arg_peri = long_peri - long_asc_node

            gm = posidonius.constants.G*( star_mass + planet_data['mass'][0] )
            ( ex , ey, ez ) = calculate_eccentricity_vector(gm, x, y, z, vx, vy, vz)
            check_long_peri = (np.arctan2(ey,ex))%(2.0*np.pi)
            # check_long_peri = (-np.arctan2(ey,ex)+2.0*np.pi )%(2.0*np.pi)
            planet_computed_data['check_longitude_of_perihelion'] = (check_long_peri)#+np.pi)%(2.0*np.pi)
            # (new_x, new_y, new_z) = ( ex , ey, ez )
            # planet_computed_data['semi-minor_axis'] = np.sqrt( np.power(planet_data['semi-major_axis'],2) *(1.0 -np.power(planet_data['eccentricity'],2) ) )

            # nom = ey * planet_data['semi-major_axis']
            # denom = ex * planet_computed_data['semi-minor_axis']
            # planet_computed_data['eccentric_anomaly'] = np.arctan2(nom,denom)
            rv = x * vx  +  y * vy  +  z * vz
            v2 = np.power(vx,2) +np.power(vy,2) +np.power(vz,2)
            nom = r*v2 - gm
            denom = eccentricity*gm
            ce = nom/denom
            change_sign = np.logical_and(eccentricity, np.abs(ce) > 1.)
            if np.any(change_sign):
                ce[change_sign] = np.sign(ce[change_sign])

            bige = np.arccos(ce)
            negative_rv = rv < 0.
            if np.any(negative_rv):
                bige[negative_rv] = 2.0*np.pi - bige[negative_rv]
                # check_long_peri[negative_rv] = (2.0*np.pi - check_long_peri[negative_rv] )%(2.0*np.pi)
            bige = (bige)%(2.0*np.pi)
            planet_computed_data['eccentric_anomaly'] = bige
            ff = planet_computed_data['eccentric_anomaly']
            m = planet_data['mass'][0]
            print(f'|\t test {ff[0]} mass star {star_mass} planet{m} position {r[0]} cos {ce[0]}')

            planet_computed_data['check_mean_anomaly'] = (planet_computed_data['eccentric_anomaly'] - planet_data['eccentricity']*np.sin(planet_computed_data['eccentric_anomaly']) )
            # sin_bige = ey / planet_computed_data['semi-minor_axis']
            # planet_computed_data['check_mean_anomaly'] = (planet_computed_data['eccentric_anomaly'] - planet_data['eccentricity']*sin_bige )

            # planet_computed_data['mean_longitude'] = (planet_computed_data['check_mean_anomaly'] + planet_data['longitude_of_perihelion'])
            planet_computed_data['mean_longitude'] = (planet_computed_data['check_mean_anomaly'] + planet_computed_data['check_longitude_of_perihelion'])
        ################################################################################

        ################################################################################
        if plot_11 :
            print('> Compute angle modes ')
            print('   is any eccentricity > 0.3 :',  np.any(eccentricity > 0.3))

            if np.any(eccentricity > 0.3):
                order = 7
            elif np.any(eccentricity > 0.25):
                order = 6
            elif np.any(eccentricity > 0.2):
                order = 5
            elif np.any(eccentricity > 0.15):
                order = 4
            elif np.any(eccentricity > 0.1):
                order = 3
            elif np.any(eccentricity <= 0.1):
                order = 2
            else:
                order = 0

            planet_computed_data['order'] = order
            
            q_min = -(7-order)
            q_max = 2*order+2
            print('> order : ', order, ' | btwn : ', q_min, q_max)

            spin = planet_data['spin_z']
            mean_motion = planet_computed_data['mean_motion']

            frequ_mode_201q = []
            frequ_mode_220q = []
            period_mode_201q = []
            period_mode_220q = []
            for q in np.arange(q_min, q_max):
                print('   Compute order : ', q)
                q_201 = calculate_tidal_excitation_frequency_mode_sigma_2mpq(0, 1, q, spin, mean_motion)
                q_220 = calculate_tidal_excitation_frequency_mode_sigma_2mpq(2, 0, q, spin, mean_motion)
                frequ_mode_201q.append(q_201)
                frequ_mode_220q.append(q_220)
                period_mode_201q.append( 2.0*np.pi / q_201 )
                period_mode_220q.append( 2.0*np.pi / q_220 )
            print('   Frequ computed len : ', len(frequ_mode_220q), len(period_mode_220q) )
            planet_computed_data['frequ_mode_201q'] = frequ_mode_201q
            planet_computed_data['frequ_mode_220q'] = frequ_mode_220q
            planet_computed_data['period_mode_201q'] = period_mode_201q
            planet_computed_data['period_mode_220q'] = period_mode_220q

        if plot_12: # Print Hansen coeff 
            print('\n--- --- --- Frequency analysis --- --- --- ')
            print('\n --- --- Hansen Coef --- --- ')
            print('|\t > Compute order ')
            print('|\t   is any ecc > 0.3 :',  np.any(eccentricity > 0.3))
            order = 0
            if np.any(eccentricity > 0.3):
                order = 7
            elif np.any(eccentricity > 0.25):
                order = 6
            elif np.any(eccentricity > 0.2):
                order = 5
            elif np.any(eccentricity > 0.15):
                order = 4
            elif np.any(eccentricity > 0.1):
                order = 3
            elif np.any(eccentricity <= 0.1):
                order = 2
            else:
                order = 0
            # order = 2
            planet_computed_data['order'] = order
            q_min = -order # -(7-order)
            q_max = order #2*order+2
            print(f'|\t > order : {order} | q from {q_min} to {q_max}  || {(7-order)} {(7+order)} ')
            # for ordr in np.arange(0, 7):
                # print(f'|\t > idx : {(7-ordr)} {(7+ordr)} ')

            # Spin and mean_motion are in rad/days
            spin = planet_computed_data['omega']
            mean_motion = planet_computed_data['mean_motion']

            g_220q = calculate_eccentricity_function_g_20q(eccentricity)
            g_201q = calculate_eccentricity_function_g_21q(eccentricity)
            amp_mode_201q = []
            amp_mode_220q = []
            for q in np.arange(7-order, 7+order+1):
                print(f"|\t take {q}")
                amp_201q = g_201q[q]
                amp_220q = g_220q[q]
                amp_mode_201q.append((3./4.)*median(amp_201q)*median(amp_201q))
                amp_mode_220q.append((3./2.)*median(amp_220q)*median(amp_220q)+(9./4.)*median(amp_220q)*median(amp_220q))

            frequ_mode_201q = []
            frequ_mode_220q = []
            for q in np.arange(q_min, q_max+1):
                print(f'|\t  > - compute order {q} ')
                w_201q = calculate_tidal_excitation_frequency_mode_sigma_2mpq(0, 1, q, spin, mean_motion)
                w_220q = calculate_tidal_excitation_frequency_mode_sigma_2mpq(2, 0, q, spin, mean_motion)
                frequ_mode_201q.append(median(w_201q))
                frequ_mode_220q.append(median(w_220q))
                
            print(f'|\t > Control frequ len | 201q : {len(frequ_mode_201q)} | 220q : {len(frequ_mode_220q)} | ')
            planet_computed_data['frequ_mode_201q'] = frequ_mode_201q
            planet_computed_data['frequ_mode_220q'] = frequ_mode_220q
            planet_computed_data['amp_mode_201q'] = amp_mode_201q
            planet_computed_data['amp_mode_220q'] = amp_mode_220q

            print('\n --- --- NAF Coef --- --- ')
            true_anomaly = planet_data['true_anomaly']
            mean_anomaly = planet_data['mean_anomaly']
            time = planet_data['current_time']
            # theta_true_ano = planet_computed_data['theta_true_ano']
            # theta_mean_ano = planet_computed_data['theta_mean_ano']
            (naff_frequs, naff_amps) = posidonius.tools.naf(time, true_anomaly, nfreqs = 29)
            # (naff_frequs, naff_amps) = posidonius.tools.naf(time, theta_true_ano, nfreqs = 18)
            # naff_periods = [2.0*np.pi / f *365.25 for f in naff_frequs]
            naff_frequs = [f /365.25 /86400. for f in naff_frequs]
            naff_periods = [2.0*np.pi / f /86400. for f in naff_frequs]
            print(f"|\t frequ {naff_frequs} \n|\t periods {naff_periods} \n|\t amp {naff_amps} \n|\t {np.sqrt(np.power(naff_amps.real,2) +np.power(naff_amps.imag,2))} ")
            planet_computed_data['naff_frequs'] = naff_frequs
            planet_computed_data['naff_periods'] = naff_periods
            planet_computed_data['naff_amps'] = np.sqrt(np.power(naff_amps.real,2) +np.power(naff_amps.imag,2))

        ################################################################################
        ## Sum over all the planets values:
        planet_angular_momentum_x = planet_data['radius_of_gyration_2'][0] * (planet_data['mass'][0]) * np.power(planet_data['radius'][0], 2) * planet_data['spin_x']
        planet_angular_momentum_y = planet_data['radius_of_gyration_2'][0] * (planet_data['mass'][0]) * np.power(planet_data['radius'][0], 2) * planet_data['spin_y']
        planet_angular_momentum_z = planet_data['radius_of_gyration_2'][0] * (planet_data['mass'][0]) * np.power(planet_data['radius'][0], 2) * planet_data['spin_z']
        planet_angular_momentum = planet_data['radius_of_gyration_2'][0] * (planet_data['mass'][0]) * np.power(planet_data['radius'][0], 2) * (planet_norm_spin)
        #planet_angular_momentum = planet_data['radius_of_gyration_2'] * (planet_data['mass']) * np.power(planet_data['radius'], 2) * (planet_norm_spin)
        total_planets_angular_momentum_x += planet_angular_momentum_x # If more than one planet is present, all of them should be added
        total_planets_angular_momentum_y += planet_angular_momentum_y # If more than one planet is present, all of them should be added
        total_planets_angular_momentum_z += planet_angular_momentum_z # If more than one planet is present, all of them should be added
        #total_planets_angular_momentum += planet_angular_momentum # If more than one planet is present, all of them should be added
        # If more than one planet is present, all of them should be added
        total_planets_angular_momentum = np.sqrt(np.power(total_planets_angular_momentum_x, 2) + np.power(total_planets_angular_momentum_y, 2) + np.power(total_planets_angular_momentum_z, 2))
        ################################################################################

        ################################################################################
        # Save computed data
        planets_computed_data[key] = planet_computed_data
        ################################################################################

    ################################################################################
    ### [start] compute total planets orbital angular momentum for one star and N planets
    factors_a = []
    cross_products = []
    r_terms = []
    v_terms = []
    accumulative_mass = star_mass
    for key in planets_keys:
        planet_mass_msun = planets_data[key]['mass'][0]
        numerator = accumulative_mass * planet_mass_msun
        accumulative_mass += planet_mass_msun
        factor_a = numerator / accumulative_mass
        factors_a.append(factor_a)

        r = np.vstack((planets_data[key]['position_x'],planets_data[key]['position_y'], planets_data[key]['position_z']))
        v = np.vstack((planets_data[key]['velocity_x'], planets_data[key]['velocity_y'], planets_data[key]['velocity_z']))
        for (r_term, v_term) in zip(r_terms, v_terms):
            r -= r_term
            v -= v_term
        r_terms.append(planet_mass_msun/accumulative_mass * r)
        v_terms.append(planet_mass_msun/accumulative_mass * v)

        x = np.cross(r.T, v.T)
        cross_products.append(x)

    lorb = np.zeros((len(star_data), 3))
    for (a, rxv) in zip(factors_a, cross_products):
        lorb += a*rxv
    total_planets_orbital_angular_momentum_x = lorb[:,0]  # Msun.AU^2.DAY-1
    total_planets_orbital_angular_momentum_y = lorb[:,1]  # Msun.AU^2.DAY-1
    total_planets_orbital_angular_momentum_z = lorb[:,2]  # Msun.AU^2.DAY-1
    total_planets_orbital_angular_momentum = np.sqrt(np.sum(np.power(lorb, 2), axis=1))  # Msun.AU^2.DAY-1
    ### [end] compute total planets orbital angular momentum for one star and N planets
    ################################################################################


    # \Delta L / L
    star_angular_momentum_x = star_data['radius_of_gyration_2'] * (star_mass) * np.power(star_data['radius'], 2) * star_data['spin_x']
    star_angular_momentum_y = star_data['radius_of_gyration_2'] * (star_mass) * np.power(star_data['radius'], 2) * star_data['spin_y']
    star_angular_momentum_z = star_data['radius_of_gyration_2'] * (star_mass) * np.power(star_data['radius'], 2) * star_data['spin_z']
    star_angular_momentum = star_data['radius_of_gyration_2'] * (star_mass) * np.power(star_data['radius'], 2) * (star_norm_spin)

    initial_total_angular_momentum_x = total_planets_orbital_angular_momentum_x[0] + total_planets_angular_momentum_x[0] + star_angular_momentum_x[0]
    initial_total_angular_momentum_y = total_planets_orbital_angular_momentum_y[0] + total_planets_angular_momentum_y[0] + star_angular_momentum_y[0]
    initial_total_angular_momentum_z = total_planets_orbital_angular_momentum_z[0] + total_planets_angular_momentum_z[0] + star_angular_momentum_z[0]
    initial_total_angular_momentum = np.sqrt(np.power(initial_total_angular_momentum_x, 2) + np.power(initial_total_angular_momentum_y, 2) + np.power(initial_total_angular_momentum_z, 2))
    total_angular_momentum_x = total_planets_orbital_angular_momentum_x + total_planets_angular_momentum_x + star_angular_momentum_x
    total_angular_momentum_y = total_planets_orbital_angular_momentum_y + total_planets_angular_momentum_y + star_angular_momentum_y
    total_angular_momentum_z = total_planets_orbital_angular_momentum_z + total_planets_angular_momentum_z + star_angular_momentum_z
    total_angular_momentum = np.sqrt(np.power(total_angular_momentum_x, 2) + np.power(total_angular_momentum_y, 2) + np.power(total_angular_momentum_z, 2))
    conservation_of_angular_momentum = np.abs((total_angular_momentum - initial_total_angular_momentum)/initial_total_angular_momentum)

    if len(conservation_of_angular_momentum) > 1:
        conservation_of_angular_momentum[0] = conservation_of_angular_momentum[1]

    print("\n< Let's plot > ")

    if plot_1:# Plot AM conservation
        fig = plt.figure(figsize=(16, 10))
        ligne = 1
        colonne = 2
        i = 0

        i=i+1
        ax = fig.add_subplot(ligne,colonne,i,)
        field = 'L'
        time = [ t*365.25 for t in planets_data['1']['current_time']]
        ax.plot(time, total_angular_momentum, '+', label=key)
        ax.set_ylabel(field)
        # ax.set_xscale('log')
        # ax.set_yscale('log')
        ax.legend(loc=0, prop={'size':8})

        i=i+1
        ax = fig.add_subplot(ligne,colonne,i,)
        field = '$\Delta L/L_{0}$'
        ax.plot(time, conservation_of_angular_momentum, label=key)
        ax.set_ylabel(field)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(loc=0, prop={'size':8})

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_AM_conservation" +".png")
        
        plt.savefig(output_figure_filename)
        # plt.show()
        plt.close()
        print("> Output figure file written to: {}".format(output_figure_filename))
    if plot_1:# Plot time step
        fig = plt.figure(figsize=(10, 8))
        ligne = 1
        colonne = 1
        i = 0
        #-----------------------------------------------------------------------------------    
        i=i+1
        ax = fig.add_subplot(ligne,colonne,i,)
        field = 'time_step'
        mu = []
        for key in planets_keys:
            planet_data = planets_data[key]
            time = planet_data['current_time']
            time_step = planet_data['time_step']
            ax.plot(time, time_step, label=key)

        ax.set_ylabel(field+' days ', fontsize=14)
        ax.set_xscale('log')
        ax.set_xlabel('years')
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        ax.legend(loc=0, prop={'size':8})

        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_time_step" +".png")
        plt.savefig(output_figure_filename)
        # plt.show()
        plt.close()
        print("> Output figure file written to: {}".format(output_figure_filename))
    
    if plot_2:# Plot data
        fig = plt.figure(figsize=(16, 22))
        ligne = 4
        colonne = 2
        i = 0
        date = 0. #5.96e3
        #-----------------------------------------------------------------------------------    
        i=i+1
        ax = fig.add_subplot(ligne,colonne,i,)
        field = 'semi-major_axis'
        mu = []
        for key in planets_keys:
            planet_data = planets_data[key]
            time = planet_data['current_time']
            sma = planet_data[field]
            ecc = planet_data['eccentricity']
            peri = sma *(1-ecc)
            apo = sma *(1+ecc)
            ls, = ax.plot(time, sma, label=field+"AU "+key)
            # ax.fill_between(time, apo, peri, color=ls.get_color(), alpha=0.1)
            planet_mass = planets_data[key]['mass'][0]
            mu.append(posidonius.constants.G_SI*( star_mass*posidonius.constants.M_SUN + planet_mass*posidonius.constants.M_SUN ))
        # planet_period = (1.510826, 2.421937, 4.049219, 6.101013, 9.207540, 12.352446, 18.772866)
        # tmp = [ np.power(t*86400. /(2.0*np.pi),2)*m for t,m in zip(planet_period,mu)]
        # planet_a = (1.154e-2, 1.580e-2, 2.227e-2, 2.925e-2, 3.849e-2, 4.683e-2, 6.189e-2)
        # planet_a = [np.power(n,1./3.) / posidonius.constants.AU for n in tmp]
        # for sma in planet_a:
        #     ax.axhline(y=sma, ls=':', color="black")
        ax.set_ylabel(field+' (AU)', fontsize=14)
        ax.set_xscale('log')
        ax.set_xlabel('years')
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        ax.legend(loc=0, prop={'size':8})

        # i=i+1
        a = plt.axes([0.12, 0.80, .18, .10],) # ( position[x,y], size[x,y] )
        field = ' $\Delta$ sma'
        for key in planets_keys:
            planet_data = planets_data[key]
            time = planet_data['current_time'][planet_data['current_time']>date]
            sma = planet_data['semi-major_axis'][planet_data['current_time']>date]
            sma0 = sma[0]
            d_sma  = [(sam - sma0) / sma0 for sam in sma]
            plt.plot(time, d_sma, label=" planet "+key)
        plt.ylabel(field, fontsize=14)
        plt.xscale('log')
        # ax.set_yscale('log')
        # plt.xlabel('years')
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        # plt.tight_layout()

        i=i+1
        ax = fig.add_subplot(ligne,colonne,i,)
        field = 'period ratio'
        time = planets_data['1']['current_time']
        period_tmp = planets_computed_data['1']['orbital_period']
        for key in planets_keys:
            if key !='1':
                planet_data = planets_data[key]
                planet_computed_data = planets_computed_data[key]
                period = planet_computed_data['orbital_period']
                ratio = period/ period_tmp
                ax.plot(time, ratio, label=" period ration $p_{i+1}$/$p_i$ ")
        ax.set_ylabel(field, fontsize=14)
        ax.set_xscale('log')
        # ax.set_yscale('log')
        ax.set_xlabel('years')
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        ax.legend(loc=0, prop={'size':12})
        ax.legend(loc=0, prop={'size':8})
        plt.tight_layout()

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'Planet rotation state Omega/n '
        for key in planets_keys:
            print('> Planet Mass', planet_data['mass'][0])
            planet_data = planets_data[key]
            print('> End Time::%.4f'%(planet_data['current_time'][len(planet_data['current_time'])-1]))
            omega_n = planets_computed_data[key]['orbital_period'] / planets_computed_data[key]['planet_rotation_period']
            omega = (2.*np.pi/ (planets_computed_data[key]['planet_rotation_period']*posidonius.constants.DAY))
            n = ( 2.*np.pi/ (planets_computed_data[key]['orbital_period']*posidonius.constants.DAY))
            frequ_on = 2.0*(omega - n)
            ax1.plot(planet_data['current_time'], omega_n, label=key +' Pinit= %10.3e'%(planets_computed_data[key]['planet_rotation_period'][0]*24.)+"hrs")
        ax1.set_ylabel(field, fontsize=14)
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':12})
        a = plt.axes([0.12, 0.52, .18, .10],) # ( position[x,y], size[x,y] )
        field = ' zoom omega_n '
        for key in planets_keys:
            planet_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]
            time = planet_data['current_time'][planet_data['current_time']>date]
            omega_n = planet_computed_data['orbital_period'] / planet_computed_data['planet_rotation_period']
            omega_n = omega_n [planet_data['current_time']>date]
            # omega_n0 = onega_n[0]
            # d_omega  = [(sam - sma0) / sma0 for sam in sma]
            plt.plot(time, omega_n, label=" planet "+key)
        plt.ylabel(field, fontsize=14)
        plt.xscale('log')
        plt.ylim(top=1.00001, bottom=0.99999)
        # ax.set_yscale('log')
        # plt.xlabel('years')
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        # plt.tight_layout()

        i=i+1
        ax = fig.add_subplot(ligne,colonne,i,)
        field = 'eccentricity'
        for key in planets_keys:
            planet_data = planets_data[key]
            # if (int(key) == 2):
            ax.plot(planet_data['current_time'], planet_data[field], label="Planet " +key)
        ax.set_ylabel(field, fontsize=14)
        ax.set_xscale('log')
        # ax.set_yscale('log')
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        ax.legend(loc=0, prop={'size':12})

        i=i+1
        ax = fig.add_subplot(ligne,colonne,i,)
        field = 'planet obliquity deg'
        for key in planets_keys:
            planet_data = planets_data[key]
            obliquity = planets_computed_data[key]['planet_obliquity']
            ax.plot(planet_data['current_time'], obliquity, label=field+" "+key +f" init : {round(obliquity[0],3)} deg| final : {round(obliquity[-1],3)} deg")
            # ax.plot(planet_data['current_time'], planet_keplerian_orbital_elements[0], label=field+" ReCAlAU")
        ax.set_ylabel(field, fontsize=14)
        ax.set_xscale('log')
        ax.set_xlabel('years')
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        ax.legend(loc=0, prop={'size':12})
        # i=i+1
        # ax = fig.add_subplot(ligne,colonne,i,)
        a = plt.axes([0.12, 0.32, .18, .10],) # ( position[x,y], size[x,y] )
        field = ' $\Delta$ Obliquity'
        for key in planets_keys:
            planet_data = planets_data[key]
            time = planet_data['current_time'][planet_data['current_time']>date]
            obliquity = planets_computed_data[key]['planet_obliquity'][planet_data['current_time']>date]
            ob0 = obliquity[0]
            d_obl  = [(ob - ob0) / ob0 for ob in obliquity]
            plt.plot(time, d_obl, label=" planet "+key)
        plt.ylabel(field, fontsize=14)
        plt.xscale('log')
        # ax.set_yscale('log')
        # plt.xlabel('years')
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        # plt.tight_layout()

        i=i+1
        ax = fig.add_subplot(ligne,colonne,i,)
        field = 'inclination'
        for key in planets_keys:
            planet_data = planets_data[key]
            inclination = [i*180./np.pi for i in planet_data[field]]
            ax.plot(planet_data['current_time'], inclination, label="Planet "+key+f" init : {round(planet_data[field][0]*180/np.pi,3)} deg | final : {round(planet_data[field][-1]*180./np.pi,3)} deg")
        ax.set_ylabel(field, fontsize=14)
        ax.set_xscale('log')
        # ax.set_yscale('log')
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        ax.legend(loc=0, prop={'size':12})

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'Precession angle'
        for key in planets_keys:
            planeplot_2t_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]
            time = planet_data['current_time']
            precession_angle = planet_computed_data['planet_precession_angle']
            ax1.plot(time, precession_angle)
        ax1.set_ylabel(field, fontsize=14)
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        # ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})
        
        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'long perihelion'
        aa=1.0
        for key in planets_keys:
            planet_data = planets_data[key]
            time = planet_data["current_time"]
            ecc = planet_data['eccentricity']
            longitude_peri = planet_data['longitude_of_perihelion'] * 180 / np.pi
            ax1.plot(time,longitude_peri, '+', alpha=aa)
            aa-=0.1
        ax1.set_ylabel('longitude peri', fontsize=14)
        # ax1.set_xscale('log')
        ax1.set_xlabel('Time years')
        # ax1.legend(loc=0, prop={'size':8})

        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_orbit_data" +".png")
        plt.savefig(output_figure_filename)
        # plt.show()
        plt.close()
        print("> Output figure file written to: {}".format(output_figure_filename))
    if plot_2b:# Plot data
        fig = plt.figure(figsize=(25, 7))
        ligne = 1
        colonne = 3
        i = 0
        #-----------------------------------------------------------------------------------    
        i=i+1
        ax = fig.add_subplot(ligne,colonne,i,)
        field = 'semi-major_axis'
        mu = []
        for key in planets_keys:
            planet_data = planets_data[key]
            time = planet_data['current_time']
            sma = planet_data[field]
            ecc = planet_data['eccentricity']
            peri = sma *(1-ecc)
            apo = sma *(1+ecc)
            ls, = ax.plot(time, sma, label=f"Planet {key} sma final: {sma[-1]:1.3e} AU ")
            ax.fill_between(time, apo, peri, color=ls.get_color(), alpha=0.1)
            planet_mass = planets_data[key]['mass'][0]
            mu.append(posidonius.constants.G_SI*( star_mass*posidonius.constants.M_SUN + planet_mass*posidonius.constants.M_SUN ))
        ax.set_ylabel(field+' (AU)', fontsize=14)
        ax.set_xscale('log')
        ax.set_xlabel('years')
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        ax.legend(loc=0, prop={'size':10})

        i=i+1
        ax = fig.add_subplot(ligne,colonne,i,sharex=ax)
        field = 'period ratio'
        time = planets_data['1']['current_time']
        period_tmp = planets_computed_data['1']['orbital_period']
        for key in planets_keys:
            if key !='1':
                planet_data = planets_data[key]
                planet_computed_data = planets_computed_data[key]
                period = planet_computed_data['orbital_period']
                ratio = period/ period_tmp
                ax.plot(time, ratio, label=" period ration $p_{i+1}$/$p_i$ : "+ str(ratio[-1]))
        ax.set_ylabel(field, fontsize=14)
        ax.set_xscale('log')
        # ax.set_yscale('log')
        ax.set_xlabel('years')
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        ax.legend(loc=2, prop={'size':12})
        # ax.legend(loc=0, prop={'size':8})
        plt.tight_layout()

        i=i+1
        ax = fig.add_subplot(ligne,colonne,i,sharex=ax)
        field = 'eccentricity'
        for key in planets_keys:
            planet_data = planets_data[key]
            # if (int(key) == 2):
            ecc = planet_data['eccentricity']
            ax.plot(planet_data['current_time'], ecc, label=f"Planet {key} ecc final : {ecc[-1]:1.3e} " )
        ax.set_ylabel(field, fontsize=14)
        ax.set_xscale('log')
        # ax.set_yscale('log')
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        ax.legend(loc=0, prop={'size':12})

        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_orbit_data" +".png")
        plt.savefig(output_figure_filename)
        # plt.show()
        plt.close()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_20:# Plot data
        fig = plt.figure(figsize=(10, 6))
        ligne = 1
        colonne = 1
        i = 0
        date = 0. #5.96e3
        #-----------------------------------------------------------------------------------    

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'Planet rotation state $\Omega$/n '
        for key in planets_keys:
            print('> Planet Mass', planet_data['mass'][0])
            planet_data = planets_data[key]
            print('> End Time::%.4f'%(planet_data['current_time'][len(planet_data['current_time'])-1]))
            omega_n = planets_computed_data[key]['orbital_period'] / planets_computed_data[key]['planet_rotation_period']
            ax1.plot(planet_data['current_time'], omega_n)
        ax1.axhline(y=1.0, ls=":", color='black')
        ax1.set_ylabel(field, fontsize=14)
        ax1.set_xscale('log')
        ax1.set_xlabel('Time (years)', fontsize=14)
        # ax1.set_xlim( right=5e2)
        y_ticks = [1-1.5e-6, 1-1e-6, 1-0.5e-6, 1., 1+0.5e-6, 1+1e-6, 1+1.5e-6 ]  # Adjust the range as needed
        y_tick_labels = ["1- 1.5e-6", "1- 1e-6", "1- 0.5e-6", "1.0", "1+ 0.5e-6", "1+ 1e-6", "1+ 1.5e-6" ]
        ax1.legend(loc=0, prop={'size':12})
        # plt.yticks(y_ticks, y_tick_labels, fontsize=14)
        plt.xticks(fontsize=14)
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_omega_data" +".png")
        plt.savefig(output_figure_filename)
        # plt.show()
        plt.close()
        print("> Output figure file written to: {}".format(output_figure_filename))
      
    if plot_21:# Plot positions
        fig = plt.figure(figsize=(8, 8))
        ligne = 1
        colonne = 1
        i = 0
        #-----------------------------------------------------------------------------------

        i=i+1
        ax = fig.add_subplot(ligne,colonne,i,)
        field = 'Planets coordinates'
        for key in planets_keys:
            planet_data = planets_data[key]

            # time = planet_data['current_time']
            
            x = planet_data['position_x'][planet_data['current_time']< 1./365.35]
            y = planet_data['position_y'][planet_data['current_time']< 1./365.35]

            ls, = ax.plot(x, y, label=key)
            ax.plot(x[0], y[0], 'x', color=ls.get_color())
        ax.plot(0., 0., '+', color='black')
        ax.set_ylabel(field+' (AU)')
        ax.set_ylim(top=0.04, bottom=-0.04)
        ax.set_xlim(right=0.04, left=-0.04)
        # ax.set_xscale('log')
        ax.set_xlabel('(AU)')
        ax.legend(loc=0, prop={'size':8})

        #-----------------------------------------------------------------------------------
        ax.legend(loc=0, prop={'size':8})
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_orbit_positions_data" +".png")
        # plt.show()
        plt.savefig(output_figure_filename)
        plt.close()
        print("> Output figure file written to: {}".format(output_figure_filename))
      
    if plot_22:# Plot precession periods
        # Parameters for b_s^m (alpha_12)
        p=1
        q=1
        s=3./2.
        m=1
        alpha = planets_data['1']['semi-major_axis']/planets_data['2']['semi-major_axis']
        n1 = planets_computed_data['1']['mean_motion']
        n2 = planets_computed_data['2']['mean_motion']
        m1 = planets_data['1']['mass']
        m2 = planets_data['2']['mass']
        
        laplace = LaplaceCoefficient()
        b = [-0.25*laplace(a, s, m, p, q) for a in alpha]
        g_ll =  b *alpha *( alpha*n1*m2/(star_mass+m1) +n2*m1/(star_mass+m2) )
        node_precession_period = 2.0*np.pi /g_ll
        print(f'|\t Compare frequencies g_LL {g_ll[0]} and period {node_precession_period[0]} ')

        fig = plt.figure(figsize=(16, 13))
        ligne = 3
        colonne = 2
        i = 0

        i=i+1
        ax = fig.add_subplot(ligne,colonne,i,)
        field = 'semi-major_axis' 
        for key in planets_keys:
            planet_data = planets_data[key]
            time = planet_data['current_time']
            sma = planet_data[field]
            ecc = planet_data['eccentricity']
            peri = sma *(1-ecc)
            apo = sma *(1+ecc)
            ls, = ax.plot(time, sma, label=field+"AU "+key)
            ax.fill_between(time, apo, peri, color=ls.get_color(), alpha=0.1)
            planet_mass = planets_data[key]['mass'][0]
        ax.set_ylabel(field+' (AU)', fontsize=16)
        ax.set_xscale('log')
        ax.set_xlabel('years')
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16)
        ax.legend(loc=0, prop={'size':8})

        i=i+1
        ax = fig.add_subplot(ligne,colonne,i,)
        field = 'inclination'
        for key in planets_keys:
            planet_data = planets_data[key]
            inclination = [i*180./np.pi for i in planet_data[field]]
            ax.plot(planet_data['current_time'], inclination, label="Planet "+key+f" init : {round(planet_data[field][0]*180/np.pi,3)} deg | final : {round(planet_data[field][-1]*180./np.pi,3)} deg")
        ax.set_ylabel(field, fontsize=10)
        ax.set_xscale('log')
        # ax.set_yscale('log')
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        ax.legend(loc=0, prop={'size':12})
        
        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'Obliquity'
        for key in planets_keys:
            planet_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]
            time = planet_data['current_time']
            orbital_inclination = [p*180./np.pi for p in planet_data['inclination']]
            obliquity = planet_computed_data['planet_obliquity']
            ls, =ax1.plot(time, obliquity, label=f'Obliqu init: {round(obliquity[0],3)} | final: {round(obliquity[-1],3)}')
            # ax1.plot(time, orbital_inclination, ls=':', color=ls.get_color(),label=f'Inclin init: {round(orbital_inclination[0],3)} | final: {round(orbital_inclination[-1],3)}')
        ax1.set_ylabel(field, fontsize=16)
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16)
        ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8}) #

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'Precession angle'
        for key in planets_keys:
            planet_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]
            time = planet_data['current_time']
            precession_angle = planet_computed_data['planet_precession_angle']
            ax1.plot(time, precession_angle)
        ax1.set_ylabel(field, fontsize=16)
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16)
        # ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        i=i+1
        ax = fig.add_subplot(ligne,colonne,i,)
        field = 'longitude_of_ascending_node'
        for key in planets_keys:
            planet_data = planets_data[key]
            long_asc_node =  planet_data['longitude_of_ascending_node'] # [i*180./np.pi for i in planet_data[field]]
            ax.plot(planet_data['current_time'], long_asc_node, label="Planet "+key)
        ax.set_ylabel(field, fontsize=10)
        ax.set_xscale('log')
        # ax.set_yscale('log')
        plt.yticks(fontsize=14)
        plt.xticks(fontsize=14)
        ax.legend(loc=0, prop={'size':12})

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'Precession periods'
        for key in planets_keys:
            planet_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]
            time = planet_data['current_time']
            spin_precession_period = planet_computed_data['spin_precession_period']
            numerical_node_precession_period = planet_computed_data['node_precession_period']
            ls,=ax1.plot(time, spin_precession_period, label='spin precession period 2  pi / a cos e '+key)
            # ax1.plot(time[:-1], numerical_node_precession_period,":", label='numerical node precession period 2 pi / g '+key, color=ls.get_color(), )
        ax1.plot(time, node_precession_period,color='black', label=' node precession period 2 pi / g')
        ax1.set_ylabel(field, fontsize=16)
        # ax1.set_ylim(top=1000, bottom=0)
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16)
        ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})


        #-----------------------------------------------------------------------------------
        ax.legend(loc=0, prop={'size':8})
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_precession_data" +".png")
        # plt.show()
        plt.savefig(output_figure_filename)
        plt.close()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_23:# Print inital conditions 
        lab={}
        lab[1]='T1-b'
        lab[2]='T1-c'
        lab[3]='T1-d'
        lab[4]='T1-e'
        lab[5]='T1-f'
        lab[6]='T1-g'
        lab[7]='T1-h'
        fig = plt.figure(figsize=(10, 10))
        ligne = 1
        colonne = 1
        #-----------------------------------------------------------------------------------

        fig, ax = plt.subplots(ligne, colonne, figsize=(15, 15))

        table_times = {}
        table_times['planet'] = []
        # table_times['time'] = []
        table_times['sma'] = []
        table_times['ecc'] = []
        table_times['incl'] = []
        table_times['long_peri'] = []
        table_times['mean_ano'] = []
        table_times['long_asc_node'] = []
        table_times['mean_long'] = []

        for key in planets_keys:
            planet_data = planets_data[key]
            # time = planet_data['current_time']
            sma = planet_data['semi-major_axis']
            ecc = planet_data['eccentricity']
            incl = planet_data['inclination']
            long_peri = planet_data['longitude_of_perihelion']
            long_asc_node = planet_data['longitude_of_ascending_node']
            mean_ano = planet_data['mean_anomaly']
            true_ano = planet_data['true_anomaly']
            mean_long = long_peri + mean_ano 

            table_times['planet'].append(key)
            # table_times['time'].append(time[0])
            table_times['sma'].append(sma[0])
            table_times['ecc'].append(ecc[0])
            table_times['incl'].append(incl[0])
            table_times['long_peri'].append(long_peri[0])
            table_times['mean_ano'].append(mean_ano[0])
            table_times['long_asc_node'].append(long_asc_node[0])
            table_times['mean_long'].append(mean_long[0])

        columns = ( 'sma', 'ecc', 'incl', 'long_peri', 'mean_ano', 'long_asc_node', 'mean_long' )
        rows = [ lab[int(in_table)] for in_table in table_times['planet'] ]
        table = np.transpose((table_times['sma'], table_times['ecc'], table_times['incl'], table_times['long_peri'], table_times['mean_ano'], table_times['long_asc_node'], table_times['mean_long']))

        the_table = ax.table(cellText=table, loc='center', colLabels=columns, rowLabels=rows)
        ax.axis('tight')
        ax.axis('off')

        # columns = ('Time', 'angle')
        # rows = [ lab[int(in_table)] for in_table in table_times['planet'] ]
        # table = np.transpose((table_times['time'], table_times['sub']))
        # the_table = ax[1].table(cellText=table, loc='center', colLabels=columns, rowLabels=rows)
        # ax[1].axis('tight')
        # ax[1].axis('off')

        
        #-----------------------------------------------------------------------------------
        # plt.tight_layout()
        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_table_data" +".png")
        plt.savefig(output_figure_filename, dpi=200)
        # plt.show()
        plt.close()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_3e:# plot ecc data
        lab={}
        lab[1]='T1-b'
        lab[2]='T1-c'
        lab[3]='T1-d'
        lab[4]='T1-e'
        lab[5]='T1-f'
        lab[6]='T1-g'
        lab[7]='T1-h'
        fig = plt.figure(figsize=(18, 10))
        ligne = 2
        colonne = 2
        i = 0
        #-----------------------------------------------------------------------------------
        
        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'ecc'
        aa=1.0
        for key in planets_keys:
            planet_data = planets_data[key]
            #tmp = 9.0
            time = planet_data['current_time']
            ecc = planet_data['eccentricity']
            ax1.plot(time, ecc,alpha=aa)
            aa-=0.1
        ax1.set_ylabel(field)
        ax1.set_xscale('log')
        # ax1.set_yscale('log')
        ax1.set_xlabel('years')
        # ax1.legend(loc=0, prop={'size':8})

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'long perihelion'
        aa=1.0
        for key in planets_keys:
            planet_data = planets_data[key]
            time = planet_data["current_time"]
            ecc = planet_data['eccentricity']
            longitude_peri = planet_data['longitude_of_perihelion']
            # longitude_peri = planets_computed_data[key]['longitude_of_perihelion']
            # ecos = ecc*np.cos(longitude_peri)
            # esin = ecc*np.sin(longitude_peri)
            # long_peri = np.arctan2(ey,ex)
            # ax1.plot(ecos,esin)
            ax1.plot(time,longitude_peri, '+', alpha=aa)
            aa-=0.1
        ax1.set_ylabel('longitude perihelion')
        ax1.set_xscale('log')
        ax1.set_xlabel('Time years')
        # ax1.legend(loc=0, prop={'size':8})

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'e sin $\omega$'
        aa=1.0
        for key in planets_keys:
            planet_data = planets_data[key]
            #tmp = 9.0
            time = planet_data['current_time']
            ecc = planet_data['eccentricity']
            longitude_peri = planet_data['longitude_of_perihelion']
            # longitude_peri = planets_computed_data[key]['check_longitude_of_perihelion']
            ecos = ecc*np.cos(longitude_peri)
            esin = ecc*np.sin(longitude_peri)
            ax1.plot(ecos, esin, alpha=aa)
            aa-=0.12
        ax1.set_ylabel(field)
        # ax1.set_xscale('log')
        # ax1.set_yscale('log')
        ax1.set_xlabel(' e cos $\omega$')
        # ax1.legend(loc=0, prop={'size':8})
        # plt.yticks(fontsize=12)
        # plt.xticks(fontsize=12)

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = '$\omega_c$ - $\omega_b$'
        time = planets_data['1']['current_time']
        omega_b = planets_data['1']['longitude_of_perihelion']
        omega_c = planets_data['2']['longitude_of_perihelion']
        tmp = omega_c - omega_b
        ax1.plot(time, tmp)
        ax1.set_ylabel(field)
        ax1.set_xscale('log')
        # ax1.set_yscale('log')
        ax1.set_xlabel('time')
        # ax1.legend(loc=0, prop={'size':8})
        # plt.yticks(fontsize=12)
        # plt.xticks(fontsize=12)

        #-----------------------------------------------------------------------------------
        plt.tight_layout()
        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_orbit_ecc_data" +".png")
        plt.savefig(output_figure_filename, dpi=200)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))
    if plot_3e:# plot ecc data
        lab={}
        lab[1]='T1-b'
        lab[2]='T1-c'
        lab[3]='T1-d'
        lab[4]='T1-e'
        lab[5]='T1-f'
        lab[6]='T1-g'
        lab[7]='T1-h'
        largeur = 8*4
        hauteur = 7*2
        fig = plt.figure(figsize=(largeur,hauteur))
        ligne = 2
        colonne = 4
        i = 0
        #-----------------------------------------------------------------------------------
        
        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'ecc'
        aa=1.0
        for key in planets_keys:
            planet_data = planets_data[key]
            #tmp = 9.0
            time = planet_data['current_time']
            ecc = planet_data['eccentricity']
            ax1.plot(time, ecc,alpha=aa, label=lab[int(key)])
            aa-=0.1
        ax1.set_ylabel(field)
        # ax1.set_xscale('log')
        # ax1.set_yscale('log')
        ax1.set_xlabel('years')
        # ax1.legend(loc=0, prop={'size':8})

        
        field = 'ecc'
        aa=1.0
        for key in planets_keys:
            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            planet_data = planets_data[key]
            time = planet_data['current_time']
            ecc = planet_data['eccentricity']
            ax1.plot(time, ecc, label=lab[int(key)])
            ax1.set_xlabel('years')
        ax1.set_ylabel(field)
        # ax1.set_xscale('log')
        # ax1.set_yscale('log')
        ax1.legend(loc=0)

        #-----------------------------------------------------------------------------------
        plt.tight_layout()
        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_orbit_ecc_eachplanet_data" +".png")
        plt.savefig(output_figure_filename, dpi=200)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_3ebis:# plot ecc data
        lab={}
        lab[1]='T1-b'
        lab[2]='T1-c'
        lab[3]='T1-d'
        lab[4]='T1-e'
        lab[5]='T1-f'
        lab[6]='T1-g'
        lab[7]='T1-h'
        fig = plt.figure(figsize=(8, 10))
        ligne = 3
        colonne = 1
        i = 0
        #-----------------------------------------------------------------------------------
        
        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'eccentricity'
        aa=1.0
        for key in planets_keys:
            planet_data = planets_data[key]
            #tmp = 9.0
            time = planet_data['current_time']
            ecc = planet_data['eccentricity']
            ax1.plot(time, ecc,label=lab[int(key)])
            aa-=0.1
        ax1.set_ylabel(field)
        ax1.set_xscale('log')
        # ax1.set_yscale('log')
        ax1.set_xlabel('time (years)')
        ax1.legend(loc=0, prop={'size':8})

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = '$\omega$ (rad)'
        aa=1.0
        for key in planets_keys:
            planet_data = planets_data[key]
            time = planet_data["current_time"]
            longitude_peri = planet_data['longitude_of_perihelion']
            ax1.plot(time,longitude_peri, '+', alpha=aa)
            aa-=0.1
        ax1.set_ylabel('$\omega$ (rad)')
        # ax1.set_xscale('log')
        y_ticks = [0., np.pi, 2.0*np.pi ,] 
        y_tick_labels = [ '0', '$\pi$', '2$\pi$'  ]
        plt.yticks(y_ticks, y_tick_labels)
        ax1.set_xlabel('time (years)')
        # ax1.legend(loc=0, prop={'size':8})

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = '$\omega_c$ - $\omega_b$'
        time = planets_data['1']['current_time']
        omega_b = (planets_data['1']['longitude_of_perihelion']+1.0*np.pi)%(2.0*np.pi)
        omega_c = (planets_data['2']['longitude_of_perihelion']+1.0*np.pi)%(2.0*np.pi)
        tmp = (omega_c - omega_b)
        ax1.plot(time, tmp)
        ax1.set_ylabel(field)
        ax1.set_xscale('log')
        # ax1.set_yscale('log')
        ax1.set_xlabel('time (years)')
        y_ticks = [-2.0*np.pi, -np.pi , 0., np.pi, 2.0*np.pi] 
        y_tick_labels = ['-2 $\pi$', '-$\pi$', '0', '$\pi$', '2$\pi$', ]
        plt.yticks(y_ticks, y_tick_labels)
        # ax1.legend(loc=0, prop={'size':8})
        # plt.yticks(fontsize=12)
        # plt.xticks(fontsize=12)

        #-----------------------------------------------------------------------------------
        plt.tight_layout()
        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_orbit_ecc_data_bis" +".png")
        plt.savefig(output_figure_filename, dpi=200)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_3s:# plot spin data
        lab={}
        lab[1]='T1-b'
        lab[2]='T1-c'
        lab[3]='T1-d'
        lab[4]='T1-e'
        lab[5]='T1-f'
        lab[6]='T1-g'
        lab[7]='T1-h'
        fig = plt.figure(figsize=(15, 14))
        ligne = 4
        colonne = 2
        i = 0
        #-----------------------------------------------------------------------------------

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'Rotation state $\Omega$/n'
        aa = 1.0
        for key in planets_keys:
            planet_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]
            time = planet_data['current_time']
            omega = planet_computed_data['omega']
            mean_motion = planet_computed_data['mean_motion'] # rad/day
            ax1.plot(time, omega/mean_motion,  label='planet '+key, alpha=aa )
            aa-=0.10
        ax1.set_ylabel(field, fontsize=16)
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16)
        ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})
        
        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'spin x'
        aa=1.0
        for key in planets_keys:
            planet_data = planets_data[key]
            time = planet_data['current_time']
            sx = planet_data['spin_x']
            ax1.plot(time, sx, alpha=aa)
            aa-=0.1
        ax1.set_ylabel(field, fontsize=16)
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16)
        ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'spin y'
        aa=1.0
        for key in planets_keys:
            planet_data = planets_data[key]
            time = planet_data['current_time']
            sx = planet_data['spin_x']
            ax1.plot(time, sx, alpha=aa)
            aa-=0.1
        ax1.set_ylabel(field, fontsize=16)
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16)
        ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'spin z'
        aa=1.0
        for key in planets_keys:
            planet_data = planets_data[key]
            time = planet_data['current_time']
            sx = planet_data['spin_x']
            ax1.plot(time, sx, alpha=aa)
            aa-=0.1
        ax1.set_ylabel(field, fontsize=16)
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16)
        ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'Obliquity'
        for key in planets_keys:
            planet_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]
            time = planet_data['current_time']
            orbital_inclination = [p*180./np.pi for p in planet_data['inclination']]
            obliquity =  planet_computed_data['planet_obliquity']
            ls, =ax1.plot(time, obliquity, alpha=aa, label=f'Obliqu init: {round(obliquity[0],3)} | final: {round(obliquity[-1],3)}')
            # ax1.plot(time, orbital_inclination, ls=':', alpha=aa, color=ls.get_color(),label=f'Inclin init: {round(orbital_inclination[0],3)} | final: {round(orbital_inclination[-1],3)}')
        ax1.set_ylabel(field, fontsize=16)
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16)
        ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8}) #

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'Orbital inclination'
        for key in planets_keys:
            planet_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]
            time = planet_data['current_time']
            orbital_inclination = [p*180./np.pi for p in planet_data['inclination']]
            obliquity =  planet_computed_data['planet_obliquity']
            # ls, =ax1.plot(time, obliquity, alpha=aa, label=f'Obliqu init: {round(obliquity[0],3)} | final: {round(obliquity[-1],3)}')
            ax1.plot(time, orbital_inclination, ls=':', alpha=aa,label=f'Inclin init: {round(orbital_inclination[0],3)} | final: {round(orbital_inclination[-1],3)}')
        ax1.set_ylabel(field, fontsize=16)
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16)
        ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8}) #

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'Precession angle'
        for key in planets_keys:
            planet_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]
            time = planet_data['current_time']
            precession_angle = planet_computed_data['planet_precession_angle']
            ax1.plot(time, precession_angle,)
        ax1.set_ylabel(field, fontsize=16)
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16)
        # ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'Longitude ascending node'
        for key in planets_keys:
            planet_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]
            time = planet_data['current_time']
            long_asc_node = planet_data['longitude_of_ascending_node']
            ax1.plot(time, long_asc_node,)
        ax1.set_ylabel(field, fontsize=16)
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16)
        # ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        #-----------------------------------------------------------------------------------
        plt.tight_layout()
        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_spin_data" +".png")
        plt.savefig(output_figure_filename, dpi=200)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_30:# plot spin/mean_motion
        lab={}
        lab[1]='T1-b'
        lab[2]='T1-c'
        lab[3]='T1-d'
        lab[4]='T1-e'
        lab[5]='T1-f'
        lab[6]='T1-g'
        lab[7]='T1-h'
        fig = plt.figure(figsize=(10, 8))
        ligne = 1
        colonne = 1
        i = 0
        #-----------------------------------------------------------------------------------

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = 'Rotation state $\Omega$/n'
        tmp = planets_computed_data['1']['planet_rotation_period']
        aa = 1.0
        for key in planets_keys:
            # if key == '4':
            planet_data = planets_data[key]
            time = planet_data['current_time']
            period = planets_computed_data[key]['planet_rotation_period']
            omega = 2.0*np.pi / period
            mean_motion = planets_computed_data[key]['mean_motion'] # rad/day
            rr0 = period[0] / tmp[0]
            rr1 = period[len(period)-1] / tmp[len(tmp)-1]
            to_label = lab[int(key)]
            # to_label='P: %5.5f'%(period[0])+'d , ratio > '+'%5.5f'%((rr0))+ ' || %5.5f'%((rr1))
            ax1.plot(time, omega/mean_motion, label=to_label , alpha=aa )
            aa-=0.10
            tmp = period
        ax1.axhline(y=1., ls=':', color='black')
        ax1.set_ylabel(field)
        ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        #-----------------------------------------------------------------------------------
        plt.tight_layout()
        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_spin_meanmotion_data" +".png")
        plt.savefig(output_figure_filename, dpi=200)
        # plt.show()
        plt.close()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_30:# plot spin/mean_motion
        lab={}
        lab[1]='T1-b'
        lab[2]='T1-c'
        lab[3]='T1-d'
        lab[4]='T1-e'
        lab[5]='T1-f'
        lab[6]='T1-g'
        lab[7]='T1-h'
        fig = plt.figure(figsize=(15, 8*len(planets_keys)))
        ligne = 7
        colonne = 2
        i = 0
        #-----------------------------------------------------------------------------------
        for key in planets_keys:
            planet_data = planets_data[key]
            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Rotation state $\Omega$/n'
            time = planet_data['current_time']
            period = planets_computed_data[key]['planet_rotation_period']
            omega = 2.0*np.pi / period
            mean_motion = planets_computed_data[key]['mean_motion'] # rad/day
            to_label = lab[int(key)]
            ax1.plot(time, omega/mean_motion, label=to_label + " spin / mean motion " )
            ax1.axhline(y=1., ls=':', color='black')
            ax1.set_ylabel(field)
            # ax1.set_yscale('log')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Rotation state $\Omega$ vs n'
            ax1.plot(time, mean_motion, label= to_label + " mean motion" )
            ax1.plot(time, omega, label=to_label + " spin " )
            ax1.set_ylabel(field)
            # ax1.set_yscale('log')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

        #-----------------------------------------------------------------------------------
        plt.tight_layout()
        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_spin_meanmotion_eachplanet_data" +".png")
        plt.savefig(output_figure_filename, dpi=200)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_31:# Plot sub-stellar point 
        lab={}
        lab[1]='T1-b'
        lab[2]='T1-c'
        lab[3]='T1-d'
        lab[4]='T1-e'
        lab[5]='T1-f'
        lab[6]='T1-g'
        lab[7]='T1-h'
        fig = plt.figure(figsize=(10, 20))
        ligne = 2
        colonne = 1
        #-----------------------------------------------------------------------------------

        fig, ax = plt.subplots(ligne, colonne, figsize=(15, 15))
        tmp = [ 00. for i in np.arange(7) ]
        print(tmp)
        table_times = {}
        table_times['planet'] = []
        table_times['time'] = []
        table_times['sub'] = []
        table_times['planet_integr'] = []
        table_times['time_integr'] = []
        table_times['sub_integr'] = []
        for key in planets_keys:
            planet_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]
            # time = [t*365.25 for t in planet_data['current_time']]
            time = planet_data['current_time']
            sub_stellar = planet_computed_data['substellar_point']
            x = time
            y = sub_stellar
            P_syn = round(median(planet_computed_data['synodic_period']),3 )
            P_syn_mean =round(median(planet_computed_data['synodic_period_mean']),3 )

            coefficients = np.polyfit(x, y, 1)
            print(' > For planet',key, f' > Computed Coef {coefficients[0]} x + {coefficients[1]}  and time scale {360./coefficients[0]}' , )
            poly_function = np.poly1d(coefficients)
            y_fit = poly_function(x)
            
            # in_label =lab[int(key)] +f" P_synodic : {round(median(planet_computed_data['synodic_period']), 3)} | end val : {round(sub_stellar[-1], 3)}"  {round(median(planet_computed_data['synodic_period']), 3)} | {round(median(planet_computed_data['synodic_period_mean']), 3)}
            in_label =lab[int(key)] + f' | y = {round(coefficients[0],3)} x + {round(coefficients[1],3)} | time-scale {round(360./coefficients[0],3)} | End value {round(sub_stellar[-1],3)}  '
            ls, = ax[0].plot(time, sub_stellar,label=in_label)
            ax[0].plot(x, y_fit, '-.', color=ls.get_color(), linewidth=2)
            for t,s in zip(time, sub_stellar):
                if np.abs(s) > tmp[int(key)-1]:
                    print("Planet: ", key, t, tmp[int(key)-1], np.abs(s) )
                    ax[0].axhline(y=s, ls=':', color='black')
                    ax[0].axvline(x=t, ls=':', color=ls.get_color())
                    tmp[int(key)-1]+=360.
                    table_times['planet'].append(key)
                    table_times['time'].append(t)
                    table_times['sub'].append(s)
        ax[0].axhline(y=0., ls=':', color='black')
        # ax[0].axhline(y=180., ls=':', color='black')
        # ax[0].axhline(y=360., ls=':', color='black')
        columns = ('Time', 'angle', 'a x', '+b')
        rows = [ lab[int(in_table)] for in_table in table_times['planet'] ]
        table = np.transpose((table_times['time'], table_times['sub']))
        the_table = ax[1].table(cellText=table, loc='center', colLabels=columns, rowLabels=rows)
        ax[1].axis('tight')
        ax[1].axis('off')

        ax[0].set_ylabel('Sub-stellar point position (degrees) ', fontsize=15)
        ax[0].set_xlabel('Time (years)', fontsize=15)
        ax[0].tick_params(axis='both', labelsize=15)
        ax[0].legend(loc=0, prop={'size':12})

        
        #-----------------------------------------------------------------------------------
        # plt.tight_layout()
        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_substellar_data" +".png")
        plt.savefig(output_figure_filename, dpi=200)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_31:# Plot sub-stellar point for paper
        lab={}
        lab[1]='T1-b'
        lab[2]='T1-c'
        lab[3]='T1-d'
        lab[4]='T1-e'
        lab[5]='T1-f'
        lab[6]='T1-g'
        lab[7]='T1-h'
        fig = plt.figure(figsize=(12, 8))
        ligne = 1
        colonne = 1
        #-----------------------------------------------------------------------------------
        i=1
        # fig, ax = plt.subplots(ligne, colonne, figsize=(15, 15))
        ax = fig.add_subplot(ligne,colonne,i,)
        tmp = [ 00. for i in np.arange(7) ]

        for key in planets_keys:
            planet_data = planets_data[key]
            # time = [t*365.25 for t in planet_data['current_time']]
            time = planet_data['current_time']
            # sub_stellar = [ ((np.pi-p*np.pi/180.)%(2.*np.pi) -np.pi ) for p in planets_computed_data[key]['substellar_point']]
            sub_stellar = planets_computed_data[key]['substellar_point']
            # x = time
            # y = sub_stellar

            # coefficients = np.polyfit(x, y, 1)
            # print(' > For planet',key, f' > Computed Coef {coefficients[0]} x + {coefficients[1]}  and time scale {360./coefficients[0]}' , )
            # poly_function = np.poly1d(coefficients)
            # y_fit = poly_function(x)
            
            # in_label =lab[int(key)] +f" P_synodic : {round(median(planets_computed_data[key]['synodic_period']), 3)} | end val : {round(sub_stellar[-1], 3)}"
            # in_label = lab[int(key)] + f' | y = {round(coefficients[0],3)} x + {round(coefficients[1],3)} | time-scale {round(360./coefficients[0],3)}'
            in_label = lab[int(key)] 
            ls, = ax.plot(time, sub_stellar,label=in_label)
            # ax.plot(x, y_fit, '-.', color=ls.get_color(), linewidth=2)
            for t,s in zip(time, sub_stellar):
                if np.abs(s) > tmp[int(key)-1]:
                    print("Planet: ", key, t, tmp[int(key)-1], np.abs(s) )
                    ax.axhline(y=s, ls=':', color='black')
                    # ax.axvline(x=t, ls=':', color=ls.get_color())
                    tmp[int(key)-1]+=360.

        ax.axhline(y=0., ls=':', color='black')
        ax.set_ylabel('Sub-stellar point position (degrees) ', fontsize=15)
        ax.set_xlabel('Time (years)', fontsize=15)
        ax.tick_params(axis='both', labelsize=15)
        ax.legend(loc=0, prop={'size':12})
        
        #-----------------------------------------------------------------------------------
        # plt.tight_layout()
        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_substellar_data_paper" +".png")
        plt.savefig(output_figure_filename, dpi=200)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_31b:# Plot sub-stellar point for b
        lab={}
        lab[1]='T1-b'
        lab[2]='T1-c'
        lab[3]='T1-d'
        lab[4]='T1-e'
        lab[5]='T1-f'
        lab[6]='T1-g'
        lab[7]='T1-h'
        fig = plt.figure(figsize=(10, 20))
        ligne = 2
        colonne = 1
        #-----------------------------------------------------------------------------------

        fig, ax = plt.subplots(ligne, colonne, figsize=(15, 15))
        tmp = [ 00. for i in np.arange(7) ]
        print(tmp)
        table_times = {}
        table_times['planet'] = []
        table_times['time'] = []
        table_times['sub'] = []
        table_times['planet_integr'] = []
        table_times['time_integr'] = []
        table_times['sub_integr'] = []
        for key in planets_keys:
            if key =='1':
                planet_data = planets_data[key]
                # time = [t*365.25 for t in planet_data['current_time']]
                time = planet_data['current_time']
                sub_stellar = planets_computed_data[key]['substellar_point']
                x = time
                y = sub_stellar

                coefficients = np.polyfit(x, y, 1)
                print(' > For planet',key, f' > Computed Coef {coefficients[0]} x + {coefficients[1]}  and time scale {360./coefficients[0]}' , )
                poly_function = np.poly1d(coefficients)
                y_fit = poly_function(x)
                
                # in_label =lab[int(key)] +f" P_synodic : {round(median(planets_computed_data[key]['synodic_period']), 3)} | end val : {round(sub_stellar[-1], 3)}"
                in_label =lab[int(key)]
                ls, = ax[0].plot(time, sub_stellar,label=in_label)
                ax[0].plot(x, y_fit, '-.', color=ls.get_color(), linewidth=2)
                for t,s in zip(time, sub_stellar):
                    if np.abs(s) > tmp[int(key)-1]:
                        print("Planet: ", key, t, tmp[int(key)-1], np.abs(s) )
                        ax[0].axhline(y=s, ls=':', color='black')
                        ax[0].axvline(x=t, ls=':', color=ls.get_color())
                        tmp[int(key)-1]+=90.
                        table_times['planet'].append(key)
                        table_times['time'].append(t)
                        table_times['sub'].append(s)
        ax[0].axhline(y=0., ls=':', color='black')
        # ax[0].axhline(y=180., ls=':', color='black')
        # ax[0].axhline(y=360., ls=':', color='black')
        columns = ('Time', 'angle')
        rows = [ lab[int(in_table)] for in_table in table_times['planet'] ]
        table = np.transpose((table_times['time'], table_times['sub']))
        the_table = ax[1].table(cellText=table, loc='center', colLabels=columns, rowLabels=rows)
        ax[1].axis('tight')
        ax[1].axis('off')

        ax[0].set_ylabel('Sub-stellar point position (degrees) ', fontsize=15)
        ax[0].set_xlabel('Time (years)', fontsize=15)
        ax[0].tick_params(axis='both', labelsize=15)
        ax[0].legend(loc=0, prop={'size':12})

        
        #-----------------------------------------------------------------------------------
        # plt.tight_layout()
        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_substellar_data_b" +".png")
        plt.savefig(output_figure_filename, dpi=200)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_31all:# Plot sub-stellar point for each planet 
        lab={}
        lab[1]='T1-b'
        lab[2]='T1-c'
        lab[3]='T1-d'
        lab[4]='T1-e'
        lab[5]='T1-f'
        lab[6]='T1-g'
        lab[7]='T1-h'

        #-----------------------------------------------------------------------------------

        fig = plt.figure(figsize=(20, 20))
        ligne = 3
        colonne = 3
        i = 0
        #-----------------------------------------------------------------------------------
        
        for key in planets_keys:
            tmp = [ 00. for i in np.arange(7) ]
            i+=1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Sub-stellar points'
            planet_data = planets_data[key]
            time = [t*365.25 for t in planet_data['current_time']]
            sub_stellar = planets_computed_data[key]['substellar_point']
            in_label =lab[int(key)] +f" P_synodic : {round(median(planets_computed_data[key]['synodic_period']), 3)} | end val : {round(sub_stellar[-1], 3)}"
            ls, = ax1.plot(planet_data['current_time'], sub_stellar,label=in_label)

            ax1.set_ylabel(' degrees')
            # ax1.set_yscale('log')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})
            ax1.xaxis.set_minor_formatter(FormatStrFormatter('% 1.3f'))

            for t,s in zip(time, sub_stellar):
                if np.abs(s) > tmp[int(key)-1]:
                    # print("Planet: ", key, t, tmp[int(key)-1], np.abs(s) )
                    ax1.axhline(y=s, ls=':', color='black')
                    ax1.axvline(x=t/365.25, ls=':', color=ls.get_color())
                    tmp[int(key)-1]+=90.
            ax1.set_ylabel('Sub-stellar point position (degrees) ', fontsize=15)
            ax1.set_xlabel('Time (years)', fontsize=15)
            ax1.tick_params(axis='both', labelsize=15)
            ax1.legend(loc=0, prop={'size':12})
        
        #-----------------------------------------------------------------------------------
        # plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_substellar_data_every_planet" +".png")
        plt.savefig(output_figure_filename)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_32:# Plot sub-stellar point 2 calculs
        lab={}
        lab[1]='T1-b'
        lab[2]='T1-c'
        lab[3]='T1-d'
        lab[4]='T1-e'
        lab[5]='T1-f'
        lab[6]='T1-g'
        lab[7]='T1-h'
        fig = plt.figure(figsize=(30, 12))
        ligne = 4
        colonne = 1
        p = 0
        #-----------------------------------------------------------------------------------

        p+=1
        # ax = fig.add_subplot(ligne,colonne,p,)
        # fig2, ax = fig.add_subplot(ligne,colonne,)
        fig, ax = plt.subplots(ligne, colonne, figsize=(15, 15))
        tmp = [ 0. for i in np.arange(7) ]
        tmp_integ = [ 0. for i in np.arange(7) ]
        print(tmp)
        # table_times = {}
        table_times = {}
        table_times['planet'] = []
        table_times['time'] = []
        table_times['sub'] = []
        table_times['planet_integr'] = []
        table_times['time_integr'] = []
        table_times['sub_integr'] = []
        for key in planets_keys:
            # if key == '2':
            # print('   evaluate ', key)
            planet_data = planets_data[key]
            time = [t*365.25 for t in planet_data['current_time']]
            sub_stellar = planets_computed_data[key]['substellar_point']
            in_label =lab[int(key)] +f" P_synodic : {round(median(planets_computed_data[key]['synodic_period']), 3)} | end val : {round(sub_stellar[-1], 3)}"
            ls, = ax[0].plot(time, sub_stellar,label=in_label)
            
            for t,s in zip(time, sub_stellar):
                if np.abs(s) > tmp[int(key)-1]:
                    # print("Planet: ", key, t, tmp[int(key)-1], np.abs(s) )
                    ax[0].axhline(y=s, ls=':', color='black')
                    ax[0].axvline(x=t, ls=':', color=ls.get_color())
                    tmp[int(key)-1]+=360.
                    table_times['planet'].append(key)
                    table_times['time'].append(t)
                    table_times['sub'].append(s)

            sub_stellar_integrated = [ p for p in planets_computed_data[key]['theta_true_ano_long_peri'] ] # [ p *180./np.pi for p in planets_computed_data[key]['theta_true_ano_long_peri'] ] 
            ls2, = ax[2].plot(time, sub_stellar_integrated, label=key) 
            for t,s in zip(time, sub_stellar_integrated):
                if np.abs(s) > tmp_integ[int(key)-1]:
                    ax[2].axhline(y=s, ls=':', color='black')
                    ax[2].axvline(x=t, ls=':', color=ls.get_color())
                    tmp_integ[int(key)-1]+=360.
                    table_times['planet_integr'].append(key)
                    table_times['time_integr'].append(t)
                    table_times['sub_integr'].append(s)

        columns = ('Time', 'angle')
        rows = [ lab[int(in_table)] for in_table in table_times['planet'] ]
        table = np.transpose((table_times['time'], table_times['sub']))
        the_table = ax[1].table(cellText=table, loc='center', colLabels=columns, rowLabels=rows)
        ax[1].axis('tight')
        ax[1].axis('off')

        rows = [ lab[int(in_table)] for in_table in table_times['planet_integr'] ]
        table = np.transpose((table_times['time_integr'], table_times['sub_integr']))
        the_table = ax[3].table(cellText=table, loc='center', colLabels=columns, rowLabels=rows)
        ax[3].axis('tight')
        ax[3].axis('off')

        ax[0].set_ylabel('Sub-stellar point position (degrees) - Cartesian ', fontsize=15)
        ax[0].set_xlabel('Time (days)', fontsize=15)
        ax[0].tick_params(axis='both', labelsize=15)
        ax[0].legend(loc=0, prop={'size':12})

        ax[2].set_ylabel('Sub-stellar point position (degrees) - True anomaly ', fontsize=15)
        ax[2].set_xlabel('Time (days)', fontsize=15)
        ax[2].tick_params(axis='both', labelsize=15)
        ax[2].legend(loc=0, prop={'size':12})
        
        #-----------------------------------------------------------------------------------
        # plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_substellar_data_2calculs" +".png")
        plt.savefig(output_figure_filename)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_33: # Check orbital angles 

        lab={}
        lab[1]='T1-b'
        lab[2]='T1-c'
        lab[3]='T1-d'
        lab[4]='T1-e'
        lab[5]='T1-f'
        lab[6]='T1-g'
        lab[7]='T1-h'

        
        fig = plt.figure(figsize=(25, 6*len(planets_keys)))
        ligne = len(planets_keys)
        colonne = 4
        i = 0
        #-----------------------------------------------------------------------------------

        for key in planets_keys:

            planet_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]

            time = [t*365.25 for t in planet_data['current_time']]
            true_ano =  planet_data['true_anomaly']
            mean_ano =  planet_data['mean_anomaly']
            long_peri = planet_data['longitude_of_perihelion']

            mean_motion = planet_computed_data['mean_motion']
            mean_longitude = planet_computed_data['mean_longitude']
            true_ano_cumul = planet_computed_data['true_anomaly_cumul']
            cartesian_angle_cumul = [p*np.pi/180. for p in planet_computed_data['cartesian_position_angle']]
            theta = planet_computed_data['theta']
            mean_ano_cumul = planet_computed_data['mean_ano_cumul']
            theta_true_ano = planet_computed_data['theta_true_ano'] 
            long_peri_rescaled = planet_computed_data['longitude_of_perihelion_rescaled'] 
            theta_true_ano_long_peri = planet_computed_data['theta_true_ano_long_peri']
            theta_mean_ano = planet_computed_data['theta_mean_ano'] 
            theta_cartesian_position_angle = planet_computed_data['theta_cartesian_position_angle']

            i+=1 # --- 1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = '  '+lab[int(key)]
            ax1.plot(time, true_ano,label=f'true anomaly')
            ax1.plot(time, true_ano_cumul,label=f'true anomaly cumul')
            ax1.plot(time, cartesian_angle_cumul, '+',label=f'cartesian angle cumul')
            ax1.axhline(y=0., ls=':', color='black')
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('days')
            ax1.legend(loc=0, prop={'size':8})

            i+=1 # --- 1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = '  '+lab[int(key)]
            ax1.plot(time, theta,label=f'theta')
            t0 = true_ano_cumul[0]
            true_ano_cumul_0 = [t-t0 for t in true_ano_cumul]
            ax1.plot(time, true_ano_cumul_0,'+',label=f'true anomaly cumul (init 0)')
            ax1.plot(time, mean_ano_cumul,label=f'mean_ano_cumul ')
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('days')
            ax1.legend(loc=0, prop={'size':8})

            i+=1 # --- 1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = '  '+lab[int(key)]
            ax1.plot(time, theta_true_ano,label=f'theta true anomaly ')
            ax1.plot(time, theta_mean_ano,'+',label=f'theta mean ano')
            ax1.plot(time, theta_cartesian_position_angle,'+',label=f'theta cartesian position angle')
            ax1.plot(time, theta_true_ano_long_peri, label=f'theta threano-longperi corrected')
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('days')
            ax1.legend(loc=0, prop={'size':8})

            i+=1 # --- 1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = '  '+lab[int(key)]
            tmp = [t-c for t,c in zip(theta_true_ano,theta_cartesian_position_angle)]
            ax1.plot(time, tmp,label=f'true anomaly - cartesian position angle  ')
            ax1.plot(time, long_peri_rescaled,label=f'long peri')
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('days')
            ax1.legend(loc=0, prop={'size':8})
        
        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_true_mean_cartesian_angles"+".png")
        plt.savefig(output_figure_filename)
        # plt.show()
        plt.close()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_33_1p: # Check true anomaly

        lab={}
        lab[1]='T1-b'
        lab[2]='T1-c'
        lab[3]='T1-d'
        lab[4]='T1-e'
        lab[5]='T1-f'
        lab[6]='T1-g'
        lab[7]='T1-h'

        
        fig = plt.figure(figsize=(20, 6*len(planets_keys)))
        ligne = len(planets_keys)
        colonne = 3
        i = 0
        #-----------------------------------------------------------------------------------

        # planet_period = (18.772866/365.25)

        for key in planets_keys:

            planet_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]

            time = [t*365.25 for t in planet_data['current_time']]
            true_ano = planet_data['true_anomaly']
            mean_ano =  planet_data['mean_anomaly']

            mean_motion = planet_computed_data['mean_motion']
            mean_longitude = planet_computed_data['mean_longitude']
            true_ano_cumul = planet_computed_data['true_anomaly_cumul']
            cartesian_angle_cumul = [p*np.pi/180. for p in planet_computed_data['cartesian_position_angle']]
            theta = planet_computed_data['theta']
            theta_true_ano = planet_computed_data['theta_true_ano'] 
            theta_mean_ano = planet_computed_data['theta_mean_ano'] 

            i+=1 # --- 1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = '  '+lab[int(key)]
            ax1.plot(time,true_ano,label=f'true anomaly')
            # ax1.plot(time,mean_ano,label=f'mean anomaly')
            # ax1.plot(time,true_ano_cumul,label=f'true anomaly cumul')
            # ax1.plot(time,cartesian_angle_cumul,label=f'cartesian anomaly cumul')
            # ax1.axvline(x=planet_period, ls=':', color='black')
            ax1.axhline(y=0., ls=':', color='black')
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('days')
            ax1.legend(loc=0, prop={'size':8})

            i+=1 # --- 1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = '  '+lab[int(key)]
            # ax1.plot(time,true_ano,label=f'true anomaly')
            # ax1.plot(time,mean_ano,label=f'mean anomaly')
            ax1.plot(time,true_ano_cumul,label=f'true anomaly cumul')
            # ax1.plot(time,cartesian_angle_cumul,label=f'cartesian anomaly cumul')
            # ax1.axvline(x=planet_period, ls=':', color='black')
            ax1.axhline(y=0., ls=':', color='black')
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('days')
            ax1.legend(loc=0, prop={'size':8})

            i+=1 # --- 1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = '  '+lab[int(key)]
            # ax1.plot(time,true_ano,label=f'true anomaly')
            # ax1.plot(time,mean_ano,label=f'mean anomaly')
            # ax1.plot(time,true_ano_cumul,label=f'true anomaly cumul')
            ax1.plot(time,cartesian_angle_cumul,label=f'cartesian anomaly cumul')
            # ax1.axvline(x=planet_period, ls=':', color='black')
            ax1.axhline(y=0., ls=':', color='black')
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('days')
            ax1.legend(loc=0, prop={'size':8})
        
        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_true_anomaly_angles"+".png")
        plt.savefig(output_figure_filename)
        # plt.show()
        plt.close()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_33_meanano :  # Check mean anomaly

        lab={}
        lab[1]='T1-b'
        lab[2]='T1-c'
        lab[3]='T1-d'
        lab[4]='T1-e'
        lab[5]='T1-f'
        lab[6]='T1-g'
        lab[7]='T1-h'

        fig = plt.figure(figsize=(30, 6*len(planets_keys)))
        ligne = len(planets_keys)
        colonne = 2
        i = 0

        for key in planets_keys:

            planet_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]
            time = [t*365.25 for t in planet_data['current_time']]
            # true_ano = planet_data['true_anomaly']
            mean_ano = planet_data['mean_anomaly']
            long_peri =  planet_data['longitude_of_perihelion']

            # mean_motion = planet_computed_data['mean_motion']
            # mean_longitude = planet_computed_data['mean_longitude']
            # true_ano_cumul = planet_computed_data['true_anomaly_cumul']
            # cartesian_angle_cumul = [p*np.pi/180. for p in planet_computed_data['cartesian_position_angle']]
            # theta = planet_computed_data['theta']
            # theta_true_ano = planet_computed_data['theta_true_ano'] 
            # theta_mean_ano = planet_computed_data['theta_mean_ano'] 
            check_mean_ano = planet_computed_data['check_mean_anomaly']
            check_long_peri = planet_computed_data['check_longitude_of_perihelion']
            print(f' |\t {check_mean_ano[0]}')

            i+=1 # --- 1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = '  '+lab[int(key)]
            ax1.plot(time,mean_ano,label=f'mean anomaly')
            ax1.plot(time,check_mean_ano,label=f'check mean anomaly')
            ax1.axhline(y=0., ls=':', color='black')
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('days')
            ax1.legend(loc=0, prop={'size':8})

            i+=1 # --- 1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = '  '+lab[int(key)]
            ax1.plot(time,long_peri,label=f'Long peri')
            ax1.plot(time,check_long_peri,label=f'check long_peri')
            ax1.axhline(y=0., ls=':', color='black')
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('days')
            ax1.legend(loc=0, prop={'size':8})
        
        plt.tight_layout()
        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_mean_anomaly_angles"+".png")
        plt.savefig(output_figure_filename)
        # plt.show()
        plt.close()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_34: # plot eccentricities
        
        fig = plt.figure(figsize=(16, 16))
        ligne = 3
        colonne = 3
        i = 0
        #-----------------------------------------------------------------------------------

        for key in planets_keys:

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Spin rates'

            planet_data = planets_data[key]
            time = planet_data['current_time']
            omega = planets_computed_data[key]['omega']

            ax1.plot(time, omega, label='Planet '+key, alpha=1.)

            ax1.set_ylabel(' rad/days ')
            ax1.set_yscale('log')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})
            ax1.xaxis.set_minor_formatter(FormatStrFormatter('% 1.3f'))

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_compare_spin_rates" +".png")
        plt.savefig(output_figure_filename, dpi=200)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_4: # compare spin mean motion
        stored_data = {
            'eMax':0., 'eMin':0., 'eMed':0., 'eMean':0., 'OMax':0., 'OMin':0., 'OMed':0., 'OMean':0. ,'einit':0.,
            }
        
        fig = plt.figure(figsize=(16, 16))
        ligne = 3
        colonne = 3
        i = 0
        #-----------------------------------------------------------------------------------

        for key in planets_keys:

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Omega n'

            planet_data = planets_data[key]
            time = planet_data['current_time']
            ecc = planet_data['eccentricity']
            period = planets_computed_data[key]['planet_rotation_period']
            omega = 2.0*np.pi / period
            mean_motion = 2.0*np.pi / planets_computed_data[key]['orbital_period'] 

            ax1.plot(time, mean_motion, label='Mean motion n', alpha=1.)
            ax1.plot(time, omega,  label='Spin rate $\Omega$', alpha=1.)
            mini = min(mean_motion)
            maxi = max(mean_motion)
            med = median(mean_motion)
            mea = mean(mean_motion)
            print(key, mini, maxi, med, mea)
            # ax1.plot(time, ecc,  label='spin', alpha=0.5)
            # ax1.set_ylabel(field+' Planet: '+ str(key))
            ax1.set_ylabel('rad/day')
            ax1.set_yscale('log')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

            # stored_data['OMax']  = round(maxi, 4)
            # stored_data['OMin']  = round(mini, 4)
            # stored_data['OMed']  = round(med, 4)
            # stored_data['OMean'] = round(mea, 4)

            # file_name = "./target/t1_all_earth_1Kyrs/t1_all_earth_1Kyrs_stat.csv"

            # df = pd.DataFrame(stored_data)
            # df.to_csv(file_name, sep="\t" )

            ax1.xaxis.set_minor_formatter(FormatStrFormatter('% 1.3f'))

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_compare_spin_meanmotion" +".png")
        plt.savefig(output_figure_filename, dpi=200)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_41: # compare mean motion

        fig = plt.figure(figsize=(16, 16))
        ligne = 3
        colonne = 3
        i = 0
        #-----------------------------------------------------------------------------------

        for key in planets_keys:

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = ' ? '

            planet_data = planets_data[key]
            time = planet_data['current_time'][planet_data['current_time'] <(1000.+10./365.25)]

            mean_motion = 2.0*np.pi / (planets_computed_data[key]['orbital_period'][planet_data['current_time'] <(1000.+10./365.25)] *86400.)
            l = planet_data['mean_anomaly'][planet_data['current_time'] <(1000.+10./365.25)]

            ax1.plot(time, mean_motion, '+', label='mean motion', )
            # ax1.plot(time, l,  label='mean anomaly', )
            # ax1.plot(time, ecc,  label='spin', alpha=0.5)
            
            ax1.set_ylabel(field+' Planet: '+ str(key))
            # ax1.set_yscale('log')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_mean_motion" +".png")
        plt.savefig(output_figure_filename)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))
    
    if plot_42: # compare mean_anomaly true anomaly

        fig = plt.figure(figsize=(16, 16))
        ligne = 7
        colonne = 2
        i = 0
        #-----------------------------------------------------------------------------------

        for key in planets_keys:
            planet_data = planets_data[key]
            time = planet_data['current_time']
            mean_motion = 2.0*np.pi / (planets_computed_data[key]['orbital_period'] )
            l = planet_data['mean_anomaly']
            mean_ano = [(m *t) %(2.*np.pi) for (m,t) in zip(mean_motion, time)]
            true_anomaly =  planet_data['true_anomaly']

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = ' ? '
            ax1.plot(time, mean_ano, label=' from mean motion', )
            ax1.set_ylabel(field+' Planet: '+ str(key))
            # ax1.set_yscale('log')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = ' ? '
            planet_data = planets_data[key]
            ax1.plot(time, l, label='mean anomaly', )
            ax1.plot(time, true_anomaly, label='True anomaly', )
            ax1.set_ylabel(field+' Planet: '+ str(key))
            # ax1.set_yscale('log')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

            # i=i+1
            # ax1 = fig.add_subplot(ligne,colonne,i,)
            # field = ' ? '
            # planet_data = planets_data[key]
            
            # ax1.set_ylabel(field+' Planet: '+ str(key))
            # # ax1.set_yscale('log')
            # ax1.set_xlabel('years')
            # ax1.legend(loc=0, prop={'size':8})

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_mean_anomaly" +".png")
        plt.savefig(output_figure_filename)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_43: # plot eccentricities
        
        fig = plt.figure(figsize=(16, 16))
        ligne = 3
        colonne = 3
        i = 0
        #-----------------------------------------------------------------------------------

        for key in planets_keys:

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Ecc'

            planet_data = planets_data[key]
            time = planet_data['current_time']
            ecc = planet_data['eccentricity']

            ax1.plot(time, ecc, label='Planet '+key, alpha=1.)

            ax1.set_ylabel(' ecc ')
            ax1.set_yscale('log')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})
            ax1.xaxis.set_minor_formatter(FormatStrFormatter('% 1.3f'))

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_compare_eccentricities" +".png")
        plt.savefig(output_figure_filename, dpi=200)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_5:

        fig = plt.figure(figsize=(16, 16))
        ligne = 3
        colonne = 3
        i = 0
        #-----------------------------------------------------------------------------------

        for key in planets_keys:

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Orbital period'

            planet_data = planets_data[key]
            time = planet_data['current_time']
            # ecc = planet_data['eccentricity']
            period = planets_computed_data[key]['planet_rotation_period']
            omega = 2.0*np.pi / period
            mean_motion = 2.0*np.pi / planets_computed_data[key]['orbital_period'] 
            if key == '1':
                ax1.axhline(y=1.510826, ls=":", color='black' )
            if key == '2':
                ax1.axhline(y=2.421937, ls=":", color='black' )
            if key == '3':
                ax1.axhline(y=4.049219, ls=":", color='black' )
            if key == '4':
                ax1.axhline(y=6.101013, ls=":", color='black' )
            if key == '5':
                ax1.axhline(y=9.20754, ls=":", color='black' )
            if key == '6':
                ax1.axhline(y=12.352446, ls=":", color='black' )
            if key == '7':
                ax1.axhline(y=18.772866, ls=":", color='black' )

            ax1.plot(time, period, label='mean motion', alpha=0.5)
            # ax1.plot(time, omega,  label='spin', alpha=0.5)
            # ax1.plot(time, ecc,  label='spin', alpha=0.5)
        
            ax1.set_ylabel(field+' Planet: '+ str(key))
            ax1.set_yscale('log')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_compare_orbital_period_with_agol" +".png")
        plt.savefig(output_figure_filename)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))
 
    if plot_6:

        fig = plt.figure(figsize=(16, 16))
        ligne = 3
        colonne = 3
        i = 0
        #-----------------------------------------------------------------------------------

        for key in planets_keys:

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Semi-Major axis'

            planet_data = planets_data[key]
            time = planet_data['current_time']
            # ecc = planet_data['eccentricity']
            sma = planet_data['semi-major_axis']
            # period = planets_computed_data[key]['planet_rotation_period']
            # omega = 2.0*np.pi / period
            # mean_motion = 2.0*np.pi / planets_computed_data[key]['orbital_period'] 

            ax1.plot(time, sma, label='sma planet '+key)
            # ax1.plot(time, omega,  label='spin', alpha=0.5)
            # ax1.plot(time, ecc,  label='spin', alpha=0.5)
        
            ax1.set_ylabel(field+' Planet: '+ str(key))
            # ax1.set_yscale('log')
            # ax1.set_xscale('log')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_compare_sma" +".png")
        plt.savefig(output_figure_filename)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_7:

        fig = plt.figure(figsize=(6, 6))
        ligne = 1
        colonne = 1
        i = 0
        #-----------------------------------------------------------------------------------

        # for key in planets_keys:
        p1 = planets_computed_data['1']['planet_rotation_period']
        p2 = planets_computed_data['2']['planet_rotation_period']
        p3 = planets_computed_data['3']['planet_rotation_period']
        p4 = planets_computed_data['4']['planet_rotation_period']
        p5 = planets_computed_data['5']['planet_rotation_period']
        p6 = planets_computed_data['6']['planet_rotation_period']
        p7 = planets_computed_data['7']['planet_rotation_period']

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = ' compare_period_ratio_and_observed_offset '

        planet_data = planets_data[key]
        t1 = planets_data['1']['current_time']
        t2 = planets_data['2']['current_time']
        t3 = planets_data['3']['current_time']
        t4 = planets_data['4']['current_time']
        t5 = planets_data['5']['current_time']
        t6 = planets_data['6']['current_time']
        t7 = planets_data['7']['current_time']
 
        ls, =ax1.plot(t1, p2/p1, label='p2/p1')
        ax1.axhline(y=1.603055, ls=":", color=ls.get_color(), label='oberved')
        ls, =ax1.plot(t2, p3/p2, label='p3/p2')
        ax1.axhline(y=1.671893, ls=":", color=ls.get_color() )
        ls, =ax1.plot(t3, p4/p3, label='p4/p5')
        ax1.axhline(y=1.506714, ls=":", color=ls.get_color() )
        ls, =ax1.plot(t4, p5/p4, label='p5/p4')
        ax1.axhline(y=1.509182, ls=":", color=ls.get_color() )
        ls, =ax1.plot(t5, p6/p5, label='p6/p5')
        ax1.axhline(y=1.341558, ls=":",color=ls.get_color() )
        ls, =ax1.plot(t6, p7/p6, label='p7/p6')
        ax1.axhline(y=1.519769, ls=":",color=ls.get_color() )

        # off = p2/p1 - 1.603055
        # ls, =ax1.plot(t1,off , label='p2/p1')
        # # ax1.axhline(y=1.603055 -3./2., ls=":", color=ls.get_color() )
        # off = p3/p2 - 1.67189
        # ls, =ax1.plot(t2, off, label='p3/p2')
        # # ax1.axhline(y=- 5./3., ls=":", color=ls.get_color() )
        # off = p4/p3 -1.506714
        # ls, =ax1.plot(t3, off, label='p4/p5')
        # # ax1.axhline(y=1.506714- 3./2., ls=":", color=ls.get_color() )
        # off = p5/p4 -1.509182
        # ls, =ax1.plot(t4, off, label='p5/p4')
        # # ax1.axhline(y=1.509182 - 3./2., ls=":", color=ls.get_color() )
        # off = p6/p5 - 1.341558
        # ls, =ax1.plot(t5, off, label='p6/p5')
        # # ax1.axhline(y=1.341558 - 4./3., ls=":",color=ls.get_color() )
        # off = p7/p6 -1.519769
        # ls, =ax1.plot(t6, off, label='p7/p6')
        # # ax1.axhline(y=1.519769 - 3./2., ls=":",color=ls.get_color() )

        # off = p2/p1 -3./2.
        # ls, =ax1.plot(t1,off , label='p2/p1')
        # ax1.axhline(y=1.603055 -3./2., ls=":", color=ls.get_color() )
        # off = p3/p2 - 5./3.
        # ls, =ax1.plot(t2, off, label='p3/p2')
        # ax1.axhline(y=1.67189- 5./3., ls=":", color=ls.get_color() )
        # off = p4/p3 - 3./2.
        # ls, =ax1.plot(t3, off, label='p4/p5')
        # ax1.axhline(y=1.506714- 3./2., ls=":", color=ls.get_color() )
        # off = p5/p4 - 3./2.
        # ls, =ax1.plot(t4, off, label='p5/p4')
        # ax1.axhline(y=1.509182 - 3./2., ls=":", color=ls.get_color() )
        # off = p6/p5 - 4./3.
        # ls, =ax1.plot(t5, off, label='p6/p5')
        # ax1.axhline(y=1.341558 - 4./3., ls=":",color=ls.get_color() )
        # off = p7/p6 - 3./2.
        # ls, =ax1.plot(t6, off, label='p7/p6')
        # ax1.axhline(y=1.519769 - 3./2., ls=":",color=ls.get_color() )

        ax1.set_ylabel(field+'  ')
        # ax1.set_yscale('log')
        # ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':10})

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_compare_period_ratio_and_observed_offset" +".png")
        plt.savefig(output_figure_filename)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_7:

        fig = plt.figure(figsize=(6, 6))
        ligne = 1
        colonne = 1
        i = 0
        #-----------------------------------------------------------------------------------

        # for key in planets_keys:
        p1 = planets_computed_data['1']['planet_rotation_period']
        p2 = planets_computed_data['2']['planet_rotation_period']
        p3 = planets_computed_data['3']['planet_rotation_period']
        p4 = planets_computed_data['4']['planet_rotation_period']
        p5 = planets_computed_data['5']['planet_rotation_period']
        p6 = planets_computed_data['6']['planet_rotation_period']
        p7 = planets_computed_data['7']['planet_rotation_period']

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = ' compare_period_ratio_and_theory '

        planet_data = planets_data[key]
        t1 = planets_data['1']['current_time']
        t2 = planets_data['2']['current_time']
        t3 = planets_data['3']['current_time']
        t4 = planets_data['4']['current_time']
        t5 = planets_data['5']['current_time']
        t6 = planets_data['6']['current_time']
        t7 = planets_data['7']['current_time']
 
        ls, =ax1.plot(t1, p2/p1, label='p2/p1')
        ax1.axhline(y=8./5., ls=":", color=ls.get_color(), label='oberved')
        ls, =ax1.plot(t2, p3/p2, label='p3/p2')
        ax1.axhline(y=5./3., ls=":", color=ls.get_color() )
        ls, =ax1.plot(t3, p4/p3, label='p4/p5')
        ax1.axhline(y=3./2., ls=":", color=ls.get_color() )
        ls, =ax1.plot(t4, p5/p4, label='p5/p4')
        ax1.axhline(y=3./2., ls=":", color=ls.get_color() )
        ls, =ax1.plot(t5, p6/p5, label='p6/p5')
        ax1.axhline(y=4./3., ls=":",color=ls.get_color() )
        ls, =ax1.plot(t6, p7/p6, label='p7/p6')
        ax1.axhline(y=3./2., ls=":",color=ls.get_color() )

        ax1.set_ylabel(field+'  ')
        # ax1.set_yscale('log')
        # ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':10})

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_compare_period_ratio_and_theory" +".png")
        plt.savefig(output_figure_filename)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_70:

        fig = plt.figure(figsize=(6, 6))
        ligne = 1
        colonne = 1
        i = 0
        #-----------------------------------------------------------------------------------

        # for key in planets_keys:
        p1 = planets_computed_data['1']['planet_rotation_period'][-1:]
        p2 = planets_computed_data['2']['planet_rotation_period'][-1:]
        p3 = planets_computed_data['3']['planet_rotation_period'][-1:]
        p4 = planets_computed_data['4']['planet_rotation_period'][-1:]
        p5 = planets_computed_data['5']['planet_rotation_period'][-1:]
        p6 = planets_computed_data['6']['planet_rotation_period'][-1:]
        p7 = planets_computed_data['7']['planet_rotation_period'][-1:]

        mean_p1 = np.mean(planets_computed_data['1']['planet_rotation_period'])
        mean_p2 = np.mean(planets_computed_data['2']['planet_rotation_period'])
        mean_p3 = np.mean(planets_computed_data['3']['planet_rotation_period'])
        mean_p4 = np.mean(planets_computed_data['4']['planet_rotation_period'])
        mean_p5 = np.mean(planets_computed_data['5']['planet_rotation_period'])
        mean_p6 = np.mean(planets_computed_data['6']['planet_rotation_period'])
        mean_p7 = np.mean(planets_computed_data['7']['planet_rotation_period'])

        # mean_p1 = planets_computed_data['1']['planet_rotation_period'][0]
        # mean_p2 = planets_computed_data['2']['planet_rotation_period'][0]
        # mean_p3 = planets_computed_data['3']['planet_rotation_period'][0]
        # mean_p4 = planets_computed_data['4']['planet_rotation_period'][0]
        # mean_p5 = planets_computed_data['5']['planet_rotation_period'][0]
        # mean_p6 = planets_computed_data['6']['planet_rotation_period'][0]
        # mean_p7 = planets_computed_data['7']['planet_rotation_period'][0]
        # print('TEST":::"', p1)

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = ' compare period ratio offsets '

        planet_data = planets_data[key]
        t1 = planets_data['1']['current_time'][-1:]
        t2 = planets_data['2']['current_time'][-1:]
        t3 = planets_data['3']['current_time'][-1:]
        t4 = planets_data['4']['current_time'][-1:]
        t5 = planets_data['5']['current_time'][-1:]
        t6 = planets_data['6']['current_time'][-1:]
        t7 = planets_data['7']['current_time'][-1:]

        # off = p2/p1 - (3./2.)
        # mean_off = mean_p2/mean_p1 - (3./2.)
        # ls, =ax1.plot(1, off,  'x',markersize=20, label='Last value')
        # ax1.plot(1, mean_off, '+',markersize=20,label='First value', color = ls.get_color())
        # ax1.plot(1, 1.603055 - (3./2.),  '.',markersize=20,label='Observed',  color = ls.get_color())

        off = p3/p2 - (5./3.)
        mean_off = mean_p3/mean_p2 - (5./3.)
        ls, =ax1.plot(2, off,  'x',markersize=20,label='Last value')
        ax1.plot(2, mean_off, '+',markersize=20, color = ls.get_color(),label='Mean value')
        ax1.plot(2, 1.671893 - (5./3.),  '.',markersize=20, color = ls.get_color(), label='Obs value')

        off = p4/p3 - (3./2.)
        mean_off = mean_p4/mean_p3 - (3./2.)
        ls, =ax1.plot(3, off, 'x',markersize=20,)
        ax1.plot(3, mean_off, '+',markersize=20, color = ls.get_color())
        ax1.plot(3, 1.506714 - (3./2.), '.',markersize=20, color = ls.get_color())

        off = p5/p4 - (3./2.)
        mean_off = mean_p5/mean_p4 - (3./2.)
        ls, =ax1.plot(4, off, 'x',markersize=20,)
        ax1.plot(4, mean_off, '+',markersize=20, color = ls.get_color())
        ax1.plot(4, 1.509182 - (3./2.), '.',markersize=20, color = ls.get_color())

        off = p6/p5 - (4./3.)
        mean_off = mean_p6/mean_p5 - (4./3.)
        ls, =ax1.plot(5, off, 'x',markersize=20,)
        ax1.plot(5, mean_off, '+',markersize=20, color = ls.get_color())
        ax1.plot(5, 1.341558 - (4./3.), '.',markersize=20, color = ls.get_color())

        off = p7/p6 - (3./2.)
        mean_off = mean_p7/mean_p6 - (3./2.)
        ls, =ax1.plot(6, off, 'x',markersize=20,)
        ax1.plot(6, mean_off, '+',markersize=20, color = ls.get_color())
        ax1.plot(6, 1.519769 - (3./2.), '.',markersize=20, color = ls.get_color())

        # plt.xticks([1,2,3,4,5,6],label=['p2/p1','p3/p2', 'p4/p3', 'p5/p4', 'p6/p5', 'p7/p6'],fontsize=16, rotation=90,)
        plt.xticks([1,2,3,4,5,6],label=['p2/p1','p3/p2', 'p4/p3', 'p5/p4', 'p6/p5', 'p7/p6'],fontsize=10, rotation=0,)

        plt.grid()
        ax1.set_ylabel(field+'  ')
        # ax1.set_yscale('log')
        # ax1.set_xscale('log')
        # ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':10})

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_compare_period_ratio_offsets" +".png")
        plt.savefig(output_figure_filename)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_70:

        fig = plt.figure(figsize=(6, 6))
        ligne = 1
        colonne = 1
        i = 0
        #-----------------------------------------------------------------------------------

        # for key in planets_keys:
        p1 = planets_computed_data['1']['planet_rotation_period'][-1:]
        p2 = planets_computed_data['2']['planet_rotation_period'][-1:]
        p3 = planets_computed_data['3']['planet_rotation_period'][-1:]
        p4 = planets_computed_data['4']['planet_rotation_period'][-1:]
        p5 = planets_computed_data['5']['planet_rotation_period'][-1:]
        p6 = planets_computed_data['6']['planet_rotation_period'][-1:]
        p7 = planets_computed_data['7']['planet_rotation_period'][-1:]

        mean_p1 = np.mean(planets_computed_data['1']['planet_rotation_period'])
        mean_p2 = np.mean(planets_computed_data['2']['planet_rotation_period'])
        mean_p3 = np.mean(planets_computed_data['3']['planet_rotation_period'])
        mean_p4 = np.mean(planets_computed_data['4']['planet_rotation_period'])
        mean_p5 = np.mean(planets_computed_data['5']['planet_rotation_period'])
        mean_p6 = np.mean(planets_computed_data['6']['planet_rotation_period'])
        mean_p7 = np.mean(planets_computed_data['7']['planet_rotation_period'])

        # mean_p1 = planets_computed_data['1']['planet_rotation_period'][0]
        # mean_p2 = planets_computed_data['2']['planet_rotation_period'][0]
        # mean_p3 = planets_computed_data['3']['planet_rotation_period'][0]
        # mean_p4 = planets_computed_data['4']['planet_rotation_period'][0]
        # mean_p5 = planets_computed_data['5']['planet_rotation_period'][0]
        # mean_p6 = planets_computed_data['6']['planet_rotation_period'][0]
        # mean_p7 = planets_computed_data['7']['planet_rotation_period'][0]
        # print('TEST":::"', p1)

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = ' compare period ratio  '

        planet_data = planets_data[key]
        t1 = planets_data['1']['current_time'][-1:]
        t2 = planets_data['2']['current_time'][-1:]
        t3 = planets_data['3']['current_time'][-1:]
        t4 = planets_data['4']['current_time'][-1:]
        t5 = planets_data['5']['current_time'][-1:]
        t6 = planets_data['6']['current_time'][-1:]
        t7 = planets_data['7']['current_time'][-1:]

        # off = p2/p1 - (3./2.)
        # mean_off = mean_p2/mean_p1 - (3./2.)
        # ls, =ax1.plot(1, off,  'x',markersize=20, label='Last value')
        # ax1.plot(1, mean_off, '+',markersize=20,label='First value', color = ls.get_color())
        # ax1.plot(1, 1.603055 - (3./2.),  '.',markersize=20,label='Observed',  color = ls.get_color())

        # off = p3/p2 - (5./3.)
        mean_off = mean_p3/mean_p2 
        # ls, =ax1.plot(2, off,  'x',markersize=20,label='Last value')
        ls, = ax1.plot(2, mean_off, '+',markersize=20, label='Mean value')
        ax1.plot(2,  (5./3.),  '.',markersize=20, color = ls.get_color(), label='Obs value')

        # off = p4/p3 - (3./2.)
        mean_off = mean_p4/mean_p3
        # ls, =ax1.plot(3, off, 'x',markersize=20,)
        ls, = ax1.plot(3, mean_off, '+',markersize=20, )
        ax1.plot(3, (3./2.), '.',markersize=20, color = ls.get_color())

        # off = p5/p4 - (3./2.)
        mean_off = mean_p5/mean_p4
        # ls, =ax1.plot(4, off, 'x',markersize=20,)
        ls, = ax1.plot(4, mean_off, '+',markersize=20, )
        ax1.plot(4, (3./2.), '.',markersize=20, color = ls.get_color())

        # off = p6/p5 - (4./3.)
        mean_off = mean_p6/mean_p5 
        # ls, =ax1.plot(5, off, 'x',markersize=20,)
        ls, = ax1.plot(5, mean_off, '+',markersize=20,)
        ax1.plot(5, (4./3.), '.',markersize=20, color = ls.get_color())

        # off = p7/p6 - (3./2.)
        mean_off = mean_p7/mean_p6
        # ls, =ax1.plot(6, off, 'x',markersize=20,)
        ls , = ax1.plot(6, mean_off, '+',markersize=20, )
        ax1.plot(6, (3./2.) , '.',markersize=20, color = ls.get_color())

        # plt.xticks([1,2,3,4,5,6],label=['p2/p1','p3/p2', 'p4/p3', 'p5/p4', 'p6/p5', 'p7/p6'],fontsize=16, rotation=90,)
        plt.xticks([1,2,3,4,5,6],label=['p2/p1','p3/p2', 'p4/p3', 'p5/p4', 'p6/p5', 'p7/p6'],fontsize=10, rotation=0,)

        plt.grid()
        ax1.set_ylabel(field+'  ')
        # ax1.set_yscale('log')
        # ax1.set_xscale('log')
        # ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':10})

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_compare_period_ratio" +".png")
        plt.savefig(output_figure_filename)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_71: # compare period ratio offset with observed - Evolution 

        fig = plt.figure(figsize=(6, 10))
        ligne = 2
        colonne = 1
        i = 0
        #-----------------------------------------------------------------------------------

        # for key in planets_keys:
        p1 = planets_computed_data['1']['planet_rotation_period']
        p2 = planets_computed_data['2']['planet_rotation_period']
        p3 = planets_computed_data['3']['planet_rotation_period']
        p4 = planets_computed_data['4']['planet_rotation_period']
        p5 = planets_computed_data['5']['planet_rotation_period']
        p6 = planets_computed_data['6']['planet_rotation_period']
        p7 = planets_computed_data['7']['planet_rotation_period']

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = ' compare period ratio offsets p$_i$/p$_{i+1} - i/j  '

        planet_data = planets_data[key]
        t1 = planets_data['1']['current_time']
        t2 = planets_data['2']['current_time']
        t3 = planets_data['3']['current_time']
        t4 = planets_data['4']['current_time']
        t5 = planets_data['5']['current_time']
        t6 = planets_data['6']['current_time']
        t7 = planets_data['7']['current_time']
 
        off = p2/p1 -8./5.
        ls, =ax1.plot(t1,off , label='p2/p1 -8/5')
        ax1.axhline(y=1.603055 -8./5., ls=":", color=ls.get_color() )

        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        # ----------
        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = ' compare period ratio offsets p$_i$/p$_{i+1} - i/j  '
        
        off = p3/p2 - 5./3.
        ls, =ax1.plot(t2, off, label='p3/p2 -5/3')
        ax1.axhline(y=1.67189- 5./3., ls=":", color=ls.get_color() )
        off = p4/p3 - 3./2.
        ls, =ax1.plot(t3, off, label='p4/p5 -3/2')
        ax1.axhline(y=1.506714- 3./2., ls=":", color=ls.get_color() )
        off = p5/p4 - 3./2.
        ls, =ax1.plot(t4, off, label='p5/p4 -3/2')
        ax1.axhline(y=1.509182 - 3./2., ls=":", color=ls.get_color() )
        off = p6/p5 - 4./3.
        ls, =ax1.plot(t5, off, label='p6/p5 -4/3')
        ax1.axhline(y=1.341558 - 4./3., ls=":",color=ls.get_color() )
        off = p7/p6 - 3./2.
        ls, =ax1.plot(t6, off, label='p7/p6 -3/2')
        ax1.axhline(y=1.519769 - 3./2., ls=":",color=ls.get_color() )

        ax1.set_ylabel(field+'  ')
        # ax1.set_yscale('log')
        # ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_compare_period_ratio_offsets_drift" +".png")
        plt.savefig(output_figure_filename)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_72: # compare period ratio offset with observed - Evolution each planets

        fig = plt.figure(figsize=(20, 15))
        ligne = 2
        colonne = 3
        i = 0
        #-----------------------------------------------------------------------------------

        # for key in planets_keys:
        p1 = planets_computed_data['1']['planet_rotation_period']
        p2 = planets_computed_data['2']['planet_rotation_period']
        p3 = planets_computed_data['3']['planet_rotation_period']
        p4 = planets_computed_data['4']['planet_rotation_period']
        p5 = planets_computed_data['5']['planet_rotation_period']
        p6 = planets_computed_data['6']['planet_rotation_period']
        p7 = planets_computed_data['7']['planet_rotation_period']

        
        field = ' compare period ratio offsets p$_i$/p$_{i+1} - i/j  '

        planet_data = planets_data[key]
        t1 = planets_data['1']['current_time']
        t2 = planets_data['2']['current_time']
        t3 = planets_data['3']['current_time']
        t4 = planets_data['4']['current_time']
        t5 = planets_data['5']['current_time']
        t6 = planets_data['6']['current_time']
        t7 = planets_data['7']['current_time']
 
        # ---------- 1
        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        ax1.axhline(y=0.0, ls=":", color='black', label='No off' )
        off = p2/p1 -3./2.
        print(f">> off init 3:2 {off[0]} ")
        ls, =ax1.plot(t1,off , label='p2/p1 -3/2')
        ax1.axhline(y=1.603055 -3./2., ls=":", color=ls.get_color(), label='off from 3:2 Obs' )
        off = p2/p1 -8./5.
        print(f">> off init 8:5 {off[0]} ")
        ls, =ax1.plot(t1,off , label='p2/p1 -8/5')
        ax1.axhline(y=1.603055 -8./5., ls=":", color=ls.get_color(), label='off from 8:5 Obs' )
        #
        xx = t1[t1<1e3]
        yy = off[t1<1e3]
        a = plt.axes([0.13, 0.78, .16, .12],) # ( position[x,y], size[x,y] )
        plt.plot(xx, yy, color=ls.get_color())
        plt.axhline(y=1.603055 -8./5., ls=":", color=ls.get_color(), label='off from 8:5 Obs' )
        xx = t1[t1>1e3]
        yy = off[t1>1e3]
        a = plt.axes([0.13, 0.63, .16, .12],) # ( position[x,y], size[x,y] )
        plt.plot(xx, yy, color=ls.get_color())
        plt.axhline(y=1.603055 -8./5., ls=":", color=ls.get_color(), label='off from 8:5 Obs' )
        #
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        # ---------- 2
        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = ' compare period ratio offsets p$_{i}$/p$_{i+1}$ - i/j  '
        off = p3/p2 - 5./3.
        ls, =ax1.plot(t2, off, label='p3/p2 -5/3')
        ax1.axhline(y=1.67189- 5./3., ls=":", color=ls.get_color() )
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})
        #
        xx = t2[t2<1e3]
        yy = off[t2<1e3]
        a = plt.axes([0.45, 0.78, .16, .12],) # ( position[x,y], size[x,y] )
        plt.plot(xx, yy, color='blue')
        plt.axhline(y=1.67189- 5./3., ls=":", color=ls.get_color(), label='Obs' )
        xx = t2[t2>1e3]
        yy = off[t2>1e3]
        a = plt.axes([0.45, 0.63, .16, .12],) # ( position[x,y], size[x,y] )
        plt.plot(xx, yy, color='blue')
        plt.axhline(y=1.67189- 5./3., ls=":", color=ls.get_color(), label='Obs' )

        # ---------- 3
        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        off = p4/p3 - 3./2.
        ls, =ax1.plot(t3, off, label='p4/p5 -3/2')
        ax1.axhline(y=1.506714- 3./2., ls=":", color=ls.get_color() )
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})
        #
        xx = t3[t3<1.5e4]
        yy = off[t3<1.5e4]
        a = plt.axes([0.8, 0.78, .16, .12],) # ( position[x,y], size[x,y] )
        plt.plot(xx, yy, color='blue')
        plt.axhline(y=1.506714- 3./2., ls=":", color=ls.get_color(), label='Obs' )
        xx = t3[t3>1.5e4]
        yy = off[t3>1.5e4]
        a = plt.axes([0.8, 0.63, .16, .12],) # ( position[x,y], size[x,y] )
        plt.plot(xx, yy, color='blue')
        plt.axhline(y=1.506714- 3./2., ls=":", color=ls.get_color(), label='Obs' )

        # ---------- 4
        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        off = p5/p4 - 3./2.
        ls, =ax1.plot(t4, off, label='p5/p4 -3/2')
        ax1.axhline(y=1.509182 - 3./2., ls=":", color=ls.get_color() )
        #
        # xx = t1[t1<1e3]
        # yy = off[t1<1e3]
        # a = plt.axes([0.13, 0.45, .16, .12],) # ( position[x,y], size[x,y] )
        # plt.plot(xx, yy, color=ls.get_color())
        # plt.axhline(y=1.509182 - 3./2., ls=":", color=ls.get_color(), label='off from 8:5 Obs' ) 
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        # ----------
        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        off = p6/p5 - 4./3.
        ls, =ax1.plot(t5, off, label='p6/p5 -4/3')
        ax1.axhline(y=1.341558 - 4./3., ls=":",color=ls.get_color() )
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        # ----------
        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        off = p7/p6 - 3./2.
        ls, =ax1.plot(t6, off, label='p7/p6 -3/2')
        ax1.axhline(y=1.519769 - 3./2., ls=":",color=ls.get_color() )
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        ax1.set_ylabel(field+'  ')
        # ax1.set_yscale('log')
        # ax1.set_xscale('log')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})
        

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_compare_period_ratio_offsets_drift_eachplanets" +".png")
        plt.savefig(output_figure_filename)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_73: # Fig-1 Charalambous et al. 2023 | Departure from 2-body MMRs

        def y_fct(x, k1, k2, k3):
            return -(k2+k1/x)/k3

        fig = plt.figure(figsize=(7, 7))
        ligne = 1
        colonne = 1
        i = 0
        #-----------------------------------------------------------------------------------

        # for key in planets_keys:
        p1 = planets_computed_data['1']['planet_rotation_period']
        p2 = planets_computed_data['2']['planet_rotation_period']
        p3 = planets_computed_data['3']['planet_rotation_period']
        p4 = planets_computed_data['4']['planet_rotation_period']
        p5 = planets_computed_data['5']['planet_rotation_period']
        p6 = planets_computed_data['6']['planet_rotation_period']
        p7 = planets_computed_data['7']['planet_rotation_period']

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = '  '

        # planet_data = planets_data[key]
        # t1 = planets_data['1']['current_time']
        # t2 = planets_data['2']['current_time']
        # t3 = planets_data['3']['current_time']
        # t4 = planets_data['4']['current_time']
        # t5 = planets_data['5']['current_time']
        # t6 = planets_data['6']['current_time']
        # t7 = planets_data['7']['current_time']
 
        x = p2/p1
        y = p3/p2
        ls, =ax1.plot(x,y, '+')

        x = p3/p2
        y = p4/p3
        ls, =ax1.plot(x,y, '+')

        x = p4/p3
        y = p5/p4
        ls, =ax1.plot(x,y, '+')

        x = p5/p4
        y = p6/p5
        ls, =ax1.plot(x,y, '+')

        x = p6/p5
        y = p7/p6
        ls, =ax1.plot(x,y, '+')

        x = np.arange(1.3,1.5,0.1)
        y = y_fct(x,3,-6,3)
        ax1.plot(x,y, '--')

        ax1.axhline(y=4/3, ls=':', color='black')
        ax1.axhline(y=3/2, ls=':', color='black')
        ax1.axhline(y=5/3, ls=':', color='black')
        ax1.axhline(y=2, ls=':', color='black')
        ax1.axvline(x=4/3, ls=':', color='black')
        ax1.axvline(x=3/2, ls=':', color='black')
        ax1.axvline(x=5/3, ls=':', color='black')
        ax1.axvline(x=2, ls=':', color='black')

        ax1.set_ylabel('n$_2$/n$_3$')
        ax1.set_xlabel('n$_1$/n$_2$')
        ax1.legend(loc=0, prop={'size':8})


        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_departure_from_2MMR" +".png")
        plt.savefig(output_figure_filename)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_8:
        fig = plt.figure(figsize=(16, 16))
        ligne = 2
        colonne = 3
        i = 0
        #-----------------------------------------------------------------------------------
        p_1 = planets_computed_data['1']['planet_rotation_period']
        t_1 = planets_data['1']['current_time']
        for key in planets_keys[1:]:

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = ' planet period ratio chain '

            planet_data = planets_data[key]
            t_2 = planet_data['current_time']
            p_2 = planets_computed_data[key]['planet_rotation_period']
            
            # if key == '1':
            #     ax1.axhline(y=1.510826, ls=":", color='black' )
            # if key == '2':
            #     ax1.axhline(y=2.421937, ls=":", color='black' )
            # if key == '3':
            #     ax1.axhline(y=4.049219, ls=":", color='black' )
            # if key == '4':
            #     ax1.axhline(y=6.101013, ls=":", color='black' )
            # if key == '5':
            #     ax1.axhline(y=9.20754, ls=":", color='black' )
            # if key == '6':
            #     ax1.axhline(y=12.352446, ls=":", color='black' )
            # if key == '7':
            #     ax1.axhline(y=18.772866, ls=":", color='black' )
            ax1.plot(t_2, p_2/p_1, label='mean motion', alpha=0.5)
            if key == '2':
                ax1.axhline(y=1.603055, ls=":", color='black' )
            if key == '3':
                ax1.axhline(y=1.671893, ls=":", color='black' )
            if key == '4':
                ax1.axhline(y=1.506714, ls=":", color='black' )
            if key == '5':
                ax1.axhline(y=1.509182, ls=":", color='black' )
            if key == '6':
                ax1.axhline(y=1.341558, ls=":", color='black' )
            if key == '7':
                ax1.axhline(y=1.519769, ls=":", color='black' )
            p_1 = p_2
            t_1 = t_2
        
            ax1.set_ylabel(field+' Planet: '+ str(key))
            ax1.set_yscale('log')
            # ax1.set_yscale('log')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_planets_period_ratio_chain" +".png")
        plt.savefig(output_figure_filename)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_9: # Resonance 3 bodies

        fig = plt.figure(figsize=(8, 16))
        fig.suptitle('3-body resonant angles')
        ligne = 5
        colonne = 1
        i = 0
        #-----------------------------------------------------------------------------------
        
        # t1 = planets_data['1']['current_time']
        # t2 = planets_data['2']['current_time']
        # t3 = planets_data['3']['current_time']
        # t4 = planets_data['4']['current_time']
        # t5 = planets_data['5']['current_time']
        # t6 = planets_data['6']['current_time']
        t7 = planets_data['7']['current_time']

        l_1 = planets_computed_data['1']['mean_longitude'] 
        l_2 = planets_computed_data['2']['mean_longitude'] 
        l_3 = planets_computed_data['3']['mean_longitude'] 
        l_4 = planets_computed_data['4']['mean_longitude'] 
        l_5 = planets_computed_data['5']['mean_longitude'] 
        l_6 = planets_computed_data['6']['mean_longitude'] 
        l_7 = planets_computed_data['7']['mean_longitude'] 

        # omega_1 =  planets_computed_data['1']['check_longitude_of_perihelion']
        # omega_2 =  planets_computed_data['2']['check_longitude_of_perihelion']
        # omega_3 =  planets_computed_data['3']['check_longitude_of_perihelion']
        # omega_4 =  planets_computed_data['4']['check_longitude_of_perihelion']
        # omega_5 =  planets_computed_data['5']['check_longitude_of_perihelion']
        # omega_6 =  planets_computed_data['6']['check_longitude_of_perihelion']
        # omega_7 =  planets_computed_data['7']['check_longitude_of_perihelion']
        omega_1 =  planets_data['1']['longitude_of_perihelion']
        omega_2 =  planets_data['2']['longitude_of_perihelion']
        omega_3 =  planets_data['3']['longitude_of_perihelion']
        omega_4 =  planets_data['4']['longitude_of_perihelion']
        omega_5 =  planets_data['5']['longitude_of_perihelion']
        omega_6 =  planets_data['6']['longitude_of_perihelion']
        omega_7 =  planets_data['7']['longitude_of_perihelion']

        theta_123 = 10.0*l_1 -25.0*l_2 +15.0*l_3
        theta_234 = 3.0*l_2 -9.0*l_3 +6.0*l_4
        theta_345 = 2.0*l_3 -5.0*l_4 +3.0*l_5
        theta_456 = 2.0*l_4 -6.0*l_5 +4.0*l_6
        theta_567 = 3.0*l_5 -6.0*l_6 +3.0*l_7

        theta_123 = theta_123%(2.0*np.pi)
        theta_234 = theta_234%(2.0*np.pi)
        theta_345 = theta_345%(2.0*np.pi)
        theta_456 = theta_456%(2.0*np.pi)
        theta_567 = theta_567%(2.0*np.pi)

        theta_123 = (theta_123)*180./np.pi
        theta_234 = (theta_234)*180./np.pi
        theta_345 = (theta_345)*180./np.pi
        theta_456 = (theta_456)*180./np.pi
        theta_567 = (theta_567)*180./np.pi

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = '$\\theta$ $_{10:25:15}$ (b+c+d)'
        ax1.plot(t7, theta_123, '+', label=f'angle t=0 {round(theta_123[0],3)}')
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = '$\\theta$ $_{3:9:6}$ (c+d+e)'
        ax1.plot(t7, theta_234, '+', label=f'angle t=0 {round(theta_234[0],3)}')
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = '$\\theta$ $_{2:5:3}$ (d+e+f)'
        ax1.plot(t7, theta_345, '+', label=f'angle t=0: {round(theta_345[0],3)}')
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = '$\\theta$ $_{2:6:4}$ (e+f+g)'
        ax1.plot(t7, theta_456, '+', label=f'angle t=0:: {round(theta_456[0],3)}')
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        # field = 'Res Angle f g h '
        field = '$\\theta$ $_{3:6:3}$ (f+g+h)'
        ax1.plot(t7, theta_567, '+', label= f'angle t=0: {round(theta_567[0],3)}')
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})

        #-----------------------------------------------------------------------------------
        plt.tight_layout()
        
        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_reso_angles_3body" +".png")
        plt.savefig(output_figure_filename)
        # plt.show()
        plt.close()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_10: # Resonance 2 bodies

        fig = plt.figure(figsize=(8, 16))
        fig.suptitle('2-body resonant angles')
        ligne = 6
        colonne = 1
        i = 0
        labels = [-180, -90, 0, 90, 180]
        #-----------------------------------------------------------------------------------

        planet_data = planets_data[key]
        # t1 = planets_data['1']['current_time']
        # t2 = planets_data['2']['current_time']
        # t3 = planets_data['3']['current_time']
        # t4 = planets_data['4']['current_time']
        # t5 = planets_data['5']['current_time']
        # t6 = planets_data['6']['current_time']
        t7 = planets_data['7']['current_time']

        # l_1 = planets_computed_data['1']['mean_longitude']
        # l_2 = planets_computed_data['2']['mean_longitude']
        # l_3 = planets_computed_data['3']['mean_longitude']
        # l_4 = planets_computed_data['4']['mean_longitude']
        # l_5 = planets_computed_data['5']['mean_longitude']
        # l_6 = planets_computed_data['6']['mean_longitude']
        # l_7 = planets_computed_data['7']['mean_longitude']

        l_b = planets_computed_data['1']['mean_longitude']
        l_c = planets_computed_data['2']['mean_longitude']
        l_d = planets_computed_data['3']['mean_longitude']
        l_e = planets_computed_data['4']['mean_longitude']
        l_f = planets_computed_data['5']['mean_longitude']
        l_g = planets_computed_data['6']['mean_longitude']
        l_h = planets_computed_data['7']['mean_longitude']


        # omega_b =  planets_computed_data['1']['check_longitude_of_perihelion']
        # omega_c =  planets_computed_data['2']['check_longitude_of_perihelion']
        # omega_d =  planets_computed_data['3']['check_longitude_of_perihelion']
        # omega_e =  planets_computed_data['4']['check_longitude_of_perihelion']
        # omega_f =  planets_computed_data['5']['check_longitude_of_perihelion']
        # omega_g =  planets_computed_data['6']['check_longitude_of_perihelion']
        # omega_h =  planets_computed_data['7']['check_longitude_of_perihelion']

        omega_b =  planets_data['1']['longitude_of_perihelion']
        omega_c =  planets_data['2']['longitude_of_perihelion']
        omega_d =  planets_data['3']['longitude_of_perihelion']
        omega_e =  planets_data['4']['longitude_of_perihelion']
        omega_f =  planets_data['5']['longitude_of_perihelion']
        omega_g =  planets_data['6']['longitude_of_perihelion']
        omega_h =  planets_data['7']['longitude_of_perihelion']

        # 8:5 | 4 angles b c
        theta1_85_bc = 8.0*l_c -5.0*l_b -3.0*omega_b
        theta2_85_bc = 8.0*l_c -5.0*l_b -3.0*omega_c
        theta3_85_bc = 8.0*l_c -5.0*l_b -2.0*omega_b -omega_c
        theta4_85_bc = 8.0*l_c -5.0*l_b -omega_b -2.0*omega_c
        theta1_85_bc = theta1_85_bc%(2.0*np.pi)
        theta2_85_bc = theta2_85_bc%(2.0*np.pi)
        theta3_85_bc = theta3_85_bc%(2.0*np.pi)
        theta4_85_bc = theta4_85_bc%(2.0*np.pi)
        theta1_85_bc = (np.pi -theta1_85_bc) *180./np.pi
        theta2_85_bc = (np.pi -theta2_85_bc) *180./np.pi
        theta3_85_bc = (np.pi -theta3_85_bc) *180./np.pi
        theta4_85_bc = (np.pi -theta4_85_bc) *180./np.pi

        # 5:3 | 3 angles c d 
        theta1_53_cd = 5.0*l_d -3.0*l_c -2.0*omega_c
        theta2_53_cd = 5.0*l_d -3.0*l_c -2.0*omega_d
        theta3_53_cd = 5.0*l_d -3.0*l_c -omega_c -omega_d
        theta1_53_cd = (theta1_53_cd%(2.0*np.pi) - np.pi) %(2.0*np.pi)
        theta2_53_cd = (theta2_53_cd%(2.0*np.pi) - np.pi) %(2.0*np.pi)
        theta3_53_cd = (theta3_53_cd%(2.0*np.pi) - np.pi) %(2.0*np.pi)
        theta1_53_cd = (np.pi-theta1_53_cd)*180./np.pi
        theta2_53_cd = (np.pi-theta2_53_cd)*180./np.pi
        theta3_53_cd = (np.pi-theta3_53_cd)*180./np.pi

        # 3:2 | 2 angles d e
        theta1_32_de = 3.0*l_e -2.0*l_d -omega_d
        theta2_32_de = 3.0*l_e -2.0*l_d -omega_e
        theta1_32_de = (theta1_32_de%(2.0*np.pi) - np.pi )%(2.0*np.pi)
        theta2_32_de = (theta2_32_de%(2.0*np.pi) - np.pi )%(2.0*np.pi)
        theta1_32_de = (np.pi-theta1_32_de)*180./np.pi
        theta2_32_de = (np.pi-theta2_32_de)*180./np.pi

        # 3:2 | 2 angles e f
        theta1_32_ef = 3.0*l_f -2.0*l_e -omega_e
        theta2_32_ef = 3.0*l_f -2.0*l_e -omega_f
        theta1_32_ef = (theta1_32_ef%(2.0*np.pi) -np.pi ) %(2.0*np.pi)
        theta2_32_ef = (theta2_32_ef%(2.0*np.pi) -np.pi ) %(2.0*np.pi)
        theta1_32_ef = (np.pi-theta1_32_ef)*180./np.pi
        theta2_32_ef = (np.pi-theta2_32_ef)*180./np.pi

        # 4:3 | 2 angles f g 
        theta1_43_fg = 4.0*l_g -3.0*l_f -omega_f
        theta2_43_fg = 4.0*l_g -3.0*l_f -omega_g
        theta1_43_fg = (theta1_43_fg%(2.0*np.pi) -np.pi ) %(2.0*np.pi)
        theta2_43_fg = (theta2_43_fg%(2.0*np.pi) -np.pi ) %(2.0*np.pi)
        theta1_43_fg = (np.pi-theta1_43_fg)*180./np.pi
        theta2_43_fg = (np.pi-theta2_43_fg)*180./np.pi

        #3:2 | 2 angles g h
        theta1_32_gh = 3.0*l_h -2.0*l_g -omega_g
        theta2_32_gh = 3.0*l_h -2.0*l_g -omega_h
        theta1_32_gh = (theta1_32_gh%(2.0*np.pi) -np.pi ) %(2.0*np.pi)
        theta2_32_gh = (theta2_32_gh%(2.0*np.pi) -np.pi ) %(2.0*np.pi)
        theta1_32_gh = (np.pi-theta1_32_gh)*180./np.pi
        theta2_32_gh = (np.pi-theta2_32_gh)*180./np.pi
        
        # print("\n --- Check resonant angles --- ")
        # print(f" b-c 8:5 | {theta1_85_bc[0]*(180./np.pi)} | {theta2_85_bc[0]*(180./np.pi)} | {theta3_85_bc[0]*(180./np.pi)}")
        # print(f" c-d 5:3 | {theta1_53_cd[0]*(180./np.pi)} | {theta2_53_cd[0]*(180./np.pi)} | {theta3_53_cd[0]*(180./np.pi)}")
        # print(f" d-e 3:2 | {theta1_32_de[0]*(180./np.pi)} | {theta2_32_de[0]*(180./np.pi)} ")
        # print(f" e-f 3:2 | {(np.pi-theta1_32_ef[0])*(180./np.pi)} | {(np.pi-theta2_32_ef[0])*(180./np.pi)} ")
        # print(f" f-g 4:3 | {(np.pi-theta1_43_fg[0])*(180./np.pi)} | {(np.pi-theta2_43_fg[0])*(180./np.pi)} ")
        # print(f" g-f 3:2 | {theta1_32_gh[0]*(180./np.pi)} | {theta2_32_gh[0]*(180./np.pi)} ")

        print("\n --- Check resonant angles --- ")
        print(f" b-c 8:5 | {theta1_85_bc[0]} | {theta2_85_bc[0]} | {theta3_85_bc[0]}")
        print(f" c-d 5:3 | {theta1_53_cd[0]} | {theta2_53_cd[0]} | {theta3_53_cd[0]}")
        print(f" d-e 3:2 | {theta1_32_de[0]} | {theta2_32_de[0]} ")
        print(f" e-f 3:2 | {(theta1_32_ef[0])} | {(theta2_32_ef[0])} ")
        print(f" f-g 4:3 | {(theta1_43_fg[0])} | {(theta2_43_fg[0])} ")
        print(f" g-f 3:2 | {theta1_32_gh[0]} | {theta2_32_gh[0]} ")

        # --- plot d
        
        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = '$\\theta$ $_{8:5}$ (b+c)'
        ax1.plot(t7, theta1_85_bc, '+', label=f'angle 1 t=0: {round(theta1_85_bc[0],3)}', alpha=1.0)
        ax1.plot(t7, theta2_85_bc, '+', label=f'angle 2 t=0: {round(theta2_85_bc[0],3)}', alpha=0.5)
        ax1.plot(t7, theta3_85_bc, '+', label=f'angle 3 t=0: {round(theta3_85_bc[0],3)}', alpha=0.4)
        ax1.plot(t7, theta4_85_bc, '+', label=f'angle 4 t=0: {round(theta4_85_bc[0],3)}', alpha=0.3)
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})
        ax1.set_yticks(labels)
        ax1.set_yticklabels(labels)
        # i=i+1
        # ax1 = fig.add_subplot(ligne,colonne,i,)
        # field = 'Res Angle 8:5 | b c angle 2'
        # ax1.plot(t7, theta2_85_bc, '+', label=f'angle 2 init: {round(theta2_85_bc[0],3)}')
        # ax1.set_ylabel(field+'  ')
        # ax1.set_xlabel('years')
        # ax1.legend(loc=0, prop={'size':8})
        # i=i+1
        # ax1 = fig.add_subplot(ligne,colonne,i,)
        # field = 'Res Angle 8:5 | b c angle 3'
        # ax1.plot(t7, theta3_85_bc, '+', label=f'angle 3 init: {round(theta3_85_bc[0],3)}')
        # ax1.set_ylabel(field+'  ')
        # ax1.set_xlabel('years')
        # ax1.legend(loc=0, prop={'size':8})
        # i=i+1
        # ax1 = fig.add_subplot(ligne,colonne,i,)
        # field = 'Res Angle 8:5 | b c angle 4'
        # ax1.plot(t7, theta4_85_bc, '+', label=f'angle 4 init: {round(theta4_85_bc[0],3)}')
        # ax1.set_ylabel(field+'  ')
        # ax1.set_xlabel('years')
        # ax1.legend(loc=0, prop={'size':8})

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = '$\\theta$ $_{5:3}$ (c+d)'
        ax1.plot(t7, theta1_53_cd, '+', label=f'angle 1 t=0: {round(theta1_53_cd[0],3)}')
        ax1.plot(t7, theta2_53_cd, '+', label=f'angle 2 t=0: {round(theta2_53_cd[0],3)}', alpha=0.5)
        ax1.plot(t7, theta3_53_cd, '+', label=f'angle 3 t=0: {round(theta3_53_cd[0],3)}', alpha=0.3)
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})
        ax1.set_yticks(labels)
        ax1.set_yticklabels(labels)
        # i=i+1
        # ax1 = fig.add_subplot(ligne,colonne,i,)
        # field = 'Res Angle 5:3 | c d angle 2'
        # ax1.plot(t7, theta2_53_cd, '+', label=f'angle 2 init: {round(theta2_53_cd[0],3)}', alpha=1.)
        # ax1.set_ylabel(field+'  ')
        # ax1.set_xlabel('years')
        # ax1.legend(loc=0, prop={'size':8})
        # i=i+1
        # ax1 = fig.add_subplot(ligne,colonne,i,)
        # field = 'Res Angle 5:3 | c d angle 3'
        # ax1.plot(t7, theta3_53_cd, '+', label=f'angle 3 init: {round(theta3_53_cd[0],3)}', alpha=1.)
        # ax1.set_ylabel(field+'  ')
        # ax1.set_xlabel('years')
        # ax1.legend(loc=0, prop={'size':8})

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = '$\\theta$ $_{3:2}$ (d+e)'
        ax1.plot(t7, theta1_32_de, '+', label=f'angle 1 t=0 {round(theta1_32_de[0],3)}')
        ax1.plot(t7, theta2_32_de, '+', label=f'angle 2 t=0 {round(theta2_32_de[0],3)}', alpha=0.4)
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})
        ax1.set_yticks(labels)
        ax1.set_yticklabels(labels)

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = '$\\theta$ $_{3:2}$ (e+f)'
        ax1.plot(t7, theta1_32_ef, '+', label=f'angle 1 t=0 {round(theta1_32_ef[0],3)}')
        ax1.plot(t7, theta2_32_ef, '+', label=f'angle 2 t=0 {round(theta2_32_ef[0],3)}', alpha=0.4)
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})
        ax1.set_yticks(labels)
        ax1.set_yticklabels(labels)

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = '$\\theta$ $_{4:3}$ (f+g)'
        ax1.plot(t7, theta1_43_fg, '+', label=f'angle 1 t=0 {round(theta1_43_fg[0],3)}')
        ax1.plot(t7, theta2_43_fg, '+', label=f'angle 2 t=0 {round(theta2_43_fg[0],3)}', alpha=0.4)
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})
        ax1.set_yticks(labels)
        ax1.set_yticklabels(labels)

        i=i+1
        ax1 = fig.add_subplot(ligne,colonne,i,)
        field = '$\\theta$ $_{3:2}$ (g+h)'
        ax1.plot(t7, theta1_32_gh, '+', label=f'angle 1 t=0 {round(theta1_32_gh[0],3)}')
        ax1.plot(t7, theta2_32_gh, '+', label=f'angle 2 t=0 {round(theta2_32_gh[0],3)}', alpha=0.4)
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')
        ax1.legend(loc=0, prop={'size':8})
        ax1.set_yticks(labels)
        ax1.set_yticklabels(labels)

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_reso_angles_2body" +".png")
        plt.savefig(output_figure_filename)
        # plt.show()
        plt.close()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_101: # Resonance 2 bodies 2 planets

        fig = plt.figure(figsize=(13, 8))
        ligne = 3
        colonne = 2
        plot = 0
        #-----------------------------------------------------------------------------------
        p1 = planets_computed_data['1']['planet_rotation_period']
        p2 = planets_computed_data['2']['planet_rotation_period']

        planet_data = planets_data[key]
        t1 = planets_data['1']['current_time']

        l_1 = planets_computed_data['1']['mean_longitude']
        l_2 = planets_computed_data['2']['mean_longitude']

        omega_1 =  planets_data['1']['longitude_of_perihelion']
        omega_2 =  planets_data['2']['longitude_of_perihelion']

        # 2:1 | 1 angles 
        theta1_21_p1p2 = 2.0*l_2 -l_1 -omega_1
        theta1_21_p1p2 = theta1_21_p1p2%(2.0*np.pi)

        # 3:2 | 2 angles 
        theta1_32_p1p2 = 3.0*l_2 -2.0*l_1 -omega_1
        theta2_32_p1p2 = 3.0*l_2 -2.0*l_1 -omega_2
        theta1_32_p1p2 = theta1_32_p1p2%(2.0*np.pi)
        theta2_32_p1p2 = theta2_32_p1p2%(2.0*np.pi)

        # 4:3 | 2 angles
        theta1_43_p1p2 = 4.0*l_2 -3.0*l_1 -omega_1
        theta2_43_p1p2 = 4.0*l_2 -3.0*l_1 -omega_2
        theta1_43_p1p2 = theta1_43_p1p2%(2.0*np.pi)
        theta2_43_p1p2 = theta2_43_p1p2%(2.0*np.pi)

        # 8:5 | 4 angles b c
        theta1_85_p1p2 = 8.0*l_2 -5.0*l_1 -3.0*omega_1
        theta2_85_p1p2 = 8.0*l_2 -5.0*l_1 -3.0*omega_2
        theta3_85_p1p2 = 8.0*l_2 -5.0*l_1 -2.0*omega_1 -omega_2
        theta4_85_p1p2 = 8.0*l_2 -5.0*l_1 -omega_1 -2.0*omega_2
        theta1_85_p1p2 = theta1_85_p1p2%(2.0*np.pi)
        theta2_85_p1p2 = theta2_85_p1p2%(2.0*np.pi)
        theta3_85_p1p2 = theta3_85_p1p2%(2.0*np.pi)
        theta4_85_p1p2 = theta4_85_p1p2%(2.0*np.pi)

        # 5:3 | 3 angles c d 
        theta1_53_p1p2 = 5.0*l_2 -3.0*l_1 -2.0*omega_1
        theta2_53_p1p2 = 5.0*l_2 -3.0*l_1 -2.0*omega_2
        theta3_53_p1p2 = 5.0*l_2 -3.0*l_1 -omega_1 -omega_2
        theta1_53_p1p2 = theta1_53_p1p2%(2.0*np.pi)
        theta2_53_p1p2 = theta2_53_p1p2%(2.0*np.pi)
        theta3_53_p1p2 = theta3_53_p1p2%(2.0*np.pi)


        # --- plot 2:1
        plot+=1
        ax1 = fig.add_subplot(ligne,colonne,plot,)
        field = 'Res Angle 2:1 | planet 1-2 '
        ax1.plot(t1, theta1_21_p1p2, 'x', label='angle 1')
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')

        # --- plot 3:2
        plot+=1
        ax1 = fig.add_subplot(ligne,colonne,plot,)
        field = 'Res Angle 3:2 | planet 1-2 '
        ax1.plot(t1, theta1_32_p1p2, 'x', label='angle 1')
        ax1.plot(t1, theta2_32_p1p2, '+', label='angle 2', alpha=.7)
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')

        # --- plot 4:3
        plot+=1
        ax1 = fig.add_subplot(ligne,colonne,plot,)
        field = 'Res Angle 4:3 | planet 1-2 '
        ax1.plot(t1, theta1_43_p1p2, 'x', label='angle 1')
        ax1.plot(t1, theta2_43_p1p2, '+', label='angle 2', alpha=.7)
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')

        # --- plot 8:5
        plot+=1
        ax1 = fig.add_subplot(ligne,colonne,plot,)
        field = 'Res Angle 8:5 | planet 1-2 '
        ax1.plot(t1, theta1_85_p1p2, 'x', label='angle 1')
        ax1.plot(t1, theta2_85_p1p2, '+', label='angle 2', alpha=.7)
        ax1.plot(t1, theta3_85_p1p2, 'x', label='angle 2', alpha=.3)
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')

        # --- plot 5:3
        plot+=1
        ax1 = fig.add_subplot(ligne,colonne,plot,)
        field = 'Res Angle 5:3 | planet 1-2 '
        ax1.plot(t1, theta1_53_p1p2, 'x', label='angle 1')
        ax1.plot(t1, theta2_53_p1p2, '+', label='angle 2', alpha=.7)
        ax1.plot(t1, theta3_53_p1p2, 'x', label='angle 2', alpha=.3)
        ax1.set_ylabel(field+'  ')
        ax1.set_xlabel('years')

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_2planet_res_angles_2body" +".png")
        plt.savefig(output_figure_filename)
        plt.close()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_102: # Look for Resonance 2 bodies by pairs

        for key1, key2 in zip(planets_keys, planets_keys[1:]):

            p1 = planets_computed_data[key1]['planet_rotation_period']
            p2 = planets_computed_data[key2]['planet_rotation_period']

            planet_data = planets_data[key]
            t1 = planets_data[key1]['current_time']

            l_1 = planets_computed_data[key1]['mean_longitude']
            l_2 = planets_computed_data[key2]['mean_longitude']

            omega_1 =  planets_data[key1]['longitude_of_perihelion']
            omega_2 =  planets_data[key2]['longitude_of_perihelion']

            # 2:1 | 1 angles 
            theta1_21_p1p2 = 2.0*l_2 -l_1 -omega_1
            theta2_21_p1p2 = 2.0*l_2 -l_1 -omega_2
            theta1_21_p1p2 = theta1_21_p1p2%(2.0*np.pi)
            theta2_21_p1p2 = theta2_21_p1p2%(2.0*np.pi)
            theta1_21_p1p2 = (np.pi-theta1_21_p1p2)*180./np.pi
            theta2_21_p1p2 = (np.pi-theta2_21_p1p2)*180./np.pi

            # 3:2 | 2 angles 
            theta1_32_p1p2 = 3.0*l_2 -2.0*l_1 -omega_1
            theta2_32_p1p2 = 3.0*l_2 -2.0*l_1 -omega_2
            theta1_32_p1p2 = theta1_32_p1p2%(2.0*np.pi)
            theta2_32_p1p2 = theta2_32_p1p2%(2.0*np.pi)
            theta1_32_p1p2 = (np.pi-theta1_32_p1p2)*180/np.pi
            theta2_32_p1p2 = (np.pi-theta2_32_p1p2)*180/np.pi

            # 4:3 | 2 angles
            theta1_43_p1p2 = 4.0*l_2 -3.0*l_1 -omega_1
            theta2_43_p1p2 = 4.0*l_2 -3.0*l_1 -omega_2
            theta1_43_p1p2 = theta1_43_p1p2%(2.0*np.pi)
            theta2_43_p1p2 = theta2_43_p1p2%(2.0*np.pi)
            theta1_43_p1p2 = (np.pi-theta1_43_p1p2)*180./np.pi
            theta2_43_p1p2 = (np.pi-theta2_43_p1p2)*180./np.pi

            # 8:5 | 4 angles b c
            theta1_85_p1p2 = 8.0*l_2 -5.0*l_1 -3.0*omega_1
            theta2_85_p1p2 = 8.0*l_2 -5.0*l_1 -3.0*omega_2
            theta3_85_p1p2 = 8.0*l_2 -5.0*l_1 -2.0*omega_1 -omega_2
            theta4_85_p1p2 = 8.0*l_2 -5.0*l_1 -omega_1 -2.0*omega_2
            theta1_85_p1p2 = theta1_85_p1p2%(2.0*np.pi)
            theta2_85_p1p2 = theta2_85_p1p2%(2.0*np.pi)
            theta3_85_p1p2 = theta3_85_p1p2%(2.0*np.pi)
            theta4_85_p1p2 = theta4_85_p1p2%(2.0*np.pi)
            theta1_85_p1p2 = (np.pi-theta1_85_p1p2)*180./np.pi
            theta2_85_p1p2 = (np.pi-theta2_85_p1p2)*180./np.pi
            theta3_85_p1p2 = (np.pi-theta3_85_p1p2)*180./np.pi
            theta4_85_p1p2 = (np.pi-theta4_85_p1p2)*180./np.pi


            # 5:3 | 3 angles c d 
            theta1_53_p1p2 = 5.0*l_2 -3.0*l_1 -2.0*omega_1
            theta2_53_p1p2 = 5.0*l_2 -3.0*l_1 -2.0*omega_2
            theta3_53_p1p2 = 5.0*l_2 -3.0*l_1 -omega_1 -omega_2
            theta1_53_p1p2 = theta1_53_p1p2%(2.0*np.pi)
            theta2_53_p1p2 = theta2_53_p1p2%(2.0*np.pi)
            theta3_53_p1p2 = theta3_53_p1p2%(2.0*np.pi)
            theta1_53_p1p2 = (np.pi-theta1_53_p1p2)*180./np.pi
            theta2_53_p1p2 = (np.pi-theta2_53_p1p2)*180./np.pi
            theta3_53_p1p2 = (np.pi-theta3_53_p1p2)*180./np.pi

            fig = plt.figure(figsize=(13, 8))
            ligne = 3
            colonne = 2
            plot = 0
            labels = [-180, -90, 0, 90, 180]
            #-----------------------------------------------------------------------------------

            # --- plot 2:1
            plot+=1
            ax1 = fig.add_subplot(ligne,colonne,plot,)
            field = 'Res Angle 2:1 | planet '+key1 +' / '+key2
            ax1.plot(t1, theta1_21_p1p2, 'x', label='angle 1')
            ax1.plot(t1, theta2_21_p1p2, 'x', label='angle 2')
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.set_yticks(labels)
            ax1.set_yticklabels(labels)

            # --- plot 3:2
            plot+=1
            ax1 = fig.add_subplot(ligne,colonne,plot,)
            field = 'Res Angle 3:2 | planet  '+key1 +' / '+key2
            ax1.plot(t1, theta1_32_p1p2, 'x', label='angle 1')
            ax1.plot(t1, theta2_32_p1p2, '+', label='angle 2', alpha=.3)
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.set_yticks(labels)
            ax1.set_yticklabels(labels)

            # --- plot 4:3
            plot+=1
            ax1 = fig.add_subplot(ligne,colonne,plot,)
            field = 'Res Angle 4:3 | planet  '+key1 +' / '+key2
            ax1.plot(t1, theta1_43_p1p2, 'x', label='angle 1')
            ax1.plot(t1, theta2_43_p1p2, '+', label='angle 2', alpha=.3)
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.set_yticks(labels)
            ax1.set_yticklabels(labels)

            # --- plot 8:5
            plot+=1
            ax1 = fig.add_subplot(ligne,colonne,plot,)
            field = 'Res Angle 8:5 | planet  '+key1 +' / '+key2
            ax1.plot(t1, theta1_85_p1p2, 'x', label='angle 1')
            ax1.plot(t1, theta2_85_p1p2, '+', label='angle 2', alpha=.4)
            ax1.plot(t1, theta3_85_p1p2, 'x', label='angle 2', alpha=.3)
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.set_yticks(labels)
            ax1.set_yticklabels(labels)

            # --- plot 5:3
            plot+=1
            ax1 = fig.add_subplot(ligne,colonne,plot,)
            field = 'Res Angle 5:3 | planet  '+key1 +' / '+key2
            ax1.plot(t1, theta1_53_p1p2, 'x', label='angle 1')
            ax1.plot(t1[::100], theta2_53_p1p2[::100], '+', label='angle 2', alpha=.4)
            ax1.plot(t1[::100], theta3_53_p1p2[::100], 'x', label='angle 2', alpha=.3)
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.set_yticks(labels)
            ax1.set_yticklabels(labels)

            #-----------------------------------------------------------------------------------
            plt.tight_layout()

            output_figure_dirname = os.path.dirname(filename)
            output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_2planet_res_angles_" +key1 +"_" +key2 +".png")
            plt.savefig(output_figure_filename)
            plt.close()
            print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_11: # Resonance modes lmpq frequency angle two bodies

        fig = plt.figure(figsize=(16, 45))
        ligne = 7
        colonne = 2
        plot = 0
        #-----------------------------------------------------------------------------------


        for key in planets_keys:

            plot+=1
            ax1 = fig.add_subplot(ligne,colonne,plot,)

            planet_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]
            time = planet_data['current_time']

            frequ_mode_220 = planet_computed_data['frequ_mode_220q']
            period_mode_220 = planet_computed_data['period_mode_220q']
            order = planet_computed_data['order']
            print('   planet ', key, ' | order ', order)
            
            field = 'omega 220q Planet '+key
            al1=1.
            for q, mode in enumerate(frequ_mode_220):
                # print('   in mode q: ', q -order)
                ax1.plot(time, mode, label='220 '+str(q-7+order), alpha=al1)
                al1-=0.05
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.set_yscale('log')
            ax1.legend(loc=0, prop={'size':8})

            plot+=1
            ax1 = fig.add_subplot(ligne,colonne,plot,)            
            field = 'omega 220q Planet '+key
            al2 = 1.
            for q, mode in enumerate(period_mode_220):
                # print('   in mode q: ', q)
                ax1.plot(time, mode, label='220 '+str(q-7+order) +' | '+str(round(median(mode),5)) +' | '+str(round(median(mode/365.25),5)), alpha=al2)
                al2-=0.05
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.set_yscale('log')
            ax1.legend(loc=0, prop={'size':8})

        #-----------------------------------------------------------------------------------
        plt.tight_layout()

        output_figure_dirname = os.path.dirname(filename)
        output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_reso_angles_spin_meanmotion" +".png")
        plt.savefig(output_figure_filename)
        plt.close()
        # plt.show()
        print("> Output figure file written to: {}".format(output_figure_filename))

    if plot_12: # Print Hansen coeff 
        
        lab={}
        lab[1]='T1-b'
        lab[2]='T1-c'
        lab[3]='T1-d'
        lab[4]='T1-e'
        lab[5]='T1-f'
        lab[6]='T1-g'
        lab[7]='T1-h'

        for key in planets_keys:
            planet_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]
            frequ_mode_201q = planet_computed_data['frequ_mode_201q']
            frequ_mode_220q = planet_computed_data['frequ_mode_220q']
            amp_mode_201q = planet_computed_data['amp_mode_201q']
            amp_mode_220q = planet_computed_data['amp_mode_220q']
            order = planet_computed_data['order']
            print("--- Control modes --- Planet ", lab[int(key)])
            print(f"|\t mode 201q | order {order} | ")
            q = (7-order)
            for fr_mode, amp_mode in zip(frequ_mode_201q, amp_mode_201q):
                print(f"|\t > {q} \n|\t frequ {(fr_mode)} \n|\t amp {(amp_mode)} ")
                q+=1

            #----------------------------------
            fig = plt.figure(figsize=(10, 10))
            ligne = 2
            colonne = 1
            fig, ax = plt.subplots(ligne, colonne, figsize=(15, 15))

            
            table_frequ = frequ_mode_201q + frequ_mode_220q
            table_amp = amp_mode_201q + amp_mode_220q
            table_period = [2.*np.pi /t for t in table_frequ]
            print("EVALUATE")
            print(table_frequ,'\n=\n', frequ_mode_201q,'\n+\n', frequ_mode_220q, len(table_frequ), len(table_period))
            columns = ( 'frequ (rad/s)', 'period (days)', 'amplitude' )
            rows1 = [ "201 "+str(idx) for idx in np.arange(-order,order+1)]
            rows2 = [ "220 "+str(idx) for idx in np.arange(-order,order+1)]
            rows = rows1 + rows2
            table = np.transpose((table_frequ, table_period, table_amp))

            the_table = ax[0].table(cellText=table, loc='center', colLabels=columns, rowLabels=rows)
            ax[0].axis('tight')
            ax[0].axis('off')

            
            naff_frequs = planet_computed_data['naff_frequs']
            naff_periods = planet_computed_data['naff_periods']
            naff_amps = planet_computed_data['naff_amps']
            rows_naff = [str(p) for p in np.arange(1,len(naff_amps)+1)]
            table = np.transpose((naff_frequs, naff_periods, naff_amps))

            print(f"| Rows naff {rows_naff}  {len(rows_naff)}")
            print(f"| naff {naff_periods}  {len(naff_periods)}")
            print(f"| naff {naff_amps}  {len(naff_amps)}")

            the_table = ax[1].table(cellText=table, loc='center', colLabels=columns, rowLabels=rows_naff)
            ax[1].axis('tight')
            ax[1].axis('off')

            #-----------------------------------------------------------------------------------
            # plt.tight_layout()
            output_figure_dirname = os.path.dirname(filename)
            output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_planet_"+str(key)+"_table_modes_data" +".png")
            plt.savefig(output_figure_filename, dpi=150)
            plt.close()
            # plt.show()
            print("> Output figure file written to: {}".format(output_figure_filename))

            table_frequ = [t for t in table_frequ]
            naff_frequs = [t for t in naff_frequs]
            results_frequ = table_frequ + naff_frequs 

            results_periods = table_period + naff_periods

            rows = rows+rows_naff

            result=np.array([
                                rows,
                                results_frequ,
                            ])

            # file_name = f'./results_Trappist1_{lab[int(key)]}_list_solutions.csv'
            # header = ['frequ', 'period', 'amplitude']
            # df = pd.DataFrame(results, columns=header)
            # df.to_csv(file_name, sep="\t" )

    if history_output :
        print(">> Preparing text output...")
        all_data = None
        # data['star_mass'] = star_data['mass']
        for key in planets_keys:
            planet_data = planets_data[key]
            planet_computed_data = planets_computed_data[key]
            data = pd.DataFrame(planet_data['current_time'], columns=['current_time'])
            data['planet'] = key
            data['planet_mass'] = planet_data['mass']
            data['planet_radius'] = planet_data['radius']
            data['semi-major_axis_AU'] = planet_data['semi-major_axis']
            data['eccentricity'] = planet_data['eccentricity']
            data['planet_mean_motion'] = planet_computed_data['mean_motion']
            data['obliquity'] = planet_computed_data['planet_obliquity']
            data['position_x'] = planet_data['position_x']
            data['position_y'] = planet_data['position_y']
            data['position_z'] = planet_data['position_z']
            data['velocity_x'] = planet_data['velocity_x']
            data['velocity_y'] = planet_data['velocity_y']
            data['velocity_z'] = planet_data['velocity_z']
            data['spin'] = planet_computed_data['omega']
            data['substellar_point'] = planet_computed_data['substellar_point']
            data['true_ano'] = planet_data['true_anomaly']
            data['mean_ano'] = planet_data['mean_anomaly']
            data['long_per'] = planet_data['longitude_of_perihelion']
            data['mean_longitude'] = planet_computed_data['mean_longitude']
            data['true_ano_cumul'] = planet_computed_data['true_anomaly_cumul']
            data['cartesian_angle_cumul'] = planet_computed_data['cartesian_position_angle']
            data['theta'] = planet_computed_data['theta']
            data['mean_ano_cumul'] = planet_computed_data['mean_ano_cumul']
            data['theta_true_ano'] = planet_computed_data['theta_true_ano'] 
            data['long_peri_rescaled'] = planet_computed_data['longitude_of_perihelion_rescaled'] 
            data['theta_true_ano_long_peri'] = planet_computed_data['theta_true_ano_long_peri']
            data['true_ano_longperi_corrected'] = planet_computed_data['true_ano_longperi_corrected']
            data['theta_mean_ano'] = planet_computed_data['theta_mean_ano'] 
            data['theta_cartesian_position_angle'] = planet_computed_data['theta_cartesian_position_angle']
            data['planet_rotation_period'] = planet_computed_data['planet_rotation_period']
            if all_data is None:
                all_data = data
            else:
                all_data = pd.concat((all_data, data))

        output_text_dirname = os.path.dirname(filename)
        output_text_filename = os.path.join(output_text_dirname, os.path.splitext(os.path.basename(filename))[0] + "_save.txt")
        all_data.to_csv(output_text_filename, sep="\t", index=False)

        print("> Output data written to plain text file: {}\n".format(output_text_filename))


