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
from sklearn.linear_model import LinearRegression
from sklearn import datasets, linear_model
# from scipy.integrate import  cumulative_trapezoid # cumtrapz
from pylaplace import LaplaceCoefficient
from numpy import diff
import sys
import time
import math
from matplotlib.ticker import MultipleLocator
from matplotlib.animation import FuncAnimation

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
    star_data, planets_data, planets_keys = posidonius.analysis.history.classify(n_particles, data, reference_particle_index=most_massive_particle_index, discard_first_years=False, first_years=4., discard_last_years=False, last_years=100.*365.25)
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
    ## Planets
    ################################################################################
    total_planets_angular_momentum_x = 0.
    total_planets_angular_momentum_y = 0.
    total_planets_angular_momentum_z = 0.
    total_planets_angular_momentum   = 0.
    total_planets_orbital_angular_momentum = 0.
    planets_computed_data = {}
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

    PI = np.pi
    # G = 1.
    G_SI = 6.67428e-11  #.   m^3.kg^-1.s^-2 (S.I. units)

    AU = 1.49597870700e11 #.   m

    DAY = 24.*3600 #.   s

    

    M_SUN =  1.9818e30 #.   Kg
    M_EARTH = 5.972e24 #.   Kg
    M_JUP = 1.898e27

    R_SUN = 696340000 # meter

    G = G_SI/(AU*AU*AU) * M_SUN * (DAY*DAY)


    # filename = "./data/r_21_a02_samem_1J1J_history.txt"
    # base_dir = os.path.dirname(filename)
    # filename_without_extension = os.path.splitext(os.path.basename(filename))[0]

    # raw_data = ascii.read(filename)
    # number = np.array(raw_data['planet'])
    # raw_t = np.array(raw_data['current_time'][np.where(number==1)])
    # raw_a = np.array(raw_data['semi-major_axis_AU'])
    # raw_e = np.array(raw_data['eccentricity'])
    # raw_w = np.array(raw_data['longitude_of_perihelion'])
    # raw_l = np.array(raw_data['mean_longitude'])
    # raw_Omega = 2*PI/np.array(raw_data['stellar_rotation_period'][np.where(number==1)]) # in DAYS^-1
    # raw_n = 2*PI/np.array(raw_data['orbital_period']) # in DAYS^-1
    # raw_delta_L = np.array(raw_data['conservation_of_angular_momentum'][np.where(number==1)])
    # raw_M = np.array(raw_data['mean_anomaly'])

    # key1 = '6'
    # key2 = '7'
    key1 = '1'
    key2 = '2'

    M = star_mass # in Msun
    m1 = planets_data[key1]['mass'][0] # in Msun
    m2 = planets_data[key2]['mass'][0] # in Msun
    # Omega_SI = raw_Omega/86400.
    # corot_r = pow((G_SI*(M*M_SUN)/(Omega_SI*Omega_SI)), 1./3.)/AU


    a_1 = planets_data[key1]['semi-major_axis'] # in AU
    e_1 = planets_data[key1]['eccentricity']
    l_1 = planets_computed_data[key1]['mean_longitude'] # in rad 
    w_1 = planets_data[key1]['longitude_of_perihelion'] # in rad
    # r_Omega_n_1 = raw_Omega/raw_n[np.where(number==1)]

    a_2 = planets_data[key2]['semi-major_axis'] # in AU
    e_2 = planets_data[key2]['eccentricity']
    l_2 = planets_computed_data[key2]['mean_longitude'] # in rad
    w_2 = planets_data[key2]['longitude_of_perihelion'] # in rad
    # r_Omega_n_2 = raw_Omega/raw_n[np.where(number==2)]

    print(f'| check a1 {a_1[0]} a2 {a_2[0]}')
    print(f'| check e1 {e_1[0]} e2 {e_2[0]}')
    print(f'| check l1 {l_1[0]} {l_1[0]%(2.0*np.pi)} l2 {l_2[0]} {l_2[0]%(2.0*np.pi)}')
    print(f'| check w1 {w_1[0]} {w_1[0]%(2.0*np.pi)} w2 {w_2[0]} {w_2[0]%(2.0*np.pi)}')

    raw_t = planets_data[key1]['current_time'] # in years

    length = len(raw_t)





    def MMR_1st_order(i,M, m1, m2, a1, a2, e1, e2, l1, l2, w1, w2, m=0):

        period_1 = planets_computed_data[key1]['orbital_period'][i] # 2*PI*np.sqrt(a1**3/(G*M))
        period_2 = planets_computed_data[key2]['orbital_period'][i] # 2*PI*np.sqrt(a2**3/(G*M))

        # A switch for which kind of first order resonance
        # simply set m = 0 and let the program decide itself    
        if m != 0:
            if round(period_1/(period_2-period_1)) != m:
                print(f"The periods do not corespond to the given resonance number: m(P)={round(period_1/(period_2-period_1))}, m={m}")
                return 1
        else:
            m = round(period_1/(period_2-period_1))
            print(f'| \t resonance order m={m} for period ratio {period_2/period_1} i.e. {(m+1)} / {m}')

            
        if m == 0:
            return (np.nan,) * 3

        eps1 = m1/M
        eps2 = m2/M
        epsp = eps1 + eps2

        p1 = -w1
        p2 = -w2

        alpha = a1/a2

        # Eq.2
        Lambda1 = m1 * np.sqrt(G*M*a1)
        Lambda2 = m2 * np.sqrt(G*M*a2)
        
        P1 = m1 *np.sqrt(G*M*a1)*(1.-np.sqrt(1.-e1**2.))
        P2 = m2 *np.sqrt(G*M*a2)*(1.-np.sqrt(1.-e2**2.))

        # Eq.9
        theta = (m+1.)*l2 - m*l1
        theta1 = l1 
        # Eq.12
        Theta = Lambda2/(m+1.)
        Theta1 = m/(m+1.)*Lambda2 + Lambda1
        # Eq.13
        x1 = np.sqrt(2.*P1)*np.cos(p1)
        x2 = np.sqrt(2.*P2)*np.cos(p2)
        y1 = np.sqrt(2.*P1)*np.sin(p1)
        y2 = np.sqrt(2.*P2)*np.sin(p2)
        # Change all action in units of \Theta_1
        x1_cap = x1/np.sqrt(Theta1)
        x2_cap = x2/np.sqrt(Theta1)
        y1_cap = y1/np.sqrt(Theta1)
        y2_cap = y2/np.sqrt(Theta1)

        # Nominal resonance location
        alpha_res = (m/(m+1.))**(2./3.)

        xi = m1/m2
        # Change all action in units of \Theta_1
        Theta_cap = Theta/Theta1
        Theta_bar = alpha_res /(m*(alpha_res + xi))

        s = np.sqrt(alpha/alpha_res)
        # Eq.21
        delta_Theta = xi/m *(1. - s) * 1./(alpha_res + xi + s*xi*(1.+xi/alpha_res))

        fm1_27 = [-1.19049, -2.02522, -2.84043 ,-3.64962 ,4.45614, -5.26125 ,-6.06552, -6.86925, -7.67261 , -8.47571 , -9.27861 , -10.0814, -10.884]
        fm1_31 = [0.42839, 2.48401, 3.28326, 4.08371, 4.88471, 5.68601, 6.48749, 7.28909, 8.09077, 8.89251, 9.69429, 10.4961, 11.2979]

        delta1 = 1./((m+1.) * Theta_cap)**2. * fm1_27[m-1]/(np.sqrt(1. - m*Theta_cap))
        delta2 = 1./((m+1.) * Theta_cap)**2. * fm1_31[m-1]/(np.sqrt((m+1.) * Theta_cap))

        delta1_bar = 1./((m+1.) * Theta_bar)**2. * fm1_27[m-1]/(np.sqrt(1. - m*Theta_bar))
        delta2_bar = 1./((m+1.) * Theta_bar)**2. * fm1_31[m-1]/(np.sqrt((1. + m)*Theta_bar))
        delta_bar = np.sqrt(delta1_bar**2 + delta2_bar**2)

        r1 = (delta1_bar*x1_cap + delta2_bar*x2_cap)/delta_bar
        s1 = (delta1_bar*y1_cap + delta2_bar*y2_cap)/delta_bar
        r2 = (delta2_bar*x1_cap - delta1_bar*x2_cap)/delta_bar
        s2 = (delta2_bar*y1_cap - delta1_bar*y2_cap)/delta_bar

        # Psi2 is not needed, noted that r2 and s2 are conserved
        # Eq.24
        Psi = 1./2. * (r1**2. + s1**2.)
        Psi2 = 1./2. * (r2**2. + s2**2.)
        psi = np.arctan2(s1, r1)
        
        phi = theta + psi
        gamma = -theta

        Phi = Psi
        Gamma = Psi - delta_Theta
        print(f'|\t Gamma {Gamma} = Psi {Psi} - dTehta {delta_Theta}')

        R = abs(fm1_27[m-1]/fm1_31[m-1])
        # Eq.31
        Q = epsp**(2./3.) * xi * (fm1_31[m-1]**2. / (9.*(m+1.)*(1.+xi)**2.))**(1./3.) * (alpha_res**(5./6.) / m) * ((R**2. + xi*np.sqrt(alpha_res))/(alpha_res+xi)**5.)**(1./3.)

        # "rescale the variables", that's what she said
        Phi_p = Phi/Q 
        # Eq.36
        Gamma_p = Gamma/Q
        X = np.sqrt(2*Phi_p) * np.cos(phi)
        Y = np.sqrt(2*Phi_p) * np.sin(phi)

        return X, Y, Gamma_p

    def Period_ratio(M, m1, m2, a1, a2):
        period_1 = 2*PI*np.sqrt(a1**3/(G*M))
        period_2 = 2*PI*np.sqrt(a2**3/(G*M))

        return period_2/period_1

    def find_X_fix(Gamma):
        # a specific function for solving X (Eq. 38 from Deck 2013)
        roots = np.roots([1., 0., -2*Gamma, 2.])
        roots = roots[np.isreal(roots)].real
        if len(roots) == 1:
            return np.array([roots[0], np.nan, np.nan])
        else:
            return roots
        
    def find_X_sol (Y, Gamma_p, H_p):
        roots = np.roots([1/4., 0., 0.5*(Y**2) - Gamma_p, 2, Y**4. + Gamma_p**2 + 2*H_p])
        roots = roots[np.isreal(roots)].real
        tem = np.full(4, np.nan)
        tem[:len(roots)] = roots
        return tem

    def Hamiltonian_prime (X, Y, Gamma_p):
        H_p = -0.5*(0.5*(X*X + Y*Y) - Gamma_p)**2 - X
        return H_p

    def func_Hamiltonian_solving_X (X, *data):
        Y, Gamma_p, H_sep = data
        return Hamiltonian_prime(X, Y, Gamma_p) - H_sep

    def print_progress(progress, i, length):
        # Construct the progress message
        message = f"Progress: {progress:.2f}%, Current frame: {i}/{length}"
        # Clear the current line
        sys.stdout.write('\r' + ' ' * (len(message) + 10))  # Overwrite with spaces
        sys.stdout.write('\r' + message)  # Carriage return and print
        sys.stdout.flush()

    X = np.zeros(length)
    Y = np.zeros(length)
    Gamma_p = np.zeros(length)
    P_ratio = np.zeros(length)
    X_fix = np.full((length, 3), np.nan)
    X_sol = np.zeros((length, 4))
    X_sol_tem = np.zeros((length, 4))
    X = np.zeros(length)
    Res_state = np.full(length, np.nan)

    separatrix_points = []
    # We manually input M=1 for stellar mass
    separatrix = False
    for i in range(length):
        X[i], Y[i], Gamma_p[i] = MMR_1st_order(i,M, m1, m2, a_1[i], a_2[i], e_1[i], e_2[i], l_1[i], l_2[i], w_1[i], w_2[i])
        print(f'| \t Gamma : {Gamma_p[i]}')
        
        if np.isnan(Gamma_p[i]):
            # Resonance = False
            P_ratio[i] = Period_ratio(M, m1, m2, a_1[i], a_2[i])
        else:
            X_fix[i] = find_X_fix(Gamma_p[i])
            X_fix[i] = np.sort(X_fix[i])
            H_p = Hamiltonian_prime (X[i], Y[i], Gamma_p[i])
            X_sol[i] = find_X_sol(0., Gamma_p[i], H_p)
            X_sol[i] = np.sort(X_sol[i])
            P_ratio[i] = Period_ratio(M, m1, m2, a_1[i], a_2[i])
            
    X3 = X_fix[:,2]
    X_cross_1 = -X3 + 2./np.sqrt(X3)
    X_cross_2 = -X3 - 2./np.sqrt(X3)
    # exit()


def resonance_series(frame_index):
    # Calculate the actual index (i) to use, scaling by 10 to get every 10th iteration
    global iteration_counter  # Reference the global counter
    i = valid_indices[frame_index]
    if (resonance_plot_real == True and Gamma_p[i] >= 3./2.): 
        ax.clear()
        Gamma_p_tem = Gamma_p[i]  # Example value, adjust as necessary
        X3_tem = X_fix[i, 2]   
        X_sep_Y0 = -X3_tem - 2./np.sqrt(X3_tem)
        X_cross_1_tem = -X3_tem + 2./np.sqrt(X3_tem)
        X_cross_2_tem = -X3_tem - 2./np.sqrt(X3_tem)
        H_sep_value = Hamiltonian_prime(X_sep_Y0, 0., Gamma_p_tem)
        
        x = np.linspace(-math.ceil(abs(X_cross_2_tem))*1.05, abs(X_cross_2_tem), 400)  # Adjust range and density as necessary
        y = np.linspace(-math.ceil(abs(X_cross_2_tem))*1.05, abs(X_cross_2_tem), 400)  # Adjust range and density as necessary
        X_plot, Y_plot = np.meshgrid(x, y)
        
        H = Hamiltonian_prime(X_plot, Y_plot, Gamma_p_tem)
    
        contour_outside = plt.contour(X_plot, Y_plot, H, levels=np.linspace(H_sep_value - np.abs(H_sep_value)*0.5, H_sep_value, 2), colors='black', linestyles='solid')
        contour_inside = plt.contour(X_plot, Y_plot, H, levels=np.linspace(H_sep_value, np.max(H), 5), colors='black', linestyles='solid')
        contour = ax.contour(X_plot, Y_plot, H, levels=[H_sep_value], colors='red', linewidths=2, linestyles='solid')
    
        # separatrix_paths = contour.collections[0].get_paths()
        separatrix_paths = contour.get_paths()  
        if separatrix_paths:
            separatrix_path = separatrix_paths[0]
            Res_state[i] = separatrix_path.contains_point((X[i], Y[i]))
            
        for x_sol in X_sol[i]:  # Iterate through solutions for this gamma_prime
            if separatrix_path.contains_point((X[i], Y[i])) == separatrix_path.contains_point((x_sol, 0.)):
                X_Y0_final[i] = x_sol
                break
        
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))  # e.g., minor ticks at 0.2 intervals
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.tick_params(axis='both', which='major', length=7, width=2)  # Major ticks
        ax.tick_params(axis='both', which='minor', length=4, width=1)   
        
        zorder_top = 10
        
        ax.scatter(X_fix[i, 0], 0., c="black", s=30, marker='o', zorder=zorder_top)
        ax.text(X_fix[i, 0], 0., 'X1', fontsize=20, verticalalignment='bottom', horizontalalignment='left')
        
        ax.scatter(X_fix[i, 1], 0., c="black", s=30, marker='o', zorder=zorder_top)
        ax.text(X_fix[i, 1], 0., 'X2', fontsize=20, verticalalignment='bottom', horizontalalignment='left')
        
        ax.scatter(X_fix[i, 2], 0., c="black", s=30, marker='o', zorder=zorder_top)
        ax.text(X_fix[i, 2], 0., 'X3', fontsize=20, verticalalignment='bottom', horizontalalignment='left')
        
        x_coords = [X_cross_1_tem, X_cross_2_tem]
        y_coords = [0, 0]
        
        ax.plot(x_coords, y_coords, linestyle='-.', color='blue', linewidth=1.5, zorder=5)
              
        time_tem = raw_t[i]
        title_text = r"$\Gamma^\prime = {:.2f}$, Time: {:.2e}".format(Gamma_p_tem, time_tem)
        ax.set_title(title_text, fontsize=22)
        
        ax.scatter(X[i], Y[i], c="purple", s=200, marker='*', zorder=10)
            
        ax.set_xlabel('X', fontsize=20)
        ax.set_ylabel('Y', fontsize=20)
        ax.tick_params(labelsize=20, direction='in', top=True, right=True, which='both')
        
        X_min = -math.ceil(abs(X_cross_2_tem))
        X_max = math.ceil(abs(X_cross_2_tem))
        ax.set_xlim(X_min, X_max)
        ax.set_ylim(X_min, X_max)
               
        total_valid_frames = len(valid_indices)
        progress = iteration_counter / total_valid_frames * 100  # Calculate progress percentage
        print_progress(progress, i, length)  
        
        iteration_counter += 1
        
    elif (resonance_plot_real == True and Gamma_p[i] <= 3./2.): 
        Res_state[iteration_counter] = False
    else:
        print("This script is switched off.")
        return

resonance_plot_real = True

if resonance_plot_real:    
    output_figure_dirname = os.path.dirname(filename)
    new_base_name = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_resonance")
    # new_base_name = filename_without_extension.replace('history', 'resonance')    
    res_gif = os.path.join(new_base_name+'_gif', new_base_name + ".gif")
    res_txt = os.path.join(new_base_name+'_text', new_base_name + ".txt")
    
    start_time = time.time()
    iteration_counter = 0
    
    valid_indices = [i for i, gamma in enumerate(Gamma_p) if gamma >= 3./2.]
    # it means we take a snapshot each 10 valid frames
    valid_indices = [i for i in valid_indices if (valid_indices.index(i) % 10 == 0)] 
    num_frames = len(valid_indices)
    
    fig, ax = plt.subplots(figsize=(10,8), dpi = 120)
    # fig.tight_layout()
    X_Y0_final = np.full(length, np.nan)
    anim = FuncAnimation(fig, resonance_series, frames=range(num_frames), repeat=False)
    
    print("Saving animation, please wait...")
    anim.save(res_gif, writer='Pillow', fps=24)
    end_time = time.time()  # Record the end time
    duration = (end_time - start_time)/60  # Calculate the duration
    print(f"\nAnimation saved. Total generation time: {duration:.2f} minutes.")
    
    data = pd.DataFrame()
    data['Gamma_p'] = Gamma_p
    data['X'] = X
    data['Y'] = Y
    data['X1_Y0'] = X_sol[:, 0]
    data['X2_Y0'] = X_sol[:, 1]
    data['X3_Y0'] = X_sol[:, 2]
    data['X4_Y0'] = X_sol[:, 3]
    data['X_Y0_final'] = X_Y0_final
    data['X1_fix'] = X_fix[:, 0]
    data['X2_fix'] = X_fix[:, 1]
    data['X3_fix'] = X_fix[:, 2]
    data['X_cross_1'] = X_cross_1
    data['X_cross_2'] = X_cross_2
    data['res_state'] = Res_state
    
    data.to_csv(res_txt, sep="\t", index=False)

resonance_plot_real = True

i = 500
for i in range(length):
    if (resonance_plot_real == True and Gamma_p[i] >= 3./2.): 
        fig, ax = plt.subplots(figsize=(10,8), dpi = 120, sharex=True)
        
        Gamma_p_tem = Gamma_p[i]  # Example value, adjust as necessary
        X3_tem = X_fix[i, 2]   
        X_sep_Y0 = -X3_tem - 2./np.sqrt(X3_tem)
        X_cross_1_tem = -X3_tem + 2./np.sqrt(X3_tem)
        X_cross_2_tem = -X3_tem - 2./np.sqrt(X3_tem)
        H_sep_value = Hamiltonian_prime(X_sep_Y0, 0., Gamma_p_tem)
        
        x = np.linspace(-math.ceil(abs(X_cross_2_tem))*1.05, abs(X_cross_2_tem), 400)  # Adjust range and density as necessary
        y = np.linspace(-math.ceil(abs(X_cross_2_tem))*1.05, abs(X_cross_2_tem), 400)  # Adjust range and density as necessary
        X_plot, Y_plot = np.meshgrid(x, y)
        
        H = Hamiltonian_prime(X_plot, Y_plot, Gamma_p_tem)

        contour_outside = ax.contour(X_plot, Y_plot, H, levels=np.linspace(H_sep_value - np.abs(H_sep_value)*0.5, H_sep_value, 2), colors='black', linestyles='solid')
        contour_inside = ax.contour(X_plot, Y_plot, H, levels=np.linspace(H_sep_value, np.max(H), 5), colors='black', linestyles='solid')
            
        contour = ax.contour(X_plot, Y_plot, H, levels=[H_sep_value], colors='red', linewidths=2, linestyles='solid')

        # separatrix_paths = contour.collections[0].get_paths()
        separatrix_paths = contour.get_paths() 
        if separatrix_paths:
            separatrix_path = separatrix_paths[0]
            
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))  # e.g., minor ticks at 0.2 intervals
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.tick_params(axis='both', which='major', length=7, width=2)  # Major ticks
        ax.tick_params(axis='both', which='minor', length=4, width=1)   
        
        ax.scatter(X_sol[i,0], 0., color = "orange", s = 150, marker = "x")
        ax.scatter(X_sol[i,1], 0., color = "orange", s = 150, marker = "x")
        ax.scatter(X_sol[i,2], 0., color = "orange", s = 150, marker = "x")
        ax.scatter(X_sol[i,3], 0., color = "orange", s = 150, marker = "x")
        
        point = (X_sol[i,0], 0.)  # The point you want to test
        inside = separatrix_path.contains_point(point)
        
        zorder_top = 10
        
        ax.scatter(X_fix[i, 0], 0., c="black", s=30, marker='o', zorder=zorder_top)
        ax.text(X_fix[i, 0], 0., 'X1', fontsize=20, verticalalignment='bottom', horizontalalignment='left')
        
        ax.scatter(X_fix[i, 1], 0., c="black", s=30, marker='o', zorder=zorder_top)
        ax.text(X_fix[i, 1], 0., 'X2', fontsize=20, verticalalignment='bottom', horizontalalignment='left')
        
        ax.scatter(X_fix[i, 2], 0., c="black", s=30, marker='o', zorder=zorder_top)
        ax.text(X_fix[i, 2], 0., 'X3', fontsize=20, verticalalignment='bottom', horizontalalignment='left')
        
        x_coords = [X_cross_1_tem, X_cross_2_tem]
        y_coords = [0, 0]
        
        ax.plot(x_coords, y_coords, linestyle='-.', color='blue', linewidth=1.5, zorder=5)
            
        time_tem = raw_t[i]
        
        title_text = r"$\Gamma^\prime = {:.2f}$, Time: {:.2e}".format(Gamma_p_tem, time_tem)
        ax.set_title(title_text, fontsize=22)
        
        ax.scatter(X[i], Y[i], c="purple", s=200, marker='*', zorder=10)
        print(i, X[i], Y[i])
        print(separatrix_path.contains_point((X[i], Y[i])))
            
        ax.set_xlabel('X', fontsize=20)
        ax.set_ylabel('Y', fontsize=20)
        ax.tick_params(labelsize=20, direction='in', top=True, right=True, which='both')
        
        X_min = -math.ceil(abs(X_cross_2_tem))
        X_max = math.ceil(abs(X_cross_2_tem))
        ax.set_xlim(X_min, X_max)
        ax.set_ylim(X_min, X_max)
        
        # ax.set_xlim(-math.ceil(abs(X_cross_2_tem))*1.1, abs(X_cross_2_tem))
        # ax.set_ylim(-math.ceil(abs(X_cross_2_tem))*1.1, math.ceil(abs(X_cross_2_tem)))

        # plt.show()
        plt.savefig(f'save_res_'+str(i)+'.png', dpi=100)
    elif Gamma_p[i] >= 3./2.:
        print("# Gamma < 3/2 = No plot") 
    elif resonance_plot_real == False:
        print("This script is switched off.")

    else:
        print("Not in resonance.")
        # plt.ioff()
