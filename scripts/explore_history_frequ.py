import os
import numpy as np
import pandas as pd
import matplotlib as mpl
# mpl.use('Agg')
# mpl.rcParams['agg.path.chunksize'] = 10000
import matplotlib.pyplot as plt
import argparse
import posidonius
import json

if __name__ == "__main__":
    print("\n \t Explore Frequency CTRL into history frequ")
    parser = argparse.ArgumentParser()
    parser.add_argument('start_case_filename', action='store', help='Filename with the initial conditions of the simulation (e.g., universe_integrator.json)')
    parser.add_argument('historic_snapshot_filename', action='store', help='Filename with the historic snapshots of the simulation (e.g., universe_integrator_history.bin)')

    args = parser.parse_args()

    universe_integrator_json = json.load(open(args.start_case_filename, "r"))

    filename = args.historic_snapshot_filename
    n_particles, data = posidonius.analysis.history.read(filename)

    # print("Tidal torque Z ", data['tidal_torque_z'][:50])
    # print("Tidal torque Z ", data['tidal_torque_z'][:50], len(data['tidal_torque_z']))
    # print( data['current_time'][:50], len(data['current_time']))
    # print("Sigma220_2", data["sigma220_2_excitative_frequency"][:50])

    most_massive_particle_index = universe_integrator_json['universe']['hosts']['index']['most_massive']
    print("Transforming positions/velocities to heliocentric coordinates using the most masssive particle at index '{}'...".format(most_massive_particle_index))
    star_data, planets_data, planets_keys = posidonius.analysis.history.classify(n_particles, data, reference_particle_index=most_massive_particle_index, discard_first_hundred_years=False)
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
        planet_data = planets_data[key]
        planet_computed_data = {}


        print("Tidal torque Z \n Max ", np.max(np.abs(planet_data['tidal_torque_z'])), '\n The min', np.min(np.abs(planet_data['tidal_torque_z'])), '\n The len ', len(planet_data['tidal_torque_z']), )#' \n The 50 first elements:', planet_data['tidal_torque_z'][:50], )

        print("orthogonal_component_of_the_tidal_force_due_to_planetary_tide Z \n Max ", np.max(np.abs(planet_data['orthogonal_component_of_the_tidal_force_due_to_planetary_tide'])), '\n The min', np.min(np.abs(planet_data['orthogonal_component_of_the_tidal_force_due_to_planetary_tide'])), '\n The len ', len(planet_data['orthogonal_component_of_the_tidal_force_due_to_planetary_tide']), )# ' \n The 50 first elements:', planet_data['orthogonal_component_of_the_tidal_force_due_to_planetary_tide'][:50], )

        print("radial_component_of_the_tidal_force \n Max ", np.max(np.abs(planet_data['radial_component_of_the_tidal_force'])), '\n The min', np.min(np.abs(planet_data['radial_component_of_the_tidal_force'])), '\n The len ', len(planet_data['radial_component_of_the_tidal_force']), )#' \n The 50 first elements:', planet_data['radial_component_of_the_tidal_force'][:50], )


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
        # #pseudo_synchronization_period  = 2.*np.pi / (pseudo_rot * (planet_data['semi-major_axis']*AU)**(-3./2.) * posidonius.constants.HOUR) # Hours
        # pseudo_synchronization_period  = 2.*np.pi / (pseudo_rot * (planet_data['semi-major_axis']*posidonius.constants.AU)**(-3./2.) * posidonius.constants.DAY) # Days
        # pseudo_synchronization_period[non_zero] = np.nan
        # planet_computed_data['pseudo_synchronization_period'] = pseudo_synchronization_period

        planet_pseudo_synchronization_period = posidonius.calculate_pseudo_synchronization_period(planet_data['semi-major_axis'], e, star_mass, planet_mass) # days
        planet_angular_frequency_pseudo_rot = posidonius.constants.TWO_PI/(planet_pseudo_synchronization_period) # days^-1

        if universe_integrator_json['universe']['consider_effects']['tides']:

            planet_computed_data['tidal_torque_z'] = planet_data['tidal_torque_z']
            planet_computed_data['orthogonal_component_of_the_tidal_force_due_to_planetary_tide'] = planet_data['orthogonal_component_of_the_tidal_force_due_to_planetary_tide']
            planet_computed_data['radial_component_of_the_tidal_force'] = planet_data['radial_component_of_the_tidal_force']

            planet_computed_data['dangular_momentum_dt_x'] = planet_data['dangular_momentum_dt_x']
            planet_computed_data['dangular_momentum_dt_y'] = planet_data['dangular_momentum_dt_y']
            planet_computed_data['dangular_momentum_dt_z'] = planet_data['dangular_momentum_dt_z']

            planet_computed_data['acceleration_x'] = planet_data['acceleration_x']
            planet_computed_data['acceleration_y'] = planet_data['acceleration_y']
            planet_computed_data['acceleration_z'] = planet_data['acceleration_z']

            planet_computed_data['spin'] = planet_data['spin']
            planet_computed_data['orbital_frequency'] = planet_data['orbital_frequency']

            planet_computed_data['sigma220_2_excitative_frequency'] = planet_data['sigma220_2_excitative_frequency']
            planet_computed_data['sigma220_1_excitative_frequency'] = planet_data['sigma220_1_excitative_frequency']
            planet_computed_data['sigma2200_excitative_frequency'] = planet_data['sigma2200_excitative_frequency']
            planet_computed_data['sigma2201_excitative_frequency'] = planet_data['sigma2201_excitative_frequency']
            planet_computed_data['sigma2202_excitative_frequency'] = planet_data['sigma2202_excitative_frequency']

            planet_computed_data['im_love_number_sigma220_2'] = planet_data['im_love_number_sigma220_2']
            planet_computed_data['im_love_number_sigma220_1'] = planet_data['im_love_number_sigma220_1']
            planet_computed_data['im_love_number_sigma2200'] = planet_data['im_love_number_sigma2200']
            planet_computed_data['im_love_number_sigma2201'] = planet_data['im_love_number_sigma2201']
            planet_computed_data['im_love_number_sigma2202'] = planet_data['im_love_number_sigma2202']

            # planet_computed_data['re_love_number_sigma220_2'] = planet_data['re_love_number_sigma220_2']
            # planet_computed_data['re_love_number_sigma220_1'] = planet_data['re_love_number_sigma220_1']
            # planet_computed_data['re_love_number_sigma2200'] = planet_data['re_love_number_sigma2200']
            # planet_computed_data['re_love_number_sigma2201'] = planet_data['re_love_number_sigma2201']
            # planet_computed_data['re_love_number_sigma2202'] = planet_data['re_love_number_sigma2202']
            
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
            dissipation_factor_scale = universe_integrator_json['universe']['particles'][int(key)]['tides']['parameters']['input']['dissipation_factor_scale']
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


        ################################################################################
        #G         =  6.6742367e-11            # m^3.kg^-1.s^-2
        #AU        =  1.49598e11               # m
        #Msun      =  1.98892e30               # kg
        planet_orbital_period = 2*np.pi*np.sqrt(np.power(planet_data['semi-major_axis']*posidonius.constants.AU, 3)/(posidonius.constants.G_SI*posidonius.constants.M_SUN*(star_data['mass']+planet_data['mass'])))
        # planet_orbital_period = 2*np.pi*np.sqrt(np.power(planet_data['semi-major_axis']*posidonius.constants.AU, 3)/(posidonius.constants.G_SI*posidonius.constants.M_SUN*(star_data['mass'])))
        planet_computed_data['orbital_period'] = planet_orbital_period/(posidonius.constants.DAY) # days
        ################################################################################

        planet_orbital_frequ = 2*np.pi / ( planet_orbital_period )
        planet_rotation_frequ = planet_norm_spin /(posidonius.constants.DAY)

        # planet_computed_data['sigma_excitative_frequency_calculated'] = 2.0*( planet_orbital_frequ -planet_rotation_frequ) + 1.0*planet_orbital_frequ

        planet_computed_data['sigma220_2_excitative_frequency_calculated'] = 2.0*( planet_orbital_frequ -planet_rotation_frequ) - 2.0*planet_orbital_frequ
        planet_computed_data['sigma220_1_excitative_frequency_calculated'] = 2.0*( planet_orbital_frequ -planet_rotation_frequ) - 1.0*planet_orbital_frequ
        planet_computed_data['sigma2200_excitative_frequency_calculated'] = 2.0*( planet_orbital_frequ -planet_rotation_frequ)
        planet_computed_data['sigma2201_excitative_frequency_calculated'] = 2.0*( planet_orbital_frequ -planet_rotation_frequ) + 1.0*planet_orbital_frequ
        planet_computed_data['sigma2202_excitative_frequency_calculated'] = 2.0*( planet_orbital_frequ -planet_rotation_frequ) + 2.0*planet_orbital_frequ

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



    print("Recover Love number")
    # planet_data = np.loadtxt('Results_Trappist_1e_05.0_150_freq_Imk2_posidonius.txt',comments='#')
    # //Display
    planet_data = np.loadtxt('Results_Trappist_1e_00.0_138_freq_Imk2_posidonius.txt',comments='#')
    planet_mass, planet_radius, planet_gyration_radius = planet_data[0,:] 

    w_lmpq = planet_data [1:,0]
    ImK2 = planet_data [1:,1]
    ReK2 = planet_data [1:,2]
    size = np.size(w_lmpq)
    print("nmbre", size)
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



# ######################################################################################################################



    print("Preparing text output...")
    all_data = None
    for key in planets_keys:
        planet_data = planets_data[key]
        data = pd.DataFrame(planet_data['current_time'], columns=['current_time'])
        data['planet_rotation_frequ'] = 2.*np.pi/ (planets_computed_data[key]['planet_rotation_period']*posidonius.constants.DAY)
        data ['planet_orbital_frequ'] = 2.*np.pi/ (planets_computed_data[key]['orbital_period']*posidonius.constants.DAY)
        data ['pseudo_synchro_frequ'] = 2.*np.pi/ (2.*np.pi / (pseudo_rot * (planet_data['semi-major_axis']*posidonius.constants.AU)**(-3./2.) * posidonius.constants.DAY)*posidonius.constants.DAY)
        data['Tidal-torque_z'] = planet_data['tidal_torque_z']
        # data['planet'] = key
        # data['semi-major_axis_AU'] = planet_data['semi-major_axis']
        # data['corotation_radius_AU'] = planets_computed_data[key]['corotation_radius']
        # data['planet_obliquity_deg'] = planets_computed_data[key]['planet_obliquity']
        # data['eccentricity'] = planet_data['eccentricity']
        # data['inclination_deg'] = planet_data['inclination'] * (180 / np.pi)
        # data['energy_lost_due_to_tides_W_per_m2'] = planets_computed_data[key]['inst_tidal_flux']
        # data['mean_energy_lost_due_to_tides_W_per_m2'] = planets_computed_data[key]['mean_tidal_flux']
        # data['planet_rotation_period_hours'] = planets_computed_data[key]['planet_rotation_period']*24
        # data['planet_Pseudo_rotation_period_hours'] = planet_pseudo_synchronization_period

        # data['planet_pseudo_synchronization_period'] = planets_computed_data[key]['pseudo_synchronization_period']
        # data['energy_lost_due_to_tides_W'] = planets_computed_data[key]['denergy_dt']
        # data['mean_energy_lost_due_to_tides_W'] = planets_computed_data[key]['gravitational_energy_lost']
        # data['star_obliquity_deg'] = planets_computed_data[key]['star_obliquity']
        # data['planet_precession_angle_deg'] = planets_computed_data[key]['planet_precession_angle']
        # data['conservation_of_angular_momentum'] = conservation_of_angular_momentum
        # data['star_rotation_period_days'] = star_rotation_period
        # data['conservation_of_energy'] = relative_energy_error
        # data['sigma220_2_excitative_frequency'] = planet_computed_data['sigma220_2_excitative_frequency']
        # data['sigma220_2_excitative_frequency_calculated'] = planet_computed_data['sigma220_2_excitative_frequency_calculated']

        # data['sigma220_1_excitative_frequency'] = planet_computed_data['sigma220_1_excitative_frequency']
        # data['sigma220_1_excitative_frequency_calculated'] = planet_computed_data['sigma220_1_excitative_frequency_calculated']

        # data['sigma2200_excitative_frequency'] = planet_computed_data['sigma2200_excitative_frequency']
        # data['sigma2200_excitative_frequency_calculated'] = planet_computed_data['sigma2200_excitative_frequency_calculated']

        # data['sigma2201_excitative_frequency'] = planet_computed_data['sigma2201_excitative_frequency']
        # data['sigma2201_excitative_frequency_calculated'] = planet_computed_data['sigma2201_excitative_frequency_calculated']

        # data['sigma2202_excitative_frequency'] = planet_computed_data['sigma2202_excitative_frequency']
        # data['sigma2202_excitative_frequency_calculated'] = planet_computed_data['sigma2202_excitative_frequency_calculated']

        # data['spin'] = planet_computed_data['spin']
        # data['spin_calculated'] = planet_rotation_frequ

        # data['orbital_frequency'] = planet_data['orbital_frequency']
        # data['orbital_frequency_calculated'] = planet_orbital_frequ

        if all_data is None:
            all_data = data
        else:
            all_data = pd.concat((all_data, data))

    output_text_dirname = os.path.dirname(filename)
    output_text_filename = os.path.join(output_text_dirname, os.path.splitext(os.path.basename(filename))[0] + "_CTRL.txt")
    all_data.to_csv(output_text_filename, sep="\t", index=False)

    print("> Output data written to plain text file: {}".format(output_text_filename))

    # print("Tidal torque Z ", planet_data['tidal_torque_z'][:50], len(planet_data['tidal_torque_z']))
    # print("Radial Force ", planet_data['radial_component_of_the_tidal_force'])
    # print("Ortho Force ",planet_data['orthogonal_component_of_the_tidal_force_due_to_planetary_tide'] )
    # print("Spin X ", planet_data['spin_x'], 'Spin Y ', planet_data['spin_y'], "Spin Z",  planet_data['spin_z'])
    # print( planet_data['current_time'][:50], len(planet_data['current_time']))

    # ===================================================================================================================================

    print("preparing plot...")
    fig = plt.figure(figsize=(25, 20))
    ligne = 3
    colonne = 2
    i=0


    i=i+1
    ax = fig.add_subplot(ligne,colonne,i)
    field = ''
    for key in planets_keys:
        planet_data = planets_data[key]
        #ax.plot(w_lmpq, ImK2, label='fichier txt', color='blue')
        ax.scatter( planet_data['current_time'],planet_data['orthogonal_component_of_the_tidal_force_due_to_planetary_tide'] , linestyle='-', label='Ortho Force', color = 'red')
    ax.set_ylabel("")
    ax.set_ylim([-1e-16, 1e-16])
    ax.set_xscale('log')
    ax.set_xlim(left=1., right=5E5)
    ax.legend(loc=0, prop={'size':8})
    #plt.setp(ax.get_xticklabels(), visible=False)


    i=i+1
    ax = fig.add_subplot(ligne,colonne,i, sharex=ax)
    field = ''
    for key in planets_keys:
        planet_data = planets_data[key]
        #ax.plot(w_lmpq, ImK2, label='fichier txt', color='blue')
        ax.scatter( planet_data['current_time'], planet_data['spin_x'] ,ls='-', label='spin_x', color = 'red')
    ax.set_ylabel("")
    # ax.set_ylim([-1e-17, 1e-17])
    ax.set_xscale('log')
    ax.set_xlim(left=1., right=5E5)
    ax.legend(loc=0, prop={'size':8})
    #plt.setp(ax.get_xticklabels(), visible=False)

###########################################################################

    i=i+1
    ax = fig.add_subplot(ligne,colonne,i, sharex=ax)
    field = ''
    for key in planets_keys:
        planet_data = planets_data[key]
        #ax.plot(w_lmpq, ImK2, label='fichier txt', color='blue')
        ax.scatter( planet_data['current_time'],planet_data['radial_component_of_the_tidal_force'], ls='-',  label='Radial Force', color = 'red')
    ax.set_ylabel("")
    ax.set_ylim([-2e-13, 2e-13])
    ax.set_xscale('log')
    ax.set_xlim(left=1., right=5E5)
    ax.legend(loc=0, prop={'size':8})
    #plt.setp(ax.get_xticklabels(), visible=False)


    i=i+1
    ax = fig.add_subplot(ligne,colonne,i, sharex=ax)
    field = ''
    for key in planets_keys:
        planet_data = planets_data[key]
        #ax.plot(w_lmpq, ImK2, label='fichier txt', color='blue')
        ax.scatter( planet_data['current_time'], planet_data['spin_y'] ,ls='-', label='spin_y', color = 'red')
    ax.set_ylabel("")
    #ax.set_ylim([0.005, 0.028])
    ax.set_xscale('log')
    ax.set_xlim(left=1., right=5E5)
    ax.legend(loc=0, prop={'size':8})
    #plt.setp(ax.get_xticklabels(), visible=False)

###########################################################################


    i=i+1
    ax = fig.add_subplot(ligne,colonne,i, sharex=ax)
    field = ''
    for key in planets_keys:
        planet_data = planets_data[key]
        #ax.plot(w_lmpq, ImK2, label='fichier txt', color='blue')
        ax.scatter( planet_data['current_time'],planet_data['tidal_torque_z'] ,ls='-', label='Tidal Torque Z', color = 'red')
    ax.set_ylabel("")
    ax.set_ylim([-1e-18, 1e-18])
    # ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(left=1., right=5E5)
    ax.legend(loc=0, prop={'size':8})
    #plt.setp(ax.get_xticklabels(), visible=False)


    i=i+1
    ax = fig.add_subplot(ligne,colonne,i, sharex=ax)
    field = ''
    for key in planets_keys:
        planet_data = planets_data[key]
        #ax.plot(w_lmpq, ImK2, label='fichier txt', color='blue')
        ax.scatter( planet_data['current_time'], planet_data['spin_z'] ,ls='-', label='spin_z', color = 'red')
    ax.set_ylabel("")
    #ax.set_ylim([0.005, 0.028])
    ax.set_xscale('log')
    ax.set_xlim(left=1., right=5E5)
    ax.legend(loc=0, prop={'size':8})
    #plt.setp(ax.get_xticklabels(), visible=False)

###########################################################################



#     i=i+1
#     ax = fig.add_subplot(ligne,colonne,i)
#     field = ''
#     for key in planets_keys:
#         planet_data = planets_data[key]
#         #ax.plot(w_lmpq, ImK2, label='fichier txt', color='blue')
#         ax.scatter( planet_data['current_time'],planets_computed_data[key]['dangular_momentum_dt_x'] , label='Da_Dt_x', color = 'red')
#     ax.set_ylabel("")
#     #ax.set_ylim([0.005, 0.028])
#     ax.set_xscale('log')
#     ax.set_xlim(left=1., right=5E5)
#     ax.legend(loc=0, prop={'size':8})
#     #plt.setp(ax.get_xticklabels(), visible=False)

#     i=i+1
#     ax = fig.add_subplot(ligne,colonne,i)
#     field = ''
#     for key in planets_keys:
#         planet_data = planets_data[key]
#         #ax.plot(w_lmpq, ImK2, label='fichier txt', color='blue')
#         ax.scatter( planet_data['current_time'],planets_computed_data[key]['acceleration_x'] , label='acc_x', color = 'red')
#     ax.set_ylabel("")
#     #ax.set_ylim([0.005, 0.028])
#     ax.set_xscale('log')
#     ax.set_xlim(left=1., right=5E5)
#     ax.legend(loc=0, prop={'size':8})
#     #plt.setp(ax.get_xticklabels(), visible=False)
# ######################################################################################################################
#     i=i+1
#     ax = fig.add_subplot(ligne,colonne,i)
#     field = ''
#     for key in planets_keys:
#         planet_data = planets_data[key]
#         #ax.plot(w_lmpq, ImK2, label='fichier txt', color='blue')
#         ax.scatter( planet_data['current_time'],planets_computed_data[key]['dangular_momentum_dt_y'] , label='Da_Dt_y', color = 'red')
#     ax.set_ylabel("")
#     #ax.set_ylim([0.005, 0.028])
#     ax.set_xscale('log')
#     ax.set_xlim(left=1., right=5E5)
#     ax.legend(loc=0, prop={'size':8})
#     #plt.setp(ax.get_xticklabels(), visible=False)

#     i=i+1
#     ax = fig.add_subplot(ligne,colonne,i)
#     field = ''
#     for key in planets_keys:
#         planet_data = planets_data[key]
#         #ax.plot(w_lmpq, ImK2, label='fichier txt', color='blue')
#         ax.scatter( planet_data['current_time'],planets_computed_data[key]['acceleration_y'] , label='acc_y', color = 'red')
#     ax.set_ylabel("")
#     #ax.set_ylim([0.005, 0.028])
#     ax.set_xscale('log')
#     ax.set_xlim(left=1., right=5E5)
#     ax.legend(loc=0, prop={'size':8})
#     #plt.setp(ax.get_xticklabels(), visible=False)
# ######################################################################################################################
#     i=i+1
#     ax = fig.add_subplot(ligne,colonne,i)
#     field = ''
#     for key in planets_keys:
#         planet_data = planets_data[key]
#         #ax.plot(w_lmpq, ImK2, label='fichier txt', color='blue')
#         ax.scatter( planet_data['current_time'],planets_computed_data[key]['dangular_momentum_dt_z'] , label='Da_Dt_z', color = 'red')
#     ax.set_ylabel("")
#     #ax.set_ylim([0.005, 0.028])
#     ax.set_xscale('log')
#     ax.set_xlim(left=1., right=5E5)
#     ax.legend(loc=0, prop={'size':8})
#     #plt.setp(ax.get_xticklabels(), visible=False)

    
#     i=i+1
#     ax = fig.add_subplot(ligne,colonne,i)
#     field = ''
#     for key in planets_keys:
#         planet_data = planets_data[key]
#         #ax.plot(w_lmpq, ImK2, label='fichier txt', color='blue')
#         ax.scatter( planet_data['current_time'],planets_computed_data[key]['acceleration_z'] , label='axx_z', color = 'red')
#     ax.set_ylabel("")
#     #ax.set_ylim([0.005, 0.028])
#     ax.set_xscale('log')
#     ax.set_xlim(left=1., right=5E5)
#     ax.legend(loc=0, prop={'size':8})
#     #plt.setp(ax.get_xticklabels(), visible=False)

    # plt.show()
#     plt.tight_layout()

    output_figure_dirname = os.path.dirname(filename)
    output_figure_filename = os.path.join(output_figure_dirname, os.path.splitext(os.path.basename(filename))[0] + "_CTRL.png")
    plt.savefig(output_figure_filename)




