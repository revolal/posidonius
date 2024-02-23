import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys, os
from matplotlib import colors
import pandas as pd
from statistics import mean, median
from scipy.interpolate import interp1d


def calculate_tidal_excitation_frequency_mode_sigma_2mpq(m, p, q, spin_angle, mean_motion):
    return (2.0 -2.0*p +q)*mean_motion -m*spin_angle


def k2_data(file_name):
    data = np.genfromtxt(file_name)
    planet_mass, planet_radius, planet_gyration_radius = data[0, :]

    ReK2 = data[1:, 2]
    size = np.size(ReK2)
    w_lmpq = data[2:size, 0]
    ImK2 = data[2:size, 1]
    ReK2 = data[2:size, 2]
    return w_lmpq, ImK2, ReK2


def calculate_kaula_numbers(wk2, nm_data, real_part_love_number, imaginary_part_love_number, love_number_excitation_frequency):
    nm_data = int(nm_data)
    line_up, column_up, line_down, column_down = 0, 0, 0, 0
    in_the_file = True

    last_line, last_column = divmod(nm_data, 32) if nm_data % 32 == 0 else (nm_data % 32 - 1, nm_data // 32)
    case_oneline_last_col = True if last_line == 0 or last_line == 1 else False
    parity, w_k2 = (True, abs(wk2)) if wk2 < 0.0 else (False, abs(wk2))

    ctrl = True
    im_k2, re_k2 = 0.0, 0.0

    if w_k2 > love_number_excitation_frequency[last_column][last_line - 1]:
        if case_oneline_last_col:
            im_k2 = imaginary_part_love_number[last_column][last_line] \
                + (w_k2 - love_number_excitation_frequency[last_column][last_line]) \
                * (imaginary_part_love_number[last_column][last_line] - imaginary_part_love_number[last_column - 1][31]) \
                / (love_number_excitation_frequency[last_column][last_line] - love_number_excitation_frequency[last_column - 1][31])
            re_k2 = real_part_love_number[last_column][last_line] \
                + (w_k2 - love_number_excitation_frequency[last_column][last_line]) \
                * (real_part_love_number[last_column][last_line] - real_part_love_number[last_column - 1][31]) \
                / (love_number_excitation_frequency[last_column][last_line] - love_number_excitation_frequency[last_column - 1][31])
        else:
            delta = (w_k2 - love_number_excitation_frequency[last_column][last_line - 1]) \
                / (love_number_excitation_frequency[last_column][last_line - 1] - love_number_excitation_frequency[last_column][last_line - 2])
            im_k2 = delta * (imaginary_part_love_number[last_column][last_line - 1] - imaginary_part_love_number[last_column][last_line - 2]) + imaginary_part_love_number[last_column][last_line - 1]
            re_k2 = delta * (real_part_love_number[last_column][last_line - 1] - real_part_love_number[last_column][last_line - 2]) + real_part_love_number[last_column][last_line - 1]
            if im_k2 < 0.:
                im_k2 = imaginary_part_love_number[last_column][last_line] + (imaginary_part_love_number[last_column][last_line] - imaginary_part_love_number[last_column][last_line - 1]) / 2.
            if re_k2 < 0.:
                re_k2 = real_part_love_number[last_column][last_line] + (real_part_love_number[last_column][last_line] - real_part_love_number[last_column][last_line - 1]) / 2.

        column_up, line_up, column_down, line_down = last_column, last_line, 0, 0
        in_the_file = False

    elif w_k2 > love_number_excitation_frequency[last_column - 1][31]:
        if w_k2 < love_number_excitation_frequency[last_column][0]:
            delta = (w_k2 - love_number_excitation_frequency[last_column - 1][31]) \
                / (love_number_excitation_frequency[last_column][0] - love_number_excitation_frequency[last_column - 1][31])
            im_k2 = (1. - delta) * imaginary_part_love_number[last_column - 1][31] \
                + delta * imaginary_part_love_number[last_column][0]
            re_k2 = (1. - delta) * real_part_love_number[last_column - 1][31] \
                + delta * real_part_love_number[last_column][0]
            column_up, line_up, column_down, line_down = last_column, 0, last_column - 1, 31
        else:
            for frequency2 in range(last_line + 1):
                if love_number_excitation_frequency[last_column][frequency2] >= w_k2:
                    if love_number_excitation_frequency[last_column][frequency2] == w_k2:
                        im_k2 = imaginary_part_love_number[last_column][frequency2]
                        re_k2 = real_part_love_number[last_column][frequency2]
                        column_up, line_up, column_down, line_down = last_column, frequency2, last_column, frequency2
                        break
                    else:
                        delta = (w_k2 - love_number_excitation_frequency[last_column][frequency2 - 1]) \
                            / (love_number_excitation_frequency[last_column][frequency2] - love_number_excitation_frequency[last_column][frequency2 - 1])
                        im_k2 = (1. - delta) * imaginary_part_love_number[last_column][frequency2 - 1] \
                            + delta * imaginary_part_love_number[last_column][frequency2]
                        re_k2 = (1. - delta) * real_part_love_number[last_column][frequency2 - 1] \
                            + delta * real_part_love_number[last_column][frequency2]
                        column_up, line_up, column_down, line_down = last_column, frequency2, last_column, frequency2 - 1
                        break

    else:
        for frequency1 in range(last_column):
            if ctrl and love_number_excitation_frequency[frequency1][31] >= w_k2:
                if ctrl and love_number_excitation_frequency[frequency1][0] > w_k2:
                    delta = (w_k2 - love_number_excitation_frequency[frequency1 - 1][31]) \
                        / (love_number_excitation_frequency[frequency1][0] - love_number_excitation_frequency[frequency1 - 1][31])
                    im_k2 = (1. - delta) * imaginary_part_love_number[frequency1 - 1][31] \
                        + delta * imaginary_part_love_number[frequency1][0]
                    re_k2 = (1. - delta) * real_part_love_number[frequency1 - 1][31] \
                        + delta * real_part_love_number[frequency1][0]
                    column_up, line_up, column_down, line_down = frequency1, 0, frequency1 - 1, 31
                    ctrl = False
                elif ctrl:
                    for frequency2 in range(len(love_number_excitation_frequency)):
                        if ctrl and love_number_excitation_frequency[frequency1][frequency2] >= w_k2:
                            if love_number_excitation_frequency[frequency1][frequency2] == w_k2:
                                im_k2 = imaginary_part_love_number[frequency1][frequency2]
                                re_k2 = real_part_love_number[frequency1][frequency2]
                                column_up, line_up,


# Plot Sub-stellar comparison
# pwd="/Volumes/Expansion/Mac_data/posidonius/from_baobab/substellar_internal"
pwd="/Volumes/Expansion/baobab/posidonius"
data_dir_path =pwd+"/substellar/"
# cases=( 
# "t1_earth",
# # "t1_earth_bis",
# "t1_intermediate",
# # "t1_intermediate_bis",
# "t1_mercury",
# # "t1_mercury_bis",
# "t1_small",) 
# # "t1_small_bis")
# output_figure_filename_each_p = 'compare_substellar_each_planets_internal.png'
# output_figure_filename_mmr = 'compare_substellar_mmr_each_planets_internal.png'
# output_figure_filename_each_case = 'compare_substellar_each_internal_every_planet.png'
# output_figure_filename_each_case_omegan = 'compare_omegam_each_internal_every_planet.png'
# output_figure_filename_t1b = 'compare_substellar_t1b_internal.png'
# output_figure_filename_t1b_positions_angle_true_ano = 'compare_position_angle_true_ano_t1b_internal.png'
# output_figure_filename_all = 'compare_substellar_all_internal.png'

output_figure_filename_each_p = 'earth_compare_substellar_each_planets_internal.png'
output_figure_filename_mmr = 'earth_compare_substellar_mmr_each_planets_internal.png'
output_figure_filename_each_case = 'earth_compare_substellar_each_internal_every_planet.png'
output_figure_filename_each_case_omegan = 'earth_compare_omegam_each_internal_every_planet.png'
output_figure_filename_each_caseplanet_omegan = 'earth_compare_omegam_each_internal_every_caseplanet.png'
output_figure_filename_t1b = 'earth_compare_substellar_t1b_internal.png'
output_figure_filename_t1b_positions_angle_true_ano = 'earth_compare_position_angle_true_ano_t1b_internal.png'
output_figure_filename_all = 'earth_compare_substellar_all_internal.png'

# ###########################################################
# # # # Have to loop over the path and cases:
# data_dir_path_cases = [ data_dir_path +p for p in cases]
# print('| Cases loaded > ', data_dir_path_cases)
# res = []
# for data_dir_path_case in data_dir_path_cases:
#     for path in os.listdir(data_dir_path_case):
#         if path.endswith('.txt'):
#             print('|\t find > ', path, ' in ', data_dir_path_case)
#             res.append(path)
# # res.sort()
# print('| files loaded > ', res)
# cases_planets_data = []
# cases = []
# for path,file in zip(data_dir_path_cases,res):
#     data = path+file
#     print('|\t Look for ', data)
#     df = pd.read_csv(data, sep='\t', header=0, )
#     case = file[11:-17]
#     print('|\t load data from ',case, ' in file ', file)
#     cases.append(case)
#     n_particles = df['planet'].max()
#     print('|\t N particle ', n_particles)
#     planets_data = {}
#     planets_keys = [] # To ensure the order
#     for i in range(n_particles):
#         planets_data["{}".format(i+1)] = df[df['planet'] == i+1]
#         planets_keys.append("{}".format(i+1))
#     cases_planets_data.append(planets_data)
# print('| N cases loaded : ', len(cases_planets_data))

###########################################################
data_dir_path_cases = data_dir_path 
print('| Cases loaded > ', data_dir_path_cases)
res = []
for path in os.listdir(data_dir_path_cases):
    print(f'| Evaluate ', path)
    if path.endswith('.txt'):
        # if not path.endswith('bis_history_save.txt'):
        if path.endswith('earth_history_save.txt'):
        # if path.endswith('intermediate_history_save.txt'):
            print('|\t find > ', path)
            res.append(path)
res.sort(reverse=True)
print(f'| Find {len(res)} cases in ', res)

cases_planets_data = []
cases = []
for file in (res):
    data = data_dir_path_cases+file
    df = pd.read_csv(data, sep='\t', header=0, )
    case = file[:-17]
    print('|\t load data from ',case, ' in file ', file)
    cases.append(case)

    n_particles = df['planet'].max()
    print('|\t N particle ', n_particles)
    planets_data = {}
    planets_keys = [] # To ensure the order
    for i in range(n_particles):
        planets_data["{}".format(i+1)] = df[df['planet'] == i+1]
        planets_keys.append("{}".format(i+1))
    cases_planets_data.append(planets_data)

first_years = 0.
last_years = 0. #50/365.25
tmp = []
for data_cases in cases_planets_data:
    planets_data = data_cases
    for key in planets_keys:
        # print(f'|\t cut time for {key}')
        planet_data = planets_data[key]
        if first_years >0.:
            planet_data = planet_data[planet_data['current_time'] >= first_years]
        if last_years >0.:
            planet_data = planet_data[planet_data['current_time'] <= last_years]
        planets_data[key] = planet_data
    data_cases = planets_data
    tmp.append(data_cases)
cases_planets_data = tmp

print('| N cases loaded : ', len(cases_planets_data))

###########################################################

print("| Pepare to plot ...")

lab={}
lab[1]='T1-b'
lab[2]='T1-c'
lab[3]='T1-d'
lab[4]='T1-e'
lab[5]='T1-f'
lab[6]='T1-g'
lab[7]='T1-h'

line = {}
line[0] = "-"
line[1] = "-"
line[2] = "-"
line[3] = "-"
line[4] = "-"
line[5] = "-"
line[6] = "-"
line[7] = "-"
line[8] = "+"
line[9] = ">"
line[10] = "<"
line[11] = "^"
line[12] = "o"
line[13] = "v"
line[14] = "x"
line[15] = "+"


plot_control_spin_meanmotion = False
plot_control_k2 = False
plot_substellar = False
plot_ecc_data = True

if plot_control_k2:
    pwd = '/Users/revolal/Desktop/posidonius/input/love_numbers/trappist1/20230309_Data_Posidonius/'
    t1 = []
    file = "Results_Trappist1b_Fe90Si02_800K_visc_freq_Imk2_posidonius.txt"
    t1.append(k2_data(pwd+file))
    file = "Results_Trappist1c_Fe90Si02_800K_visc_freq_Imk2_posidonius.txt"
    t1.append(k2_data(pwd+file))
    file = "Results_Trappist1d_Fe90Si02_600K_visc_freq_Imk2_posidonius.txt"
    t1.append(k2_data(pwd+file))
    file = "Results_Trappist1e_Fe90Si02_300K_visc_freq_Imk2_posidonius.txt"
    t1.append(k2_data(pwd+file))
    file = "Results_Trappist1f_Fe90Si02_300K_visc_freq_Imk2_posidonius.txt"
    t1.append(k2_data(pwd+file))
    file = "Results_Trappist1g_Fe90Si02_300K_visc_freq_Imk2_posidonius.txt"
    t1.append(k2_data(pwd+file))
    file = "Results_Trappist1h_Fe90Si02_300K_visc_freq_Imk2_posidonius.txt"
    t1.append(k2_data(pwd+file))

    for case, data_cases in zip(cases,cases_planets_data):
        kk = 0
        
        for key in planets_keys:
            print(f'|\t for planet {key}')
            ecc = planet_data['eccentricity']

            frequ, Im, Re = t1[kk]
            x = frequ
            y = Im
            f = interp1d(x, y, kind='linear', fill_value="extrapolate")

            order = np.zeros(len(ecc))
            order7 = ecc > 0.3
            order6 = ecc > 0.25
            order5 = ecc > 0.2
            order4 = ecc > 0.15
            order3 = ecc > 0.1
            order2 = ecc <= 0.1
            if np.any(order7):
                order[order7] = 7.
            if np.any(order6):
                order[order6] = 6.
            if np.any(order5):
                order[order5] = 5.
            if np.any(order4):
                order[order4] = 4.
            if np.any(order3):
                order[order3] = 3.
            if np.any(order2):
                order[order2] = 2.
            print(order)
            print(order[order2])
            eccentricity = ecc
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

            q_min = -order
            q_max = order
            print(f'|\t > order : {order} | q from {q_min} to {q_max}  || {(7-order)} {(7+order)} ')
            # for ordr in np.arange(0, 7):
                # print(f'|\t > idx : {(7-ordr)} {(7+ordr)} ')

            # Spin and mean_motion are in rad/days
            spin = planet_data['spin']/86400.
            mean_motion = planet_data['planet_mean_motion']/86400.
            frequ_201q = []
            imk_201q = []
            frequ_220q = []
            imk_220q = []
            frequ_mode_201q = []
            frequ_mode_220q = []
            for q in np.arange(q_min, q_max+1):
                print(f'|\t  > - compute order {q} ')
                w_201q = calculate_tidal_excitation_frequency_mode_sigma_2mpq(0, 1, q, spin, mean_motion)
                w_220q = calculate_tidal_excitation_frequency_mode_sigma_2mpq(2, 0, q, spin, mean_motion)
                
                frequ_mode_201q.append(median(w_201q))
                frequ_mode_220q.append(median(w_220q))

                frequ_201q.append((w_201q))
                frequ_220q.append((w_220q))
                interpol_201 = f(w_201q)
                interpol_220 = f(w_220q)
                imk_201q.append((interpol_201))
                imk_220q.append((interpol_220))

            fig = plt.figure(figsize=(10,6))
            ligne = 1
            colonne = 1
            plot = 0

            plot+=1
            ax1 = fig.add_subplot(ligne,colonne,plot)
        
            print(f'|\t Plot | {case}')
            to_label = f' {lab[int(key)]} ' 
            planets_data = data_cases
            planet_data = planets_data[key]

            
            
            ax1.plot(frequ, Im, label = 'planet ' +key)
            for (fr, imk2) in zip(frequ_201q, imk_201q):
                ax1.plot(fr,imk2 ,'+', color='red')
            for (fr, imk2) in zip(frequ_220q, imk_220q):
                ax1.plot(fr,imk2 ,'x', color='green')
            ax1.set_xlabel("Excitation frequency rad.s$^{-1}$", fontsize=15)
            ax1.set_ylabel(f" T1 {key}  $\Im (k_2)$", fontsize=15)
            for (q, f_201) in zip(np.arange(q_min, q_max+1),frequ_mode_201q):
                ax1.axvline(x = np.abs(f_201), ls=':', color='black', label=f' frequ {q} {f_201} ')
            for (q,f_220) in zip(np.arange(q_min, q_max+1),frequ_mode_220q):
                ax1.axvline(x = np.abs(f_220), ls='-.', color='black',  label=f' frequ {q} {f_220} ')
            ax1.set_yscale('log')
            ax1.set_xscale('log')
            plt.yticks(fontsize=12)
            plt.xticks(fontsize=12)
            ax1.legend(loc=0, prop={'size':5})
            kk+=1

            plt.tight_layout()
            # plt.show()
            save_file_name = f'evaluate_imk2_case_{case}_planet_{key}.png'
            plt.savefig(save_file_name, dpi=200)
            plt.close()
            print("| Output figure file written to: {}".format(save_file_name))


if plot_control_spin_meanmotion:
    fig = plt.figure(figsize=(10,6))#*len(planets_keys)))
    ligne = 1#len(planets_keys) 
    colonne = 1
    plot = 0
    for case, data_cases in zip(cases,cases_planets_data):
        c=0
        al=1
        plot+=1
        ax1 = fig.add_subplot(ligne,colonne,plot)
        for key in planets_keys:
            
            if key!='9':
            
                print(f'|\t Plot / label {case}')
                to_label = f' {lab[int(key)]} ' 
                planets_data = data_cases
                planet_data = planets_data[key]
                time = [t*365.25 for t in planet_data['current_time']]#[planet_data['current_time']<50./365.25]]
                spin = planet_data['spin']#[planet_data['current_time']<50./365.25]
                mean_motion = planet_data['planet_mean_motion']#[planet_data['current_time']<50./365.25]
                omega_n = spin / mean_motion
                field = 'rad/day' #case
                # ax1.plot(time, omega_n, line[c], label=to_label, alpha=al)
                ax1.plot(time, mean_motion, line[c], label='Mean motion n', alpha=al)
                ax1.plot(time, spin, line[c], label='Spin $\Omega$', alpha=al)
                print(f'|\t CASE {to_label} MIN {min(omega_n)} \t MAX {max(omega_n)} \t MED {median(omega_n)}')
                # al-=0.05
                c+=1
                ax1.set_ylabel(field, fontsize=16)
                ax1.set_xlabel('time (days)', fontsize=16)
                plt.yticks(fontsize=16)
                plt.xticks(fontsize=16)
                # if key=='7':
                ax1.legend(loc=0, prop={'size':12})
                # ax1.axhline(y=1., ls=':', color='black')
                # y_ticks = [0.9985, 0.999, 0.9995 ,1., 1.0005, 1.001, 1.0015]  # Adjust the range as needed
                # y_tick_labels = ['0.9985', '0.999', '0.9995', '1.0', '1.0005', '1.001', '1.0015']
                # plt.yticks(y_ticks, y_tick_labels)
            # ax1.axhline(y=180., ls=':', color='black')
            # ax1.axhline(y=360., ls=':', color='black')
            # ax1.axhline(y=540., ls=':', color='black')
    # ax1.legend(loc=0, prop={'size':10})
    plt.tight_layout()
    # plt.show()
    plt.savefig('paper_substellar_omega_n_t1.png', dpi=200)
    plt.close()
    print("| Output figure file written to: {}".format(output_figure_filename_each_caseplanet_omegan))

    fig = plt.figure(figsize=(10,6*1))#len(planets_keys)))
    ligne = 1#len(planets_keys) 
    colonne = 1
    plot = 0
    for case, data_cases in zip(cases,cases_planets_data):
        c=0
        al=1
        for key in planets_keys:
            if key =='1':
                plot+=1
                ax1 = fig.add_subplot(ligne,colonne,plot)
                print(f'|\t Plot / label {case}')
                to_label = f' {lab[int(key)]} ' 
                planets_data = data_cases
                planet_data = planets_data[key]
                time = [t*365.25 for t in planet_data['current_time'][planet_data['current_time']<150./365.25]]
                spin = planet_data['spin'][planet_data['current_time']<150./365.25]
                mean_motion = planet_data['planet_mean_motion'][planet_data['current_time']<150./365.25]
                omega_n = spin / mean_motion
                field = 'Spin state $\Omega$ /n' #case
                ax1.plot(time, mean_motion, line[c], label=to_label, alpha=al)
                ax1.plot(time, spin, '+', label=to_label, alpha=al)
                al-=0.08
                c+=1
                ax1.set_ylabel(field, fontsize=16)
                ax1.set_xlabel('days', fontsize=16)
                plt.yticks(fontsize=16)
                plt.xticks(fontsize=16)
                ax1.legend(loc=0, prop={'size':10})
                # ax1.axhline(y=1., ls=':', color='black')
            # ax1.axhline(y=180., ls=':', color='black')
            # ax1.axhline(y=360., ls=':', color='black')
            # ax1.axhline(y=540., ls=':', color='black')
    # plt.show()
    ax1.legend(loc=0, prop={'size':10})
    plt.tight_layout()
    # plt.savefig(output_figure_filename_each_caseplanet_omegan)
    plt.close()
    print("| Output figure file written to: {}".format(output_figure_filename_each_caseplanet_omegan))

    fig = plt.figure(figsize=(10,6*1))#len(planets_keys)))
    ligne = 1#len(planets_keys) 
    colonne = 1
    plot = 0
    for case, data_cases in zip(cases,cases_planets_data):
        c=0
        al=1
        for key in planets_keys:
            if key =='1':
                plot+=1
                ax1 = fig.add_subplot(ligne,colonne,plot)
                print(f'|\t Plot / label {case}')
                to_label = f' {lab[int(key)]} ' 
                planets_data = data_cases
                planet_data = planets_data[key]
                time = [t*365.25 for t in planet_data['current_time'][planet_data['current_time']<150./365.25]]
                true_ano_cumul = planet_data['true_ano_cumul'][planet_data['current_time']<150./365.25]
                cartesian_angle_cumul = [p*np.pi/180. for p in planet_data['cartesian_angle_cumul'][planet_data['current_time']<150./365.25]]
                field = 'angles' +case
                ax1.plot(time, true_ano_cumul, line[c], label=to_label, alpha=al)
                ax1.plot(time, cartesian_angle_cumul, '+', label=to_label, alpha=al)
                al-=0.08
                c+=1
                ax1.set_ylabel(field, fontsize=16)
                ax1.set_xlabel('days', fontsize=16)
                plt.yticks(fontsize=16)
                plt.xticks(fontsize=16)
                ax1.legend(loc=0, prop={'size':10})
                # ax1.axhline(y=1., ls=':', color='black')
            # ax1.axhline(y=180., ls=':', color='black')
            # ax1.axhline(y=360., ls=':', color='black')
            # ax1.axhline(y=540., ls=':', color='black')
    plt.show()
    ax1.legend(loc=0, prop={'size':10})
    plt.tight_layout()
    # plt.savefig(output_figure_filename_each_caseplanet_omegan)
    plt.close()
    print("| Output figure file written to: {}".format(output_figure_filename_each_caseplanet_omegan))


    fig = plt.figure(figsize=(10,6))#*len(cases)))
    ligne = 1 #len(cases) 
    colonne = 1
    plot = 0
    for case, data_cases in zip(cases,cases_planets_data):
        c=0
        print(f'| Plot planet {key}')
        plot+=1
        ax1 = fig.add_subplot(ligne,colonne,plot)
        for key in planets_keys:
            # if key=='1':
            
            print(f'|\t Plot / label {case}')
            to_label = f' {lab[int(key)]} ' 
            planets_data = data_cases
            planet_data = planets_data[key]
            time = [t*365.25 for t in planet_data['current_time'][planet_data['current_time']<50./365.25]]
            # time = planet_data['current_time']
            substellar_point = [p*np.pi/180. for p in planet_data['substellar_point'][planet_data['current_time']<50./365.25]]
            field = 'Sub-stellar point position (rad) '
            ax1.plot(time, substellar_point, line[c], label=to_label)
            c+=1
            ax1.set_ylabel(field, fontsize=16)
            ax1.set_xlabel('time (years)', fontsize=16)
            plt.yticks(fontsize=16)
            plt.xticks(fontsize=16)
        ax1.legend(loc=0, prop={'size':10})
        ax1.axhline(y=0., ls=':', color='black')
        # ax1.axhline(y=np.pi, ls=':', color='black')
        # ax1.axhline(y=2.0*np.pi, ls=':', color='black')
        # ax1.axhline(y=4.0*np.pi, ls=':', color='black')
        # ax1.axhline(y=180., ls=':', color='black')
        # ax1.axhline(y=360., ls=':', color='black')
        # ax1.axhline(y=540., ls=':', color='black')
        # y_ticks = [0., np.pi, 2.0*np.pi, 3.0*np.pi, 4.0*np.pi, 5.0*np.pi]  # Adjust the range as needed
        # y_tick_labels = ['0', '$\pi$', '$2\pi$', '$3\pi$', '$4\pi$', '$5\pi$']
        # plt.yticks(y_ticks, y_tick_labels)
    plt.show()
    ax1.legend(loc=0, prop={'size':10})
    plt.tight_layout()
    # plt.savefig(output_figure_filename_each_case)
    plt.close()
    print("| Output figure file written to: {}".format(output_figure_filename_each_case))

if plot_ecc_data:
    fig = plt.figure(figsize=(10*len(planets_keys),6*len(cases_planets_data)))
    ligne = len(cases_planets_data)
    colonne = len(planets_keys) 
    plot = 0
    for case, data_cases in zip(cases,cases_planets_data):
        print(f'|\t Plot {case}')
        planets_data = data_cases
        
        for key in planets_keys:
            plot+=1
            ax1 = fig.add_subplot(ligne,colonne,plot)
            planet_data = planets_data[key]
            # time = [t*365.25 for t in planet_data['current_time']]
            time = planet_data['current_time']
            # substellar_point = planet_data['substellar_point']
            ecc = planet_data['eccentricity']
            long_peri = planet_data['long_per']
            field = f' {case} | {lab[int(key)]} '
            ax1.plot(time, long_peri)
            ax1.axhline(y=0., ls=':', color='black')
            ax1.set_ylabel(field, fontsize=16)
            ax1.set_xlabel('years', fontsize=16)
            plt.yticks(fontsize=16)
            plt.xticks(fontsize=16)
            ax1.legend(loc=0, prop={'size':8})
    plt.tight_layout()
    plt.savefig(output_figure_filename_all)
    plt.close()
    print("| Output figure file written to: {}".format(output_figure_filename_all))


# fig = plt.figure(figsize=(10*len(planets_keys),6*len(cases_planets_data)))
# ligne = len(cases_planets_data)
# colonne = len(planets_keys) 
# plot = 0
# for case, data_cases in zip(cases,cases_planets_data):
#     print(f'|\t Plot {case}')
#     planets_data = data_cases
#     for key in planets_keys:
#         planet_data = planets_data[key]
#         time = [t*365.25 for t in planet_data['current_time']]
#         substellar_point = planet_data['substellar_point']
#         plot+=1
#         field = f' {case} | {lab[int(key)]} '
#         ax1 = fig.add_subplot(ligne,colonne,plot)
#         ax1.plot(time, substellar_point)
#         ax1.axhline(y=0., ls=':', color='black')
#         ax1.set_ylabel(field, fontsize=16)
#         ax1.set_xlabel('days', fontsize=16)
#         # ax1.set_xlim(right=50.)
#         plt.yticks(fontsize=16)
#         plt.xticks(fontsize=16)
#         ax1.legend(loc=0, prop={'size':8})
# plt.tight_layout()
# plt.savefig(output_figure_filename_all)
# plt.close()
# print("| Output figure file written to: {}".format(output_figure_filename_all))

# fig = plt.figure(figsize=(10,6*len(planets_keys)))
# ligne = len(planets_keys) 
# colonne = 1
# plot = 0
# for key in planets_keys:
#     plot+=1
#     ax1 = fig.add_subplot(ligne,colonne,plot)
#     c=0
#     print(f'| Plot planet {key}')
#     for case, data_cases in zip(cases,cases_planets_data):
#         print(f'|\t Plot / label {case}')
#         to_label = case
#         planets_data = data_cases
#         planet_data = planets_data[key]
#         time = [t for t in planet_data['current_time']]
#         substellar_point = planet_data['substellar_point']
#         field = f' {lab[int(key)]} '
#         ax1.plot(time, substellar_point, line[c], label=to_label)
#         c+=1
#         ax1.set_ylabel(field, fontsize=16)
#         ax1.set_xlabel('days', fontsize=16)
#         plt.yticks(fontsize=16)
#         plt.xticks(fontsize=16)
#     ax1.axhline(y=0., ls=':', color='black')
#     # ax1.axhline(y=180., ls=':', color='black')
#     # ax1.axhline(y=360., ls=':', color='black')
#     # ax1.axhline(y=540., ls=':', color='black')
#     ax1.legend(loc=0, prop={'size':10})
# plt.tight_layout()
# plt.savefig(output_figure_filename_each_p)
# plt.close()
# print("| Output figure file written to: {}".format(output_figure_filename_each_p))

if plot_substellar:
    fig = plt.figure(figsize=(10,6*len(cases)))
    ligne = len(cases) 
    colonne = 1
    plot = 0
    for case, data_cases in zip(cases,cases_planets_data):
        plot+=1
        ax1 = fig.add_subplot(ligne,colonne,plot)
        c=0
        print(f'| Plot planet {key}')
        # tmp = [ 00. for i in np.arange(7) ]
        planet_period =(1.5112868036098193, 2.4227038555269917, 4.0504368087076257, 6.1029264019452052, 9.2103976604134150, 12.3562321332659266, 18.7787920964476669)
        for key in planets_keys:
            print(f'|\t Plot / label {case}')
            to_label = f' {lab[int(key)]} ' 
            planets_data = data_cases
            planet_data = planets_data[key]
            # time = [t*365.25 for t in planet_data['current_time'][planet_data['current_time']<40./365.25]]
            spin = planet_data['spin'] #[planet_data['current_time']<50./365.25]
            mean_motion = planet_data['planet_mean_motion'] #[planet_data['current_time']<50./365.25]
            omega_n = spin / mean_motion  
            print(f'|\t CASE {to_label} MIN {min(omega_n)} \t MAX {max(omega_n)} \t MED {median(omega_n)} MEAN {mean(omega_n)}')
            time = planet_data['current_time']
            substellar_point = [p*np.pi/180. for p in planet_data['substellar_point']]#[planet_data['current_time']<40./365.25]]
            field = 'Sub-stellar point position (rad) '
            ls, =ax1.plot(time, substellar_point, line[c], label=to_label)
            
            ax1.set_ylabel(field, fontsize=16)
            ax1.set_xlabel('time (years)', fontsize=16)
            plt.yticks(fontsize=16)
            plt.xticks(fontsize=16)
            # for t,s in zip(time, substellar_point):
            #     if np.abs(s) > tmp[int(key)-1]:
            #         # print("Planet: ", key, t, tmp[int(key)-1], np.abs(s) )
            #         # ax1.axhline(y=s, ls=':', color='black')
            #         ax1.axvline(x=t, ls=':', color=ls.get_color())
            #         tmp[int(key)-1]+=np.pi
            # ax1.axvline(x=planet_period[c], ls=':', color=ls.get_color())
            c+=1
        ax1.legend(loc=0, prop={'size':10})
        ax1.axhline(y=0., ls=':', color='black')
        # ax1.axhline(y=np.pi, ls=':', color='black')
        ax1.axhline(y=2.0*np.pi, ls=':', color='black')
        ax1.axhline(y=4.0*np.pi, ls=':', color='black')
        ax1.axhline(y=6.0*np.pi, ls=':', color='black')
        ax1.axhline(y=8.0*np.pi, ls=':', color='black')
        ax1.axhline(y=10.0*np.pi, ls=':', color='black')
        # ax1.axhline(y=180., ls=':', color='black')
        # ax1.axhline(y=360., ls=':', color='black')
        # ax1.axhline(y=540., ls=':', color='black')
        y_ticks = [0., 2.0*np.pi, 4.0*np.pi, 6.0*np.pi, 8.0*np.pi, 10.0*np.pi]  # Adjust the range as needed
        y_tick_labels = ['0', '$2\pi$', '$4\pi$', '$6\pi$', '$8\pi$', '$10\pi$', ]
        plt.yticks(y_ticks, y_tick_labels)
    ax1.legend(loc=0, prop={'size':10})
    plt.tight_layout()
    # plt.show()
    plt.savefig(output_figure_filename_each_case)
    plt.close()
    print("| Output figure file written to: {}".format(output_figure_filename_each_case))




# fig = plt.figure(figsize=(10,6*len(cases)))
# ligne = len(cases) 
# colonne = 1
# plot = 0
# for case, data_cases in zip(cases,cases_planets_data):
#     plot+=1
#     ax1 = fig.add_subplot(ligne,colonne,plot)
#     c=0
#     print(f'| Plot planet {key}')
#     al=1
#     for key in planets_keys:
#         print(f'|\t Plot / label {case}')
#         to_label = f' {lab[int(key)]} ' 
#         planets_data = data_cases
#         planet_data = planets_data[key]
#         time = [t*365.25 for t in planet_data['current_time'][planet_data['current_time']<150./365.25]]
#         spin = planet_data['spin'][planet_data['current_time']<150./365.25]
#         mean_motion = planet_data['planet_mean_motion'][planet_data['current_time']<150./365.25]
#         omega_n = spin / mean_motion
#         field = 'Spin state $\Omega$ /n' #case
#         ax1.plot(time, omega_n, line[c], label=to_label, alpha=al)
#         al-=0.08
#         c+=1
#         ax1.set_ylabel(field, fontsize=16)
#         ax1.set_xlabel('days', fontsize=16)
#         plt.yticks(fontsize=16)
#         plt.xticks(fontsize=16)
#         ax1.legend(loc=0, prop={'size':10})
#         ax1.axhline(y=1., ls=':', color='black')
#     # ax1.axhline(y=180., ls=':', color='black')
#     # ax1.axhline(y=360., ls=':', color='black')
#     # ax1.axhline(y=540., ls=':', color='black')
# ax1.legend(loc=0, prop={'size':10})
# plt.tight_layout()
# plt.savefig(output_figure_filename_each_case_omegan)
# plt.close()
# print("| Output figure file written to: {}".format(output_figure_filename_each_case_omegan))


# fig = plt.figure(figsize=(10,6))
# ligne = 1
# colonne = 1
# plot=0
# plot+=1
# field = 'Compare T1-b substellar'
# ax1 = fig.add_subplot(ligne,colonne,plot)
# l=0
# for case, data_cases in zip(cases,cases_planets_data):
#     print(f'|\t Plot {case}')
#     planets_data = data_cases
#     key='1'
#     lab = f' {case} | {lab[int(key)]} '
#     planet_data = planets_data[key]
#     time = [t*365.25 for t in planet_data['current_time']]
#     substellar_point = planet_data['substellar_point']
#     ax1.plot(time, substellar_point, line[l], label=lab)
#     l+=1
# ax1.legend(loc=0, prop={'size':8})
# ax1.axhline(y=0., ls=':', color='black')
# ax1.set_ylabel(field, fontsize=16)
# ax1.set_xlabel('days', fontsize=16)
# # ax1.set_xlim(right=50.)
# plt.yticks(fontsize=16)
# plt.xticks(fontsize=16)
# plt.tight_layout()
# plt.savefig(output_figure_filename_t1b)
# plt.close()
# print("| Output figure file written to: {}".format(output_figure_filename_t1b))




# fig = plt.figure(figsize=(10,6))
# ligne = 1
# colonne = 1
# plot = 0
# plot+=1
# ax1 = fig.add_subplot(ligne,colonne,plot)
# c=0

# for case, data_cases in zip(cases,cases_planets_data):
#     print(f'|\t Plot / label {case}')

#     to_label = f' {lab[int(key)]} ' 
#     planets_data = data_cases

#     mean_p1 = np.mean(planets_data['1']['planet_rotation_period'])
#     mean_p2 = np.mean(planets_data['2']['planet_rotation_period'])
#     mean_p3 = np.mean(planets_data['3']['planet_rotation_period'])
#     mean_p4 = np.mean(planets_data['4']['planet_rotation_period'])
#     mean_p5 = np.mean(planets_data['5']['planet_rotation_period'])
#     mean_p6 = np.mean(planets_data['6']['planet_rotation_period'])
#     mean_p7 = np.mean(planets_data['7']['planet_rotation_period'])

#     mean_off = mean_p3/mean_p2 
#     ls, = ax1.plot(2, mean_off, '+',markersize=20, label='Mean value')
#     ax1.plot(2,  (5./3.),  '.',markersize=20, color = ls.get_color(), label='Obs value')

#     mean_off = mean_p4/mean_p3
#     ax1.plot(3, mean_off, '+',markersize=20, )
#     ax1.plot(3, (3./2.), '.',markersize=20, color = ls.get_color())

#     mean_off = mean_p5/mean_p4
#     ax1.plot(4, mean_off, '+',markersize=20, )
#     ax1.plot(4, (3./2.), '.',markersize=20, color = ls.get_color())

#     mean_off = mean_p6/mean_p5 
#     ax1.plot(5, mean_off, '+',markersize=20,)
#     ax1.plot(5, (4./3.), '.',markersize=20, color = ls.get_color())

#     mean_off = mean_p7/mean_p6
#     ax1.plot(6, mean_off, '+',markersize=20, )
#     ax1.plot(6, (3./2.) , '.',markersize=20, color = ls.get_color())

#     plt.xticks([1,2,3,4,5,6],label=['p2/p1','p3/p2', 'p4/p3', 'p5/p4', 'p6/p5', 'p7/p6'],fontsize=10, rotation=0,)

#     field = case
#     c+=1
#     # ax1.axhline(y=0., ls=':', color='black')
#     ax1.set_ylabel(field, fontsize=16)
#     ax1.set_xlabel('days', fontsize=16)
#     plt.yticks(fontsize=16)
#     plt.xticks(fontsize=16)
#     ax1.legend(loc=0, prop={'size':8})
#     plt.grid()
# ax1.legend(loc=0, prop={'size':8})
# plt.tight_layout()
# plt.savefig(output_figure_filename_mmr)
# plt.close()
# print("| Output figure file written to: {}".format(output_figure_filename_mmr))

# // ///////////////   ///////////////////


# fig = plt.figure(figsize=(10,6))
# ligne = 1
# colonne = 1
# plot=0
# plot+=1
# field = 'Compare T1-b theta/true-ano _long-peri_corrected'
# ax1 = fig.add_subplot(ligne,colonne,plot)
# l=0
# for case, data_cases in zip(cases,cases_planets_data):
#     print(f'|\t Plot {case}')
#     planets_data = data_cases
#     key='1'
#     lab = f' {case} | {lab[int(key)]} '
#     planet_data = planets_data[key]
#     time = [t*365.25 for t in planet_data['current_time']]
#     to_plot = planet_data['theta_true_ano_long_peri']
#     ax1.plot(time, to_plot, line[l], label=lab, )
#     ax1.legend(loc=0, prop={'size':8})
#     l+=1
#     print(f'|\t\t last time {round(time[-1],3)} days')
# # ax1.axhline(y=0., ls=':', color='black')
# ax1.set_ylabel(field, fontsize=16)
# ax1.set_xlabel('days', fontsize=16)
# # ax1.set_xscale('log')
# # ax1.set_xlim(right=50.)
# plt.yticks(fontsize=16)
# plt.xticks(fontsize=16)
# plt.tight_layout()
# plt.savefig(output_figure_filename_t1b_positions_angle_true_ano)
# plt.close()
# print("| Output figure file written to: {}".format(output_figure_filename_t1b_positions_angle_true_ano))


# data['planet_mean_motion'] = planet_computed_data['mean_motion']
# data['position_x'] = planet_data['position_x']
# data['position_y'] = planet_data['position_y']
# data['position_z'] = planet_data['position_z']
# data['spin'] = planet_computed_data['omega']
# data['substellar_point'] = planet_computed_data['substellar_point']
# data['true_ano'] = planet_data['true_anomaly']
# data['mean_ano'] = planet_data['mean_anomaly']
# data['long_per'] = planet_data['longitude_of_perihelion']
# data['mean_longitude'] = planet_computed_data['mean_longitude']
# data['true_ano_cumul'] = planet_computed_data['true_anomaly_cumul']
# data['cartesian_angle_cumul'] = planet_computed_data['cartesian_position_angle']
# data['theta'] = planet_computed_data['theta']
# data['mean_ano_cumul'] = planet_computed_data['mean_ano_cumul']
# data['theta_true_ano'] = planet_computed_data['theta_true_ano'] 
# data['long_peri_rescaled'] = planet_computed_data['longitude_of_perihelion_rescaled'] 
# data['theta_true_ano_long_peri'] = planet_computed_data['theta_true_ano_long_peri']
# data['theta_mean_ano'] = planet_computed_data['theta_mean_ano'] 
# data['theta_cartesian_position_angle'] = planet_computed_data['theta_cartesian_position_angle']