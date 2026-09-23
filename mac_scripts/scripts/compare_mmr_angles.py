import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys, os
from matplotlib import colors
import pandas as pd
from statistics import mean, median
from scipy.interpolate import interp1d



# Plot Sub-stellar comparison
# pwd="/Volumes/Expansion/Mac_data/posidonius/from_baobab/substellar_internal"
pwd="/Volumes/Expansion/baobab/posidonius"
data_dir_path =pwd+"/mmr_test/"

cases=(
# "t1_AS-t1a-mass-l-sma_res_allearth" ,
# "t1_AS-t1a-mass-l_res_allearth" ,
# "t1_AS-t1a-mass_res_allearth" ,
"t1_AS_allearth" ,
# "t1_AS_allearth_O20" ,
# "t1_AS_gr" ,
# "t1_AS_gr_allearth" ,
# "t1_AS_gr_allearth_O20" ,
# "t1_AS_notides",
)




###########################################################
# # # Have to loop over the path and cases:
data_dir_path_cases = [ data_dir_path +p +'/' for p in cases]
print('| Cases loaded > ', data_dir_path_cases)
res = []
for data_dir_path_case in data_dir_path_cases:
    for path in os.listdir(data_dir_path_case):
        if path.endswith('history_save.txt'):
            print('|\t find > ', path, ' in ', data_dir_path_case)
            res.append(path)
# res.sort()
print('| files loaded > ', res)
cases_planets_data = []
cases = []
for path,file in zip(data_dir_path_cases,res):
    data = path+file
    print('|\t Look for ', data)
    df = pd.read_csv(data, sep='\t', header=0, )
    case = file[11:-17]
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
print('| N cases loaded : ', len(cases_planets_data))



# ###########################################################
# data_dir_path_cases = data_dir_path 
# print('| Cases loaded > ', data_dir_path_cases)
# res = []
# for path in os.listdir(data_dir_path_cases):
#     print(f'| Evaluate ', path)
#     if path.endswith('.txt'):
#         # if not path.endswith('bis_history_save.txt'):
#         if path.endswith('_history_save.txt'):
#         # if path.endswith('intermediate_history_save.txt'):
#             print('|\t find > ', path)
#             res.append(path)
# res.sort(reverse=True)
# print(f'| Find {len(res)} cases in ', res)

# cases_planets_data = []
# cases = []
# for file in (res):
#     data = data_dir_path_cases+file
#     df = pd.read_csv(data, sep='\t', header=0, )
#     case = file[:-17]
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

# first_years = 0.
# last_years = 0. #50/365.25
# tmp = []
# for data_cases in cases_planets_data:
#     planets_data = data_cases
#     for key in planets_keys:
#         # print(f'|\t cut time for {key}')
#         planet_data = planets_data[key]
#         if first_years >0.:
#             planet_data = planet_data[planet_data['current_time'] >= first_years]
#         if last_years >0.:
#             planet_data = planet_data[planet_data['current_time'] <= last_years]
#         planets_data[key] = planet_data
#     data_cases = planets_data
#     tmp.append(data_cases)
# cases_planets_data = tmp

# print('| N cases loaded : ', len(cases_planets_data))

# ###########################################################


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


plot_mmr_angles = True


if plot_mmr_angles:
    
    if True : # Resonance 3 bodies

        for case, data_cases in zip(cases,cases_planets_data):
            planets_data = data_cases

            fig = plt.figure(figsize=(8, 20))
            ligne = 5
            colonne = 1
            i = 0

            # plot+=1
            # ax1 = fig.add_subplot(ligne,colonne,plot)
            # for key in planets_keys:

            # t1 = planets_data['1']['current_time']
            # t2 = planets_data['2']['current_time']
            # t3 = planets_data['3']['current_time']
            # t4 = planets_data['4']['current_time']
            # t5 = planets_data['5']['current_time']
            # t6 = planets_data['6']['current_time']
            t7  = planets_data['7']['current_time']

            time = t7[t7<1600./365.25]

            for key in planets_keys:
                planet_data = planets_data[key]
                time_filter = 1600. /365.25
                tmp = planet_data['mean_longitude']
                planet_data['mean_longitude'] = tmp[t7<time_filter]
                tmp = planet_data['long_per']
                planet_data['long_per'] = tmp[t7<time_filter]
                planets_data[key] = planet_data


            l_1 = planets_data['1']['mean_longitude']
            l_2 = planets_data['2']['mean_longitude']
            l_3 = planets_data['3']['mean_longitude']
            l_4 = planets_data['4']['mean_longitude']
            l_5 = planets_data['5']['mean_longitude']
            l_6 = planets_data['6']['mean_longitude']
            l_7 = planets_data['7']['mean_longitude']

            omega_1 =  planets_data['1']['long_per']
            omega_2 =  planets_data['2']['long_per']
            omega_3 =  planets_data['3']['long_per']
            omega_4 =  planets_data['4']['long_per']
            omega_5 =  planets_data['5']['long_per']
            omega_6 =  planets_data['6']['long_per']
            omega_7 =  planets_data['7']['long_per']

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
            field = 'Res Angle b c d'
            ax1.plot(time, theta_123, '+', label=f'angle 1 init {round(theta_123[0],3)}')
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Res Angle c d e'
            ax1.plot(time, theta_234, '+', label=f'angle 2 init {round(theta_234[0],3)}')
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Res Angle d e f'
            ax1.plot(time, theta_345, '+', label=f'angle 3 init: {round(theta_345[0],3)}')
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Res Angle e f g'
            ax1.plot(time, theta_456, '+', label=f'angle 4 init: {round(theta_456[0],3)}')
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Res Angle f g h '
            ax1.plot(time, theta_567, '+', label= f'angle 5 init: {round(theta_567[0],3)}')
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

            #-----------------------------------------------------------------------------------
            plt.tight_layout()
            filename = 'target/mmr/' +case +'_2_body_res_angles.png'
            plt.savefig(filename)
            # plt.show()
            plt.close()
            print("> Output figure file written to: {}".format(filename))

    if True: # Resonance 2 bodies

        for case, data_cases in zip(cases,cases_planets_data):
            planets_data = data_cases

            #-----------------------------------------------------------------------------------

            fig = plt.figure(figsize=(8, 20))
            ligne = 5
            colonne = 1
            plot = 0

            plot+=1
            ax1 = fig.add_subplot(ligne,colonne,plot)

            # t1 = planets_data['1']['current_time']
            # t2 = planets_data['2']['current_time']
            # t3 = planets_data['3']['current_time']
            # t4 = planets_data['4']['current_time']
            # t5 = planets_data['5']['current_time']
            # t6 = planets_data['6']['current_time']
            t7  = planets_data['7']['current_time']

            time = t7[t7<1600./365.25]

            l_b = planets_data['1']['mean_longitude'][t7<1600./365.25]
            l_c = planets_data['2']['mean_longitude'][t7<1600./365.25]
            l_d = planets_data['3']['mean_longitude'][t7<1600./365.25]
            l_e = planets_data['4']['mean_longitude'][t7<1600./365.25]
            l_f = planets_data['5']['mean_longitude'][t7<1600./365.25]
            l_g = planets_data['6']['mean_longitude'][t7<1600./365.25]
            l_h = planets_data['7']['mean_longitude'][t7<1600./365.25]

            omega_b =  planets_data['1']['long_per'][t7<1600./365.25]
            omega_c =  planets_data['2']['long_per'][t7<1600./365.25]
            omega_d =  planets_data['3']['long_per'][t7<1600./365.25]
            omega_e =  planets_data['4']['long_per'][t7<1600./365.25]
            omega_f =  planets_data['5']['long_per'][t7<1600./365.25]
            omega_g =  planets_data['6']['long_per'][t7<1600./365.25]
            omega_h =  planets_data['7']['long_per'][t7<1600./365.25]

            # 8:5 | 4 angles b c
            theta1_85_bc = 8.0*l_c -5.0*l_b -3.0*omega_b
            theta2_85_bc = 8.0*l_c -5.0*l_b -3.0*omega_c
            theta3_85_bc = 8.0*l_c -5.0*l_b -2.0*omega_b -omega_c
            theta4_85_bc = 8.0*l_c -5.0*l_b -omega_b -2.0*omega_c
            theta1_85_bc = theta1_85_bc%(2.0*np.pi)
            theta2_85_bc = theta2_85_bc%(2.0*np.pi)
            theta3_85_bc = theta3_85_bc%(2.0*np.pi)
            theta4_85_bc = theta4_85_bc%(2.0*np.pi)
            theta1_85_bc = (theta1_85_bc) *180./np.pi
            theta2_85_bc = (theta2_85_bc) *180./np.pi
            theta3_85_bc = (theta3_85_bc) *180./np.pi
            theta4_85_bc = (theta4_85_bc) *180./np.pi

            # 5:3 | 3 angles c d 
            theta1_53_cd = 5.0*l_d -3.0*l_c -2.0*omega_c
            theta2_53_cd = 5.0*l_d -3.0*l_c -2.0*omega_d
            theta3_53_cd = 5.0*l_d -3.0*l_c -omega_c -omega_d
            theta1_53_cd = theta1_53_cd%(2.0*np.pi)
            theta2_53_cd = theta2_53_cd%(2.0*np.pi)
            theta3_53_cd = theta3_53_cd%(2.0*np.pi)
            theta1_53_cd = (np.pi-theta1_53_cd)*180./np.pi
            theta2_53_cd = (np.pi-theta2_53_cd)*180./np.pi
            theta3_53_cd = (np.pi-theta3_53_cd)*180./np.pi

            # 3:2 | 2 angles d e
            theta1_32_de = 3.0*l_e -2.0*l_d -omega_d
            theta2_32_de = 3.0*l_e -2.0*l_d -omega_e
            theta1_32_de = theta1_32_de%(2.0*np.pi)
            theta2_32_de = theta2_32_de%(2.0*np.pi)
            theta1_32_de = (np.pi-theta1_32_de)*180./np.pi
            theta2_32_de = (np.pi-theta2_32_de)*180./np.pi

            # 3:2 | 2 angles e f
            theta1_32_ef = 3.0*l_f -2.0*l_e -omega_e
            theta2_32_ef = 3.0*l_f -2.0*l_e -omega_f
            theta1_32_ef = theta1_32_ef%(2.0*np.pi)
            theta2_32_ef = theta2_32_ef%(2.0*np.pi)
            theta1_32_ef = (np.pi-theta1_32_ef)*180./np.pi
            theta2_32_ef = (np.pi-theta2_32_ef)*180./np.pi

            # 4:3 | 2 angles f g 
            theta1_43_fg = 4.0*l_g -3.0*l_f -omega_f
            theta2_43_fg = 4.0*l_g -3.0*l_f -omega_g
            theta1_43_fg = theta1_43_fg%(2.0*np.pi)
            theta2_43_fg = theta2_43_fg%(2.0*np.pi)
            theta1_43_fg = (np.pi-theta1_43_fg)*180./np.pi
            theta2_43_fg = (np.pi-theta2_43_fg)*180./np.pi

            #3:2 | 2 angles g h
            theta1_32_gh = 3.0*l_h -2.0*l_g -omega_g
            theta2_32_gh = 3.0*l_h -2.0*l_g -omega_h
            theta1_32_gh = theta1_32_gh%(2.0*np.pi)
            theta2_32_gh = theta2_32_gh%(2.0*np.pi)
            theta1_32_gh = (np.pi-theta1_32_gh)*180./np.pi
            theta2_32_gh = (np.pi-theta2_32_gh)*180./np.pi
            
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
            field = 'Res Angle 8:5 | b c angle 1'
            ax1.plot(t7, theta1_85_bc, '+', label=f'angle 1 init: {round(theta1_85_bc[0],3)}')
            ax1.plot(t7, theta2_85_bc, '+', label=f'angle 2 init: {round(theta2_85_bc[0],3)}')
            ax1.plot(t7, theta3_85_bc, '+', label=f'angle 3 init: {round(theta3_85_bc[0],3)}')
            ax1.plot(t7, theta4_85_bc, '+', label=f'angle 4 init: {round(theta4_85_bc[0],3)}')
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Res Angle 5:3 | c d angle 1'
            ax1.plot(t7, theta1_53_cd, '+', label=f'angle 1 init: {round(theta1_53_cd[0],3)}')
            ax1.plot(t7, theta2_53_cd, '+', label=f'angle 2 init: {round(theta2_53_cd[0],3)}', alpha=1.)
            ax1.plot(t7, theta3_53_cd, '+', label=f'angle 3 init: {round(theta3_53_cd[0],3)}', alpha=1.)
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})


            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Res Angle 3:2 | d e '
            ax1.plot(t7, theta1_32_de, '+', label=f'angle 1 init: {round(theta1_32_de[0],3)}')
            ax1.plot(t7, theta2_32_de, '+', label=f'angle 2 init: {round(theta2_32_de[0],3)}', alpha=.7)
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Res Angle 3:2 | e f '
            ax1.plot(t7, theta1_32_ef, '+', label=f'angle 1 init: {round(theta1_32_ef[0],3)}')
            ax1.plot(t7, theta2_32_ef, '+', label=f'angle 2 init: {round(theta2_32_ef[0],3)}', alpha=.7)
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Res Angle 4:3 | f g '
            ax1.plot(t7, theta1_43_fg, '+', label=f'angle 1 init: {round(theta1_43_fg[0],3)}')
            ax1.plot(t7, theta2_43_fg, '+', label=f'angle 2 init: {round(theta2_43_fg[0],3)}', alpha=.7)
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

            i=i+1
            ax1 = fig.add_subplot(ligne,colonne,i,)
            field = 'Res Angle 3:2 | g h '
            ax1.plot(t7, theta1_32_gh, '+', label=f'angle 1 init: {round(theta1_32_gh[0],3)}')
            ax1.plot(t7, theta2_32_gh, '+', label=f'angle 2 init: {round(theta2_32_gh[0],3)}', alpha=.7)
            ax1.set_ylabel(field+'  ')
            ax1.set_xlabel('years')
            ax1.legend(loc=0, prop={'size':8})

            #-----------------------------------------------------------------------------------
            plt.tight_layout()
            filename = 'target/mmr/' +case +'_3_body_res_angles.png'
            plt.savefig(filename)
            # plt.show()
            plt.close()
            print("> Output figure file written to: {}".format(filename))
