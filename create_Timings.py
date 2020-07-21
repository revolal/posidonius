import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
#mpl.rcParams['agg.path.chunksize'] = 10000
import matplotlib.pyplot as plt
import argparse
import posidonius
import json
from sklearn.linear_model import LinearRegression

print("")
# path = "/hpcstorage/bolmonte/Posidonius/20200207_TTV_TRAPPIST1"
# path_plot = "/home/spectro/bolmonte/Posidonius/tests/20200207_TTV_TRAPPIST1/plots/"
path = "../posidonius/target"
path_plot = "../posidonius/target"
#which_planet = 1
#block = 1


# filename = path+"/PLANET1.aei"
# #filename = "test%i/PLANET%i.aei" %(i, which_planet)
# print(filename)
# data = np.loadtxt (filename, skiprows = 4)

# time = data[:,0]
# x    = data[:,1]

# # Interolate to have time at which x is 0:
# time_transit = []
# epoch = []

# ep_tmp = 0
# time_transit_tmp = []
# epoch_tmp        = []
# for j in np.arange(0, len(time)-1):
#     if (x[j]<0 and x[j+1]>0):
#         time_transit_tmp.append(np.interp(0.0, [x[j], x[j+1]], [time[j], time[j+1]]))
#         epoch_tmp.append(ep_tmp)
#         ep_tmp = ep_tmp + 1
# time_transit = time_transit_tmp
# epoch = epoch_tmp


# filename_final_value_TTVdiff = "../posidonius/target/data_transit_timings_T1/transit_timing_planet_1_test_1.txt"
# f = open(filename_final_value_TTVdiff, "w+")

# for j in np.arange(0, len(epoch)):
#     f.write("%.0f \t %.20e \n" %(epoch[j], time_transit[j]))
# f.close


# #print("PLOT HERE")

# #fig = plt.figure(figsize=(8, 8))
# #ax = fig.add_subplot(1,1,1)

# #line, = ax.plot(epoch, time_transit, '-', label = "No noise")
# #for j in np.arange(0, len(noise_level)):
#     #i = len(noise_level)-1 -j
#     #ax.plot(epoch, noisy_time_transit[i], '.', label = "noise level = %.1e s" %noise_level[i])
# #ax.plot(epoch, time_transit, '-', c=line.get_color())


# #ax.set_ylabel("Transit Time (s)", fontsize = 16)
# #ax.legend(loc=0, prop={'size':12})
# #ax.tick_params(axis='both', labelsize=14)

# #plt.tight_layout()
# ##plt.show()

# #filename_plot = "plots/Noisy_transit_time_k2_50_dt_1_k2f_1_planet_%i_mass_ref.pdf" %which_planet
# #print(filename_plot)
# #output_figure_filename = filename_plot
# #plt.savefig(output_figure_filename, transparent=True)
# #plt.close(fig)





#for which_planet in np.arange(1, 8):
for which_planet in np.arange(1, 8):

    for i in [0,1,2]:
        filename = path+"/test%i/PLANET%i.aei"  %(i, which_planet)
        #filename = "test%i/PLANET%i.aei" %(i, which_planet)
        print(filename)
        data = np.loadtxt (filename, skiprows = 4)

        time = data[:,0]
        x    = data[:,1]

        # Interolate to have time at which x is 0:
        time_transit = []
        epoch = []

        ep_tmp = 0
        time_transit_tmp = []
        epoch_tmp        = []
        for j in np.arange(0, len(time)-1):
            if (x[j]<0 and x[j+1]>0):
                time_transit_tmp.append(np.interp(0.0, [x[j], x[j+1]], [time[j], time[j+1]]))
                epoch_tmp.append(ep_tmp)
                ep_tmp = ep_tmp + 1
        time_transit = time_transit_tmp
        epoch = epoch_tmp


        # filename_final_value_TTVdiff = "../posidonius/target/data_transit_timings_T1/transit_timing_planet_%i_test_.txt" %(i)
        filename_final_value_TTVdiff = "../posidonius/target/data_transit_timings_T1/transit_timing_planet_%i_test_%i.txt" %(which_planet, i)
        f = open(filename_final_value_TTVdiff, "w+")

        for j in np.arange(0, len(epoch)):
            f.write("%.0f \t %.20e \n" %(epoch[j], time_transit[j]))
        f.close


        #print("PLOT HERE")

        #fig = plt.figure(figsize=(8, 8))
        #ax = fig.add_subplot(1,1,1)

        #line, = ax.plot(epoch, time_transit, '-', label = "No noise")
        #for j in np.arange(0, len(noise_level)):
            #i = len(noise_level)-1 -j
            #ax.plot(epoch, noisy_time_transit[i], '.', label = "noise level = %.1e s" %noise_level[i])
        #ax.plot(epoch, time_transit, '-', c=line.get_color())


        #ax.set_ylabel("Transit Time (s)", fontsize = 16)
        #ax.legend(loc=0, prop={'size':12})
        #ax.tick_params(axis='both', labelsize=14)

        #plt.tight_layout()
        ##plt.show()

        #filename_plot = "plots/Noisy_transit_time_k2_50_dt_1_k2f_1_planet_%i_mass_ref.pdf" %which_planet
        #print(filename_plot)
        #output_figure_filename = filename_plot
        #plt.savefig(output_figure_filename, transparent=True)
        #plt.close(fig)

