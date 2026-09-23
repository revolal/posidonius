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
import sys
# from sklearn.linear_model import LinearRegression

print("> Create Timing")
# path = "/hpcstorage/bolmonte/Posidonius/20200207_TTV_TRAPPIST1"
# path_plot = "/home/spectro/bolmonte/Posidonius/tests/20200207_TTV_TRAPPIST1/plots/"

# path = "../scratch/posidonius/t1_all_earth_res_1Kyrs_ctl/"
# path_plot = "../scratch/posidonius/t1_all_earth_res_1Kyrs_ctl/"
# path = "../scratch/posidonius/t1_AS_res_allearth/t1_AS_res_allearth_history_"
# path = "../scratch/posidonius/t1_AS_res_lowfrequ/t1_AS_res_lowfrequ_history_"
#which_planet = 1
#block = 1


# Check if the command-line argument is provided
if len(sys.argv) > 1:
    # Access the argument using sys.argv[1]
    mypath = sys.argv[1]
    print(f"\t Case :{mypath}")
else:
    print("No path provided as a command-line argument.")
    mypath = "ERROR"

path = mypath

#for which_planet in np.arange(1, 8):
for which_planet in np.arange(1, 2):

    # for i in [0,1,2,3,4,5]:
    filename = path+"PLANET%i.aei" %(which_planet)
    #filename = "test%i/PLANET%i.aei" %(i, which_planet)
    print(filename)
    data = np.loadtxt(filename, skiprows = 4)

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


    filename_final_value_TTVdiff = path+"transit_timing_planet_%i.txt" %(which_planet)
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

