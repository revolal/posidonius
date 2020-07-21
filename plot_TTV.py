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
# #path = "/hpcstorage/bolmonte/Posidonius/TTV_TRAPPIST1/data_for_BOD_posidonius"
# path = "/hpcstorage/bolmonte/Posidonius/20200207_TTV_TRAPPIST1"
# path_plot = "/home/spectro/bolmonte/Posidonius/tests/20200207_TTV_TRAPPIST1/plots/"
path = "../posidonius/target"
path_plot = "../posidonius/target"

#which_planet = 1
#block = 1

masses_BO = [0, 49, 98, 147, 196, 245, 294, 343, 392, 441, 490]
masses_BO_label = ["mass_ref", "mass_3sig_neg", "mass_3sig_pos", "mass_1sig_neg", "mass_1sig_pos", "mass_05sig_neg", "mass_05sig_pos", "mass_01sig_neg", "mass_01sig_pos", "mass_001sig_neg", "mass_001sig_pos"]

#for which_mass in np.arange(0,len(masses_BO)):
for which_mass in np.arange(0,1):
    for which_planet in np.arange(1, 8):
    #for which_planet in np.arange(4, 8):
        filename_final_value_TTVdiff = "final_value_TTV_diff_planet_%i_%s_plot2_article.txt" %(which_planet, masses_BO_label[which_mass])
        f = open(filename_final_value_TTVdiff, "w+")

        #for block in np.arange(21, 22):
        # for block in np.arange(22, 23):
        for block in np.arange(23, 24):

            # vary k2Dt
            if block == 0:
                # k2 = 0
                test_compare = [0, 5, 6, 7]
                k2, k2f = 0, 0
                label = [0.01, 1.00, 100.]
            if block == 1:
                # k2 = 1
                test_compare = [0, 3, 8, 9, 10]
                k2, k2f = 1, 0
                label = [0.0, 0.01, 1.00, 100.]
            if any(x == block for x in [0,1]):
                name = "k$_{2,p}$ = %i, k$_{2f,p}$ = %i, $\Delta \\tau_p$ = " %(k2, k2f)
                filename_plot = "%sTTV_effect_dt_for_k2_%i_k2f_%i_GR_0_planet_%i_%s.pdf" %(path_plot, k2, k2f, which_planet, masses_BO_label[which_mass])


            # vary k2, k2f = 0, GR = 0
            if block == 2:
                test_compare = [0, 2, 3, 4]
                dt, k2f = 0, 0
                label = [0.01, 1.00, 100.]
            if block == 3:
                test_compare = [0, 6, 18, 9, 19, 20]
                dt, k2f = 1, 0
                label = [0.0, 0.01, 1.00, 10., 100.]
            if block == 4:
                test_compare = [0, 21, 25, 29, 33]
                dt, k2f = 1, 1
                label = [1, 10, 25, 50]
            if block == 5:
                test_compare = [0, 22, 26, 30, 34]
                dt, k2f = 1, 10
                label = [1, 10, 25, 50]
            if block == 6:
                test_compare = [0, 23, 27, 31, 35]
                dt, k2f = 1, 25
                label = [1, 10, 25, 50]
            if block == 7:
                test_compare = [0, 24, 28, 32, 36]
                dt, k2f = 1, 50
                label = [1, 10, 25, 50]
            if any(x == block for x in np.arange(2,8)):
                name = "$\Delta \\tau_p$ = %i, k$_{2f,p}$ = %i, k$_{2,p}$ =" %(dt, k2f)
                filename_plot = "%sTTV_effect_k2_for_dt_%i_k2f_%i_GR_0_planet_%i_%s.pdf" %(path_plot, dt, k2f, which_planet, masses_BO_label[which_mass])

            # vary k2f
            if block == 8:
                # k2 = 0, dt = 0
                test_compare = [0, 12, 13, 14]
                k2, dt = 0, 0
                label = [0.01, 1.00, 100.]
            if block == 9:
                test_compare = [0, 21, 22, 23, 24]
                k2, dt = 1, 1
                label = [1, 10, 25, 50]
            if block == 10:
                test_compare = [0, 25, 26, 27, 28]
                k2, dt = 10, 1
                label = [1, 10, 25, 50]
            if block == 11:
                test_compare = [0, 29, 30, 31, 32]
                k2, dt = 25, 1
                label = [1, 10, 25, 50]
            if block == 12:
                test_compare = [0, 33, 34, 35, 36]
                k2, dt = 50, 1
                label = [1, 10, 25, 50]
            if any(x == block for x in np.arange(8,13)):
                name = "$\Delta \\tau_p$ = %i, k$_{2,p}$ = %i, k$_{2f,p}$ = " %(dt, k2)
                filename_plot = "%sTTV_effect_k2f_for_k2_%i_dt_%i_GR_0_planet_%i_%s.pdf" %(path_plot, k2, dt, which_planet, masses_BO_label[which_mass])

            # general relativity
            if block == 13:
                test_compare = [0, 15, 16, 17]
                name = ""
                label = ["Kidder et al. 1995", "Anderson et al. 1975", "Newhall et al. 1983"]
                filename_plot = "%sTTV_effect_GR_for_k2_0_dt_0_k2f_0_planet_%i_%s.pdf" %(path_plot, which_planet, masses_BO_label[which_mass])


            # STELLAR k2 and ROTATIONAL FLATENNING
            # Vary k2fs
            if block == 14:
                test_compare = [0, 37, 38, 39, 40]
                k2s, dt = 1, 0
                label = [0, 1, 10, 100]
            if block == 15:
                test_compare = [0, 41, 42, 43, 44]
                k2s, dt = 10, 0
                label = [0, 1, 10, 100]
            if block == 16:
                test_compare = [0, 45, 46, 47, 48]
                k2s, dt = 100, 0
                label = [0, 1, 10, 100]
            if any(x == block for x in np.arange(14,17)):
                name = "$\Delta \\tau_\star$ = %i, k$_{2,\star}$ = %i, k$_{2f,\star}$ =" %(dt, k2s)
                filename_plot = "%sTTV_effect_k2fs_for_k2s_%i_dt_%i_planet_%i_%s.pdf" %(path_plot, k2s, dt, which_planet, masses_BO_label[which_mass])

            # Vary k2s
            if block == 17:
                test_compare = [0, 37, 41, 45]
                k2fs, dt = 0, 0
                label = [1, 10, 100]
            if block == 18:
                test_compare = [0, 38, 42, 46]
                k2fs, dt = 1, 0
                label = [1, 10, 100]
            if block == 19:
                test_compare = [0, 39, 43, 47]
                k2fs, dt = 10, 0
                label = [1, 10, 100]
            if block == 20:
                test_compare = [0, 40, 44, 48]
                k2fs, dt = 100, 0
                label = [1, 10, 100]
            if any(x == block for x in np.arange(17,21)):
                name = "$\Delta \\tau_\star$ = %i, k$_{2f,\star}$ = %i, k$_{2,\star}$ =" %(dt, k2fs)
                filename_plot = "%sTTV_effect_k2s_for_k2fs_%i_dt_%i_planet_%i_%s.pdf" %(path_plot, k2fs, dt, which_planet, masses_BO_label[which_mass])

            # All effects, for reference value
            if block == 21:
                #test_compare = [0, 6, 3, 13, 37, 1038, 15, 2000, 2001]
                #label = ["$\Delta \\tau_p$ = 1", "k$_{2,p}$ = 1", "k$_{2f,p}$ = 1", "k$_{2,\star}$ = 1", "k$_{2f,\star}$ = 1", "GR", "All effects", "All effects, k$_{2,p}$ = k$_{2f,p}$ = 50"]
                test_compare = [0, 1, 2, 3, 4, 5]
                label = ["GR" \
                         , "k$_{2,p}$ = 1 k$_{2,p}^{ref}$, k$_{2f,p}$ = 1 k$_{2f,p}^{ref}$" \
                         , "k$_{2,p}$ = 50 k$_{2,p}^{ref}$, k$_{2f,p}$ = 1 k$_{2f,p}^{ref}$" \
                         , "k$_{2,p}$ = 1 k$_{2,p}^{ref}$, k$_{2f,p}$ = 50 k$_{2f,p}^{ref}$" \
                         , "k$_{2,p}$ = 50 k$_{2,p}^{ref}$, k$_{2f,p}$ = 50 k$_{2f,p}^{ref}$"]
                name = "All effects"
                #filename_plot = "%sTTV_all_effects_planet_%i_%s.pdf" %(path_plot, which_planet, masses_BO_label[which_mass])
                filename_plot = "%sTTV_all_effects_planet_%i_%s_test.pdf" %(path_plot, which_planet, masses_BO_label[which_mass])

            if block == 22:
                test_compare = [0, 1, 2, 6, 7, 8]
                label = ["GR" \
                         , "k$_{2,p}$ = 1 k$_{2,p}^{ref}$, k$_{2f,p}$ = 1 k$_{2f,p}^{ref}$" \
                         , "k$_{2,p}$ = 10 k$_{2,p}^{ref}$, k$_{2f,p}$ = 1 k$_{2f,p}^{ref}$" \
                         , "k$_{2,p}$ = 1 k$_{2,p}^{ref}$, k$_{2f,p}$ = 10 k$_{2f,p}^{ref}$" \
                         , "k$_{2,p}$ = 10 k$_{2,p}^{ref}$, k$_{2f,p}$ = 10 k$_{2f,p}^{ref}$"]
                name = "All effects"
                #filename_plot = "%sTTV_all_effects_planet_%i_%s_test_zoom_age_%.0f.pdf" %(path_plot, which_planet, masses_BO_label[which_mass], time[0][ind_time])

            #if block == 14:
                #test_compare = []

            if block == 23:
                test_compare = [0,1,2]
                label = ["CTL"
                            , "Kaula" ]
                name = "Test No-effect CTL Kaula"
                filename_plot = "%sTTV_7_planets_%i_%s_test.jpg" %(path_plot, which_planet, masses_BO_label[which_mass])


            time = []
            x = []

            #masses_BO = [0, 49, 98]
            test_compare_mass = [test_compare[i] + masses_BO[which_mass] for i in np.arange(len(test_compare))]

            for i in test_compare_mass:
                # filename = path+"/PLANET1.aei"
                filename = path + "/test%i/PLANET%i.aei" %(i, which_planet)
                # filename = path +"/PLANET%i.aei" %( which_planet)
                print(filename)
                data = np.loadtxt (filename, skiprows = 4)

                time.append(data[:,0])
                x.append(data[:,1])

            nmax_time = len(time[0])-1
            nmid_time = int((len(time[0])-1)/2.)
            print("nmax_time = %i" %nmax_time)
            print("time[nmax_time] = %.2e" %time[0][nmax_time])
            print("nmid_time = %i" %nmid_time)
            print("time[nmid_time] = %.2e" %time[0][nmid_time])
            ind_time = nmid_time

            # if block == 22:
            #     filename_plot = "%sTTV_all_effects_planet_%i_%s_test_zoom_age_%.0f.pdf" %(path_plot, which_planet, masses_BO_label[which_mass], time[0][ind_time])

            # raise(Exception)
            # Interolate to have time at which x is 0:
            time_transit = []
            epoch = []

            for i in np.arange(0, len(test_compare_mass)):
                ep_tmp = 0
                time_transit_tmp = []
                epoch_tmp        = []
                for j in np.arange(0, len(time[i])-1):
                    if (x[i][j]<0 and x[i][j+1]>0):
                        time_transit_tmp.append(np.interp(0.0, [x[i][j], x[i][j+1]], [time[i][j], time[i][j+1]]))
                        epoch_tmp.append(ep_tmp)
                        ep_tmp = ep_tmp + 1
                time_transit.append(time_transit_tmp)
                epoch.append(epoch_tmp)

            print("len(time_transit) = ", len(time_transit[0]))

            ## PLOT X
            #fig = plt.figure(figsize=(12, 8))
            #ax = fig.add_subplot(1,1,1)

            #for i in np.arange(0, len(test_compare)):
                #line, = ax.plot(time[i], x[i], '-', label = "k$_2$ = %.0e" %label[i])
                #ax.plot(time_transit[i], np.zeros(len(time_transit[i])), 'o', c=line.get_color())

            ##ax.set_title(title_name)
            #ax.legend(loc=0, prop={'size':8})
            #ax.set_ylabel("Position x (AU)")
            #ax.set_xlabel("time (day)")
            #plt.tight_layout()
            #plt.show()

            # linear regression

            TTV = []
            TTV_diff = []

            for i in np.arange(0, len(test_compare_mass)):

                print("%i, len(time_transit) = %i" %(i, len(epoch[i])))

                #model = LinearRegression().fit(np.array(epoch[0]).reshape((-1,1)), time_transit[0])
                model = LinearRegression().fit(np.array(epoch[i]).reshape((-1,1)), time_transit[i])

                T_initial = model.intercept_
                mean_period = model.coef_
                if i == 0:
                    T_initial_0 = T_initial
                    mean_period_0 = mean_period

                TTV_tmp = []
                TTV_diff_tmp = []
                #for j in np.arange(0, int(len(epoch[i])/2.)):
                for j in np.arange(0, len(epoch[i])):

                    #print "T_initial   = ", T_initial
                    #print "mean_period = ", mean_period[0]
                    #print "T_initial   = ", float(T_initial)
                    #print "mean_period = ", float(mean_period[0])

                    tmp   = time_transit[i][j] - (T_initial   + epoch[i][j]*mean_period[0])
                    tmp_0 = time_transit[0][j] - (T_initial_0 + epoch[0][j]*mean_period_0[0])

                    TTV_tmp.append(posidonius.constants.DAY * tmp)
                    TTV_diff_tmp.append(posidonius.constants.DAY * (tmp-tmp_0))

                TTV.append(TTV_tmp)
                TTV_diff.append(TTV_diff_tmp)

            f.write("%s \n" %name)
            f.write("test, TTV_diff at 1500 day \n")
            for i in np.arange(1, len(test_compare_mass)):
                f.write("%i \t %.2e \n" %(test_compare_mass[i], TTV_diff[i][len(TTV_diff[i])-1]))


            print("PLOT HERE")

            fig = plt.figure(figsize=(8, 8))
            ax = fig.add_subplot(2,1,1)

            line, = ax.plot(epoch[0], TTV[0], '-', label = "No effects")
            for i in np.arange(1, len(test_compare_mass)):
                if (block == 13 or block == 21 or  block == 22 or block == 23):
                    line, = ax.plot(epoch[i], TTV[i], '-', label = label[i-1])
                else:
                    line, = ax.plot(epoch[i], TTV[i], '-', label = name+"%.1e" %label[i-1])


            ax.set_ylabel("TTV (s)", fontsize = 16)
            #ax.legend(loc=0, prop={'size':12})
            ax.legend(loc=0, prop={'size':16})
            ax.tick_params(axis='both', labelsize=16)

            ax = fig.add_subplot(2,1,2, sharex=ax)

            line, = ax.plot(epoch[0], [1.0*x for x in TTV_diff[0]], '-', label = "No effects")
            for i in np.arange(1, len(test_compare_mass)):
                if (block == 13 or block == 21 or block == 22 or block == 23):
                    line, = ax.plot(epoch[i], [1.0*x for x in TTV_diff[i]], '-', label = label[i-1])
                else:
                    line, = ax.plot(epoch[i], [1.0*x for x in TTV_diff[i]], '-', label = name+"%.1e" %label[i-1])
            #sum_all_effects = np.array(TTV_diff[1]) + np.array(TTV_diff[2]) + np.array(TTV_diff[3]) + np.array(TTV_diff[4]) + np.array(TTV_diff[5]) + np.array(TTV_diff[6])
            #ax.plot(epoch[0], sum_all_effects, '--', c=line.get_color())
            #sum_important_effects = np.array(TTV_diff[2]) + np.array(TTV_diff[3]) + np.array(TTV_diff[5]) + np.array(TTV_diff[6])
            #ax.plot(epoch[0], sum_important_effects, '-.', c=line.get_color())


            #ax.set_ylim([-30, 30])
            ax.tick_params(axis='both', labelsize=16)
            ax.set_ylabel("TTV difference (s)", fontsize = 16)
            ax.set_xlabel("Epoch", fontsize = 16)
            plt.tight_layout()
            #plt.show()

            print("filename plot", filename_plot)
            output_figure_filename = filename_plot
            plt.savefig(output_figure_filename, transparent=True)
            plt.close(fig)
        f.close
