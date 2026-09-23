import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys, os
from matplotlib import colors
import pandas as pd
import posidonius
import sklearn
from sklearn.linear_model import LinearRegression


def take_data(filename):
    
    with open(filename, 'r') as file_obj:
        # Read the first line (header)
        header = file_obj.readline()
        
        time = []
        ttv= []

        # Read the remaining lines (columns)
        for line in file_obj:
            columns = line.split()
            columns = [float(value) for value in columns if value != '']
            time.append(columns[0])
            ttv.append(columns[1])
                     
    return time, ttv



# scratch_path = "../scratch/posidonius/bolmont_ttv/AS_ci/"

# list_cases = [ "t1_all_earth_res_1Kyrs_notides", 
#          "t1_all_earth_res_1Kyrs_ctl", 
#          "t1_all_earth_res_1Kyrs", 
#          "t1_highfrequ_res_1Kyrs", 
#          "t1_lowfrequ_res_1Kyrs" ]
# list_cases = ["t1_ttv_1Kyrs_all_earth_notides", "t1_ttv_1Kyrs_all_earth_notides_test"]
# list_cases = ["t1_ttv_1Kyrs_all_earth_notides"]
# list_cases = ["t1_ttv_1Kyrs_all_earth_ctl"]
# list_cases = ["t1_ttv_1Kyrs_all_earth_ref"]
# list_cases = ["t1_ttv_1Kyrs_high_frequ"]
# list_cases = ["t1_ttv_1Kyrs_all_earth_notides", "t1_ttv_1Kyrs_all_earth_ctl", "t1_ttv_1Kyrs_high_frequ", "t1_ttv_1Kyrs_low_frequ_3", "t1_ttv_1Kyrs_all_earth_ref_3"]
# list_cases = ['t1_AS_res_allearth', 't1_AS_res_highfrequ', 't1_AS_res_lowfrequ']
# list_cases = [ 't1_ttv_AS_notides', 't1_ttv_AS_ctl', 't1_ttv_AS_allearth', 't1_ttv_AS_highfrequ', 't1_ttv_AS_lowfrequ' ]
# list_cases = [ 't1_ttv_AS_notides', 't1_ttv_AS_ctl', 't1_ttv_AS_allearth' ]

# # Case Bolmont +20 : 10yrs 
# scratch_path = "../scratch/posidonius/bolmont_ttv/AS-ci/start_0_10yrs/"
# list_cases = [ 't1_ttv_AS_notides', 't1_ttv_AS_ctl', 't1_ttv_AS_allearth', 't1_ttv_AS_highfrequ', 't1_ttv_AS_lowfrequ' ] #, "t1_ttv_AS_rotflat", "t1_ttv_AS_tide-rotflat"]
# plot_target = "./target/paper_substellar/start_0_10yrs/"

# # Case Bolmont +20 : 1500days 
# scratch_path = "../scratch/posidonius/bolmont_ttv/AS-ci/"
# list_cases = [ 't1_ttv_AS_notides', 't1_ttv_AS_ctl', 't1_ttv_AS_allearth', 't1_ttv_AS_highfrequ', 't1_ttv_AS_lowfrequ' ] #, "t1_ttv_AS_rotflat", "t1_ttv_AS_tide-rotflat"]
# plot_target = "./target/bolmont_ttv/AS_ci/"

# # Case For AS : 30000 years
# scratch_path = "../scratch/posidonius/bolmont_ttv/AS-ci/long_run/1e5yrs/"
# list_cases = [ 't1_ttv_AS_notides', 't1_ttv_AS_ctl', 't1_ttv_AS_allearth', 't1_ttv_AS_highfrequ', 't1_ttv_AS_lowfrequ' ] #, "t1_ttv_AS_rotflat", "t1_ttv_AS_tide-rotflat"]
# plot_target = "./target/bolmont_ttv/AS_ci/long_run"

# Case For AS in res : 
# scratch_path = "../scratch/posidonius/mmr_test/saved/"
# list_cases = [ "t1_AS_notides", "t1_AS_allearth", "t1_AS_allearth_O20", "t1_AS_gr" , "t1_AS_gr_allearth", "t1_AS_gr_allearth_O20"] 
# plot_target = "../scratch/posidonius/mmr_test/"
# scratch_path = "../scratch/posidonius/paper_substellar/AS_ci/start_0_10yrs/"
# list_cases = ["t1_sub_AS_rotflat_allearth", "t1_sub_AS_tides-rotflat_allearth"]

scratch_path = "/Volumes/Expansion/ygg/posidonius/substellar/thesis/"
plot_target =  "../My_thesis/tex_file/jupyter_plots/"
list_cases=[
# "t1_kaula_ctl_earth10pc"
# "t1_kaula_ctl_earth"
# "t1_ctl_earth10pc"
# "t1_ctl_earth"
"t1_maxwell_biger",
"t1_maxwell_inter",
"t1_maxwell_earth",
"t1_sundbergcooper_biger",
"t1_sundbergcooper_inter",
"t1_sundbergcooper_earth",
"t1_andrade_biger",
"t1_andrade_inter",
"t1_andrade_earth",
]

case_time_transit = {}
case_TTV = []
case_TTV_diff = []
case_epoch = []
ttv_tmp_0 = []
for index, case in enumerate(list_cases):
    dir_path = scratch_path +'/' +case+'_history'
    print('> For ', dir_path)
    # res = []
    # for path in os.listdir(dir_path):
    #     if path.endswith('_history.txt'):
    #         # print(path)
    #         res.append(path)
    # res.sort()

    time = []
    position_x = []
    for planet in np.arange(1, 8):
        filename = dir_path+"_PLANET%i.aei" %(planet)
        #filename = "test%i/PLANET%i.aei" %(i, which_planet)
        print(filename)
        data = np.loadtxt(filename, skiprows = 4)
        time.append(data[:,0])
        position_x.append(data[:,1])

    nmax_time = len(time[0])-1
    nmid_time = int((len(time[0])-1)/2.)
    print("nmax_time = %i" %nmax_time)
    print("time[nmax_time] = %.2e" %time[0][nmax_time])
    print("nmid_time = %i" %nmid_time)
    print("time[nmid_time] = %.2e" %time[0][nmid_time])
    ind_time = nmid_time

    time_transit = []
    epoch = []

    for planet in np.arange(0, 7):
        ep_tmp = 0
        time_transit_tmp = []
        epoch_tmp = []
        for j in np.arange(0, len(time[planet])-1):
            if (position_x[planet][j]<0 and position_x[planet][j+1]>0):
                time_transit_tmp.append(np.interp(0.0, [position_x[planet][j], position_x[planet][j+1]], [time[planet][j], time[planet][j+1]]))
                epoch_tmp.append(ep_tmp)
                ep_tmp = ep_tmp + 1
        # print(len(time_transit_tmp))
        time_transit.append(time_transit_tmp)
        epoch.append(epoch_tmp)

    print(f"|\t for {case} | len(time_transit) = {len(time_transit)} ")
    case_time_transit[case] = time_transit

    # # -------------------------------
    # fig = plt.figure(figsize=(18, 9))
    # ligne = 3
    # colonne = 3
    # plot = 0

    # for key in np.arange(0, 7):
    #     plot+=1
    #     ax1 = fig.add_subplot(ligne,colonne,plot,)
    #     time = [t/365.25 for t in time_transit[key]]
    #     ax1.plot(epoch[key], time , label='planet '+str(key))
    #     ax1.set_ylabel('date of transit (year)')
    #     ax1.set_xlabel('Epoch Transit Number')
    #     ax1.legend(loc=0, prop={'size':8})
    # plt.savefig(dir_path +"_planets_transit_epoch.jpg", dpi=150)
    # plt.close()

    # -------------------------------
    fig = plt.figure(figsize=(18, 9))
    ligne = 3
    colonne = 3
    plot = 0

    planets_TTV = []
    planets_TTV_diff = []
    planets_epoch = []


    for planet in np.arange(0, 7):

        print("%i, len(time_transit) = %i" %(planet, len(epoch[planet])))

        model = LinearRegression().fit(np.array(epoch[planet]).reshape((-1,1)), time_transit[planet])

        T_initial = model.intercept_
        mean_period = model.coef_

        
        if index ==0: #if case=='t1_ttv_AS_notides':
            print(' Re-compute tmp_0 in ',case)
            T_initial_0 = T_initial
            mean_period_0 = mean_period
            epoch_0 = epoch[planet]
            # time_transit_0 = time_transit[planet]
            time_transit_N0 = case_time_transit[case][planet]
            tmp = [ ( tt - (T_initial_0 + ep *mean_period_0[0])) for tt, ep in zip( time_transit_N0, epoch_0 ) ]
            # ttv_tmp_0.append(time_transit_N0 - (T_initial_0 + epoch_0[0]*mean_period_0[0]))
            ttv_tmp_0.append(tmp)
        # print('check', planet)
        # print('check', ttv_tmp_0[planet])
        print(f'\t for {case} | ttv_tmp_0 used t=0 {ttv_tmp_0[planet][0]} | t=f {ttv_tmp_0[planet][-1]} | ( {len(ttv_tmp_0[planet])} for {len(ttv_tmp_0)} planets) ')

        TTV = []
        TTV_diff = []
        cdt = True
        for t_epoch in np.arange(0, len(epoch[planet])):
            
            ttv_tmp = time_transit[planet][t_epoch] -(T_initial +epoch[planet][t_epoch] *mean_period[0])
            # tmp_0 = time_transit[0][t_epoch] - (T_initial_0 + epoch[0][t_epoch]*mean_period_0[0])
            # tmp_0 = time_transit_0[t_epoch] - (T_initial_0 + epoch_0[t_epoch]*mean_period_0[0]) 
            # tmp_0 = time_transit_N0 - (T_initial_0 + epoch_0[0]*mean_period_0[0])

            if cdt :
                # print (f"\t T_initial = {T_initial} \t| mean_period = {mean_period[0]} \t| T_initial = {float(T_initial)} " )
                print(f'Compute for {case} \t the time transit {len(time_transit[planet])} | {time_transit[planet][0]} --- {time_transit[planet][-1]}  ')
                print(f' {case} || tmp {ttv_tmp} | ttv_tmp_0 {ttv_tmp_0[planet][t_epoch]} ')
                print(f' Check TTV diff 0 = {ttv_tmp} - {ttv_tmp_0[planet][t_epoch]}  = {(ttv_tmp-ttv_tmp_0[planet][t_epoch])}')
                cdt = False
            if t_epoch == len(epoch[planet])-1:
                print(f' Check TTV diff F = {ttv_tmp} - {ttv_tmp_0[planet][t_epoch]}  = {(ttv_tmp-ttv_tmp_0[planet][t_epoch])} ')
            TTV.append(posidonius.constants.DAY *ttv_tmp)
            TTV_diff.append(posidonius.constants.DAY * (ttv_tmp-ttv_tmp_0[planet][t_epoch]))

        planets_TTV.append(TTV)
        planets_TTV_diff.append(TTV_diff)
        planets_epoch.append(epoch[planet])

        TTV_min = [t/60. for t in TTV]

        plot+=1
        ax1 = fig.add_subplot(ligne,colonne,plot,)
        ax1.plot(epoch[planet], TTV_min , label='planet '+str(planet+1))
        ax1.set_ylabel('TTV (min)')
        ax1.set_xlabel('Epoch Transit Number')
        ax1.legend(loc=0, prop={'size':8})

        # plot+=1
        # ax1 = fig.add_subplot(ligne,colonne,plot,)
        # # ax1.plot(epoch[i], TTV , label='planet '+str(i+1))
        # # ax1.plot(epoch[i], [1.0*x for x in TTV_diff[i]], '-', label = 'Planet: '+str(i+1))
        # # print(epoch[i])
        # # print(TTV_diff)
        # ax1.plot(epoch[i], TTV_diff, '-', label = 'Planet: '+str(i+1))
        # ax1.set_ylabel('TTV diff ')
        # ax1.set_xlabel('Epoch Transit Number')
        # ax1.legend(loc=0, prop={'size':8})

    case_TTV.append(planets_TTV)
    case_TTV_diff.append(planets_TTV_diff)
    case_epoch.append(planets_epoch)

    plt.savefig(dir_path +"_planets_ttv.jpg", dpi=150)
    plt.close()

# compare --::

if False:
    # -------------------------------
    fig = plt.figure(figsize=(10, 40))
    ligne = 7
    colonne = 1
    plot = 0

    for plot in np.arange(1,8):
        # plot+=1
        ax1 = fig.add_subplot(ligne,colonne,plot,)
        cc = 0
        al=1.
        for epoch, ttv in zip(case_epoch,case_TTV):
            # for planet_epoch, planet_ttv in zip(epoch, ttv):
            (planet_epoch, planet_ttv) = (epoch[plot-1], ttv[plot-1])
            planet_ttv = [t /60. for t in planet_ttv]
            if plot == 7:
                ax1.plot(planet_epoch, planet_ttv ,label=list_cases[cc], alpha=al)
                # ax1.legend(bbox_to_anchor=(1.05, 0.5))
                ax1.legend(loc=0, prop={'size':8})
            else:
                ax1.plot(planet_epoch, planet_ttv, alpha=al)
            al -=0.1
            ax1.set_ylabel('TTV (min)')
            ax1.set_xlabel('Epoch Transit Number Planet '+str(plot))
            # ax1.legend(loc=0, prop={'size':8})
            cc+=1
    plt.tight_layout()
    plt.savefig( plot_target +"compare_planets_ttv_oneplot.jpg", dpi=150)
    # plt.show()
    plt.close()


# -------------------------------
ligne = 2
colonne = 4
largeur = 6*colonne
hauteur = 4*ligne
fig = plt.figure(figsize=(largeur,hauteur))

plot = 0

for plot in np.arange(1,8):
    # plot+=1
    ax1 = fig.add_subplot(ligne,colonne,plot,)
    cc = 0
    al=1.
    for epoch, ttv in zip(case_time_transit,case_TTV):
        # for planet_epoch, planet_ttv in zip(epoch, ttv):
        # print(f"|\t  list {len(case_time_transit)} for {list_cases[cc]} take {len(case_time_transit[list_cases[cc]])} for planet {plot-1} : {len(case_time_transit[list_cases[cc]][plot-1])} ")
        (planet_epoch, planet_ttv) = (epoch[plot-1], ttv[plot-1])
        planet_epoch = case_time_transit[list_cases[cc]][plot-1]
        planet_ttv = [t /60. for t in planet_ttv]
        ax1.plot(planet_epoch, planet_ttv ,label=list_cases[cc], alpha=al)
        #     # ax1.legend(bbox_to_anchor=(1.05, 0.5))
        ax1.legend(loc=0, prop={'size':8})
        # if plot == 7:
        #     ax1.plot(planet_epoch, planet_ttv ,label=list_cases[cc], alpha=al)
        #     # ax1.legend(bbox_to_anchor=(1.05, 0.5))
        #     ax1.legend(loc=0, prop={'size':8})
        # else:
        #     ax1.plot(planet_epoch, planet_ttv, alpha=al)
        al -=0.1
        ax1.set_ylabel('TTV (min)')
        ax1.set_xlabel('Time Transit '+str(plot))
        # ax1.legend(loc=0, prop={'size':8})
        cc+=1
plt.tight_layout()
plt.savefig( plot_target +"compare_planets_ttv_vs_time_oneplot.jpg", dpi=150)
# plt.show()
plt.close()

if False:
    # -------------------------------
    fig = plt.figure(figsize=(20, 30))
    ligne = 7
    colonne = len(list_cases) #colonne = 5
    plot = 0


    for planet in np.arange(1,8):
        cc = 0
        for epoch, ttv in zip(case_epoch,case_TTV):
            plot+=1
            ax1 = fig.add_subplot(ligne,colonne,plot,)

            # for planet_epoch, planet_ttv in zip(epoch, ttv):
            (planet_epoch, planet_ttv) = (epoch[planet-1], ttv[planet-1])
            if plot == 7:
                ax1.plot(planet_epoch, planet_ttv ,label=list_cases[cc])
                ax1.legend(bbox_to_anchor=(1.05, 0.5))
            else:
                ax1.plot(planet_epoch, planet_ttv, label=list_cases[cc])
            ax1.set_ylabel('TTV ')
            ax1.set_xlabel('Epoch Transit Number Planet '+str(plot))
            # ax1.legend(loc=0, prop={'size':8})
            cc+=1
    plt.tight_layout()
    plt.savefig( plot_target +"compare_planets_ttv.jpg", dpi=150)
    plt.close()

if False:
    # -------------------------------
    fig = plt.figure(figsize=(10, 40))
    ligne = 7
    colonne = 1
    plot = 0

    for plot in np.arange(1,8):
        # plot+=1
        ax1 = fig.add_subplot(ligne,colonne,plot,)
        cc = 0
        al=1.
        for epoch, ttv_diff in zip(case_epoch,case_TTV_diff):
            # for planet_epoch, planet_ttv in zip(epoch, ttv):
            (planet_epoch, planet_ttv_diff) = (epoch[plot-1], ttv_diff[plot-1])
            if plot == 7:
                ax1.plot(planet_epoch, planet_ttv_diff ,label=list_cases[cc], alpha=al)
                # ax1.legend(bbox_to_anchor=(1.05, 0.5))
                ax1.legend(loc=0, prop={'size':8})
            else:
                ax1.plot(planet_epoch, planet_ttv_diff, alpha=al)
            ax1.set_ylabel('TTV Diff (s)')
            ax1.set_xlabel('Epoch Transit Number Planet '+str(plot))
            # ax1.legend(loc=0, prop={'size':8})
            cc+=1
            al-=0.1
    plt.tight_layout()
    plt.savefig( plot_target +"compare_planets_ttv_diff_oneplot.jpg", dpi=150)
    plt.close()
# -------------------------------

ligne = 2
colonne = 4
largeur = 6*colonne
hauteur = 4*ligne
fig = plt.figure(figsize=(largeur,hauteur))

plot = 0

for plot in np.arange(1,8):
    # plot+=1
    ax1 = fig.add_subplot(ligne,colonne,plot,)
    cc = 0
    al=1.
    for epoch, ttv_diff in zip(case_epoch,case_TTV_diff):
        # for planet_epoch, planet_ttv in zip(epoch, ttv):
        (planet_epoch, planet_ttv_diff) = (epoch[plot-1], ttv_diff[plot-1])
        planet_epoch = case_time_transit[list_cases[cc]][plot-1]
        ax1.plot(planet_epoch, planet_ttv_diff ,label=list_cases[cc], alpha=al)
            # ax1.legend(bbox_to_anchor=(1.05, 0.5))
        ax1.legend(loc=0, prop={'size':8})
        # if plot == 7:
        #     ax1.plot(planet_epoch, planet_ttv_diff ,label=list_cases[cc], alpha=al)
        #     # ax1.legend(bbox_to_anchor=(1.05, 0.5))
        #     ax1.legend(loc=0, prop={'size':8})
        # else:
        #     ax1.plot(planet_epoch, planet_ttv_diff, alpha=al)
        ax1.set_ylabel('TTV Diff (s)')
        ax1.set_xlabel('Epoch Transit Number Planet '+str(plot))
        # ax1.legend(loc=0, prop={'size':8})
        cc+=1
        al-=0.1
plt.tight_layout()
plt.savefig( plot_target +"compare_planets_ttv_diff_vs_time_oneplot.jpg", dpi=150)
# plt.show()
plt.close()

if False:
    # -------------------------------
    fig = plt.figure(figsize=(20, 30))
    ligne = 7
    colonne = len(list_cases) #colonne = 5
    plot = 0

    for planet in np.arange(1,8):
        cc = 0
        for epoch, ttv_diff in zip(case_epoch,case_TTV_diff):
            plot+=1
            ax1 = fig.add_subplot(ligne,colonne,plot,)
        
            (planet_epoch, planet_ttv_diff) = (epoch[planet-1], ttv_diff[planet-1])

            # for planet_epoch, planet_ttv in zip(epoch, ttv):
            if plot == 7:
                ax1.plot(planet_epoch, planet_ttv_diff ,label=list_cases[cc])
                # ax1.legend(bbox_to_anchor=(1.05, 0.5))
            else:
                ax1.plot(planet_epoch, planet_ttv_diff, label=list_cases[cc])
            ax1.set_ylabel('TTV Diff ')
            ax1.set_xlabel('Epoch Transit Number Planet '+str(planet) +' ' +list_cases[cc])
            ax1.legend(loc=0, prop={'size':8})
            cc+=1
    plt.tight_layout()
    plt.savefig( plot_target +"compare_planets_ttv_diff.jpg", dpi=150)
    plt.close()

