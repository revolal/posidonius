import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['agg.path.chunksize'] = 10000
import matplotlib.pyplot as plt
import argparse
import posidonius
import json
import re

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('start_case_filename', action='store', help='Filename with the initial conditions of the simulation (e.g., universe_integrator.json)')
    parser.add_argument('historic_snapshot_filename', action='store', help='Filename with the historic snapshots of the simulation (e.g., universe_integrator_history.bin)')

    args = parser.parse_args()

    universe_integrator_json = json.load(open(args.start_case_filename, "r"))

    filename = args.historic_snapshot_filename
    n_particles, data = posidonius.analysis.history.read(filename)
    # star_data, planets_data, planets_keys = posidonius.analysis.history.classify(n_particles, data, discard_first_hundred_years=True)
    star_data, planets_data, planets_keys = posidonius.analysis.history.classify(n_particles, data, discard_first_years=False, first_years=1e4, discard_last_years=False, last_years=1e4 + 100. )
    star_mass = star_data['mass'][0]

    print(filename)
    # test_number_1 = re.split("test", filename)[1][0]
    # test_number_2 = re.split("test", filename)[1][1]
    # test_number_3 = re.split("test", filename)[1][2]
    # test_number_4 = re.split("test", filename)[1][3]
    # if str.isdigit(test_number_4):
    #     test_number = test_number_1+test_number_2+test_number_3+test_number_4
    # elif (str.isdigit(test_number_3) and str.isdigit(test_number_4) == False):
    #     test_number = test_number_1+test_number_2+test_number_3
    # elif (str.isdigit(test_number_3) == False and str.isdigit(test_number_2)):
    #     test_number = test_number_1+test_number_2
    # else:
    #     test_number = test_number_1


    # Data duration with the observing run of autumn
    end_time = 1050.*365.25 # day

    indice_time = len(planets_data[planets_keys[0]]['current_time'][:])- 1
    for i in np.arange(0, len(planets_data[planets_keys[0]]['current_time'])):
        if planets_data[planets_keys[0]]['current_time'][i]*365.25e0 > end_time:
            indice_time = i+1
            break

    print("planets_data[planets_keys[0]]['current_time'][i]*365.25e0 = %.2e" %(planets_data[planets_keys[0]]['current_time'][indice_time]*365.25e0))

    # path_PLANET = "/hpcstorage/bolmonte/Posidonius/20200207_TTV_TRAPPIST1/test%s/"  %(test_number)
    output_dirname = os.path.dirname(filename)
    output_filename = os.path.join(output_dirname, os.path.splitext(os.path.basename(filename))[0])
    # plt.savefig(output_figure_filename)
    path_PLANET = output_filename # "../scratch/posidonius/t1_all_earth_res_1Kyrs_ctl/"
    #path_PLANET = "/hpcstorage/bolmonte/Posidonius/20191219_TTV_TRAPPIST1/test%s/"  %("0_dt_1e-3d")
    print ("path_PLANET = ", path_PLANET)
    #raise(Exception)

    for key in planets_keys:
        planet_data = planets_data[key]

        filename_PLANET = path_PLANET+"_PLANET"+str(key)+".aei"
        #filename_PLANET = "PLANET"+str(key)+".aei"
        print(filename_PLANET)
        file = open(filename_PLANET,"w")
        file.write("\n")
        file.write("PLANET"+str(key)+"\n")
        file.write("\n")
        file.write("%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n" % ("time (days)", "x", "y", "z", "vx", "vy", "vz", "mass"))

        #for i in np.arange(0, len(planet_data['current_time'])):
        for i in np.arange(0, indice_time):
            file.write("%.13e \t %.13e \t %.13e \t %.13e \t %.13e \t %.13e \t %.13e \t %.13e \n" % (planet_data['current_time'][i]*365.25e0 \
                                                                                        , planet_data['position_x'][i], planet_data['position_y'][i], planet_data['position_z'][i] \
                                                                                        , planet_data['velocity_x'][i], planet_data['velocity_y'][i], planet_data['velocity_z'][i] \
                                                                                        , planet_data['mass'][i]))
        file.close()

        #output_text_filename = os.path.join(output_text_dirname, os.path.splitext(os.path.basename(filename_PLANET))[0])
        #all_data.to_csv(output_text_filename, sep="\t", index=False)

        #print("Output plain text file written to: {}".format(output_text_filename))



