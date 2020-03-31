import posidonius
import numpy as np
import argparse

from scipy import stats
import matplotlib.pyplot as plt
from pylab import plot, show, errorbar


Trappist_1e_1 = np.loadtxt('Results_Trappist_1e_00.0_138_freq_Imk2_posidonius.txt',comments='#')
Trappist_1e_2 = np.loadtxt('Results_Trappist_1e_05.0_150_freq_Imk2_posidonius.txt',comments='#')

print("The first one \n ", Trappist_1e_1)
# print("The second one \n", Trappist_1e_2)

print("OK then\n", Trappist_1e_1[0,:])


Mass1, Radius1, Gyration_radius1 = Trappist_1e_1[0,:] 
print("\nMass1 ", Mass1, "\nRad ", Radius1, "\nGyration ", Gyration_radius1)
w_lmpq_1 = Trappist_1e_1 [1:,0]
ImK2_1 = Trappist_1e_1 [1:,1]
ReK2_1 = Trappist_1e_1 [1:,2]

Mass2, Radius2, Gyration_radius2 = Trappist_1e_2[0,:] 
print("\nMass1 ", Mass2, "\nRad ", Radius2, "\nGyration ", Gyration_radius2)
w_lmpq_2 = Trappist_1e_2 [1:,0]
ImK2_2 = Trappist_1e_2 [1:,1]
ReK2_2 = Trappist_1e_2 [1:,2]


# slope_Im1, intercept_Im1, r_value_Im1, p_value_Im1, std_err_Im1 = stats.linregress(ImK2_1, w_lmpq_1)
# print("\nslope_Im1 %f     intercept_Im1 %f        r_value_Im1 %f      p_value_Im1 %f       std_err_Im1 %f" %(slope_Im1, intercept_Im1, r_value_Im1, p_value_Im1, std_err_Im1) )

# slope_Re1, intercept_Re1, r_value_Re1, p_value_Re1, std_err_Re1 = stats.linregress(ReK2_1, w_lmpq_1)
# print("\nslope_Re1 %f       intercept_Re1 %f        r_value_Re1 %f      p_value_Re1 %f      std_err_Re1 %f" % (slope_Re1, intercept_Re1, r_value_Re1, p_value_Re1, std_err_Re1))


plt.plot(w_lmpq_1, ReK2_1, 'o', label='original data', color='red')
plt.plot(w_lmpq_1, ImK2_1, 'o', label='original data', color='red')

plt.plot(w_lmpq_2, ReK2_2, 'o', label='original data')
plt.plot(w_lmpq_2, ImK2_2, 'o', label='original data')

plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.show()

def easy_regress(X,Y):
    return