import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.interpolate
import sys, os

Rsun=695800e3 # (m)
Msun=1.99e30 # (kg)
AU=149597870700 # (m)
Grav_const = 6.67E-11 # Gravitational constant, SI units
Pi = 3.14159265359 # The number Pi




#==================================================================================================================
print("Recover K2 data")



#===============================================================================================
planet_data_SE1 = np.loadtxt('Results_Model_Super_Earth_SOL_res05_00_100_freq_Imk2_posidonius.txt',comments='#')
planet_mass_SE1, planet_radius_SE1, planet_gyration_radius_squared_SE1 = planet_data_SE1[0,:] 

w_lm_SE1   = planet_data_SE1[1:,0]
ImK2_SE1   = planet_data_SE1[1:,1]
ReK2_SE1   = planet_data_SE1[1:,2]
size_SE1   = np.size(w_lm_SE1)


#===============================================================================================
planet_data_SE2 = np.loadtxt('Results_Model_Super_Earth_SOL_res05_00_150_freq_Imk2_posidonius.txt',comments='#')
planet_mass_SE2, planet_radius_SE2, planet_gyration_radius_squared_SE2 = planet_data_SE2[0,:] 

w_lm_SE2   = planet_data_SE2[1:,0]
ImK2_SE2   = planet_data_SE2[1:,1]
ReK2_SE2   = planet_data_SE2[1:,2]
size_SE2   = np.size(w_lm_SE2)

#===============================================================================================
planet_data_SE3 = np.loadtxt('Results_Model_Super_Earth_SOL_res10_00_050_freq_Imk2_posidonius.txt',comments='#')
planet_mass_SE3, planet_radius_SE3, planet_gyration_radius_squared_SE3 = planet_data_SE3[0,:] 

w_lm_SE3   = planet_data_SE3[1:,0]
ImK2_SE3   = planet_data_SE3[1:,1]
ReK2_SE3   = planet_data_SE3[1:,2]
size_SE3   = np.size(w_lm_SE3)

#===============================================================================================
planet_data_SE4 = np.loadtxt('Results_Model_Super_Earth_SOL_res10_00_100_freq_Imk2_posidonius.txt',comments='#')
planet_mass_SE4, planet_radius_SE4, planet_gyration_radius_squared_SE4 = planet_data_SE4[0,:] 

w_lm_SE4   = planet_data_SE4[1:,0]
ImK2_SE4   = planet_data_SE4[1:,1]
ReK2_SE4   = planet_data_SE4[1:,2]
size_SE4   = np.size(w_lm_SE4)


#===============================================================================================
planet_data_b = np.loadtxt('Results_Trappist_1b_20.5_100_freq_Imk2_posidonius.txt',comments='#')
# planet_data_b = np.loadtxt('Results_Trappist_1b_24.0_150_freq_Imk2_posidonius.txt',comments='#')
planet_mass_b, planet_radius_b, planet_gyration_radius_squared_b = planet_data_b[0,:] 

w_lm_b   = planet_data_b[1:,0]
ImK2_b   = planet_data_b[1:,1]
ReK2_b   = planet_data_b[1:,2]
size_b   = np.size(w_lm_b)

#===============================================================================================
planet_data_c = np.loadtxt('Results_Trappist_1c_07.5_100_freq_Imk2_posidonius.txt',comments='#')
# planet_data_c = np.loadtxt('Results_Trappist_1c_11.0_150_freq_Imk2_posidonius.txt',comments='#')
planet_mass_c, planet_radius_c, planet_gyration_radius_squared_c = planet_data_c[0,:] 

w_lm_c   = planet_data_c[1:,0]
ImK2_c   = planet_data_c[1:,1]
ReK2_c   = planet_data_c[1:,2]
size_c   = np.size(w_lm_c)

#===============================================================================================
planet_data_d = np.loadtxt('Results_Trappist_1d_19.5_100_freq_Imk2_posidonius.txt',comments='#')
# planet_data_d = np.loadtxt('Results_Trappist_1d_22.5_150_freq_Imk2_posidonius.txt',comments='#')
planet_mass_d, planet_radius_d, planet_gyration_radius_squared_d = planet_data_d[0,:] 

w_lm_d   = planet_data_d[1:,0]
ImK2_d   = planet_data_d[1:,1]
ReK2_d   = planet_data_d[1:,2]
size_d   = np.size(w_lm_d)

#===============================================================================================
planet_data_e = np.loadtxt('Results_Trappist_1e_00.0_138_freq_Imk2_posidonius.txt',comments='#')
# planet_data_e = np.loadtxt('Results_Trappist_1e_05.0_150_freq_Imk2_posidonius.txt',comments='#')
planet_mass_e, planet_radius_e, planet_gyration_radius_squared_e = planet_data_e[0,:] 

w_lm_e   = planet_data_e[1:,0]
ImK2_e   = planet_data_e[1:,1]
ReK2_e   = planet_data_e[1:,2]
size_e   = np.size(w_lm_e)

#===============================================================================================
planet_data_f = np.loadtxt('Results_Trappist_1f_11.0_100_freq_Imk2_posidonius.txt',comments='#')
# planet_data_f = np.loadtxt('Results_Trappist_1f_14.5_150_freq_Imk2_posidonius.txt',comments='#')
planet_mass_f, planet_radius_f, planet_gyration_radius_squared_f = planet_data_f[0,:] 

w_lm_f   = planet_data_f[1:,0]
ImK2_f   = planet_data_f[1:,1]
ReK2_f   = planet_data_f[1:,2]
size_f   = np.size(w_lm_f)

#===============================================================================================
planet_data_g = np.loadtxt('Results_Trappist_1g_19.0_100_freq_Imk2_posidonius.txt',comments='#')
# planet_data_g = np.loadtxt('Results_Trappist_1g_22.5_150_freq_Imk2_posidonius.txt',comments='#')
planet_mass_g, planet_radius_g, planet_gyration_radius_squared_g = planet_data_g[0,:] 

w_lm_g   = planet_data_g[1:,0]
ImK2_g   = planet_data_g[1:,1]
ReK2_g   = planet_data_g[1:,2]
size_g   = np.size(w_lm_g)

#===============================================================================================
planet_data_h = np.loadtxt('Results_Trappist_1h_10.0_100_freq_Imk2_posidonius.txt',comments='#')
# planet_data_h = np.loadtxt('Results_Trappist_1h_13.5_150_freq_Imk2_posidonius.txt',comments='#')
planet_mass_h, planet_radius_h, planet_gyration_radius_squared_h = planet_data_h[0,:] 

w_lm_h   = planet_data_h[1:,0]
ImK2_h   = planet_data_h[1:,1]
ReK2_h   = planet_data_h[1:,2]
size_h   = np.size(w_lm_h)


# #===============================================================================================
# planet_data_venus = np.loadtxt('results_V1_Armann2_atm_freq_Imk2_posidonius.txt',comments='#')
# planet_mass_venus, planet_radius_venus, planet_gyration_radius_squared_venus = planet_data_venus[0,:] 

# w_lm_venus   = planet_data_venus[1:,0]
# ImK2_venus   = planet_data_venus[1:,1]
# ReK2_venus   = planet_data_venus[1:,2]
# size_venus   = np.size(w_lm_venus)


# Constants
Msun = 1.9818e30
Rsun = 6.957e8
Lsun = 3.8275e26
OmSun = 2.6e-6 # Solar rotation rate in rad.s-1
G = 6.67384e-11
jour = 86400. # [s]
an   = 3.15576e7
AU   = 149597870700
Mearth = 5.9736E+24
Rearth = 6378137.
Mjup = 1.8986E+27 # Jupiter mass (kg)
Rjup = 6.991E+07 # Jupiter radius (m)

Mstar   = 0.08 *Msun
Mplanet = 0.772 *Mearth
Rplanet = 0.910 * 6371E+3


#==================================================================================================================

print("Preparing plot...")

# win = Gtk.Window()
# win.connect("destroy", lambda x: Gtk.main_quit())
# win.set_default_size(400,300)
# win.set_title("Embedding in GTK")

# vbox = Gtk.VBox()
# win.add(vbox)

fig = plt.figure(figsize=(13, 5))
ligne = 1
colonne = 2
i = 0 

i=i+1
ax = fig.add_subplot(ligne,colonne,i)
field = 'Imaginary part $k_2$'
ax.plot(w_lm_b ,ImK2_b , label='T-1b', color = 'red')
ax.plot(w_lm_c ,ImK2_c , label='T-1c', color = 'blue')
ax.plot(w_lm_d ,ImK2_d , label='T-1d', color = 'green')
ax.plot(w_lm_e ,ImK2_e , label='T-1e', color = 'yellow')
ax.plot(w_lm_f ,ImK2_f , label='T-1f', color = 'orange')
ax.plot(w_lm_g ,ImK2_g , label='T-1g', color = 'cyan')
ax.plot(w_lm_h ,ImK2_h , label='T-1h', color = 'brown')
# ax.plot(w_lm_SE1 ,ImK2_SE1 , label='SE1', color = 'red')
# ax.plot(w_lm_SE2 ,ImK2_SE2 , label='SE2', color = 'blue')
# ax.plot(w_lm_SE3 ,ImK2_SE3 , label='SE3', color = 'green')
# ax.plot(w_lm_SE4 ,ImK2_SE4 , label='SE4', color = 'orange')
# ax.plot(w_lm_venus ,ImK2_venus , label='Venus', color = 'pink')
ax.set_ylabel(field+" ")
# ax.set_yscale('log')
ax.set_xlabel('Excitative frequency $\sigma$  ($s^{-1}$)')
ax.set_xscale('log')
ax.set_xlim(right=1E-2, left=1E-5)
ax.set_ylim(top=0.006, bottom=0)
# ax.legend(loc=0, prop={'size':8})
#plt.setp(ax.get_xticklabels(), visible=False)


i=i+1
ax = fig.add_subplot(ligne,colonne,i)
field = 'Real part $k_2$'
ax.plot(w_lm_b ,ReK2_b , label='T-1b', color = 'red')
ax.plot(w_lm_c ,ReK2_c , label='T-1c', color = 'blue')
ax.plot(w_lm_d ,ReK2_d , label='T-1d', color = 'green')
ax.plot(w_lm_e ,ReK2_e , label='T-1e', color = 'yellow')
ax.plot(w_lm_f ,ReK2_f , label='T-1f', color = 'orange')
ax.plot(w_lm_g ,ReK2_g , label='T-1g', color = 'cyan')
ax.plot(w_lm_h ,ReK2_h , label='T-1h', color = 'brown')
# ax.plot(w_lm_SE1 ,ReK2_SE1 , label='SE1', color = 'red')
# ax.plot(w_lm_SE2 ,ReK2_SE2 , label='SE2', color = 'blue')
# ax.plot(w_lm_SE3 ,ReK2_SE3 , label='SE3', color = 'green')
# ax.plot(w_lm_SE4 ,ReK2_SE4 , label='SE4', color = 'orange')
# ax.plot(w_lm_venus ,ReK2_venus , label='Venus', color = 'pink')
ax.set_ylabel(field+" ")
ax.set_yscale('log')
ax.set_xlabel('Excitative frequency $\sigma$  ($s^{-1}$)')
# ax.set_xscale('log')
ax.set_xlim(right=1E-2)#, left=1E-15)
# ax.set_ylim(top=2., bottom=0.)
# ax.legend(loc=0, prop={'size':8})
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., prop={'size':12})
#plt.setp(ax.get_xticklabels(), visible=False)

# #==============================================================================

# i=i+1
# ax = fig.add_subplot(ligne,colonne,i)
# field = 'Imaginary part $k_2$'
# # ax.plot(w_lm_b ,ImK2_b , label='b', color = 'red')
# ax.plot(w_lm_c ,ImK2_c , label='c', color = 'b')
# # ax.plot(w_lm_d ,ImK2_d , label='d', color = 'green')
# # ax.plot(w_lm_e ,ImK2_e , label='e', color = 'red')
# # ax.plot(w_lm_f ,ImK2_f , label='f', color = 'red')
# # ax.plot(w_lm_g ,ImK2_g , label='g', color = 'red')
# # ax.plot(w_lm_h ,ImK2_h , label='h', color = 'red')
# ax.set_ylabel(field+" ")
# ax.set_xlabel('Excitative frequency $\sigma$  ($s^{-1}$)')
# ax.set_xscale('log')
# ax.set_xlim(right=1E-2, left=1E-14)
# ax.set_ylim(top=0.5, bottom=0)
# ax.legend(loc=0, prop={'size':8})
# #plt.setp(ax.get_xticklabels(), visible=False)


# i=i+1
# ax = fig.add_subplot(ligne,colonne,i)
# field = 'Real part $k_2$'
# # ax.plot(w_lm_b ,ReK2_b , label='b', color = 'red')
# ax.plot(w_lm_c ,ReK2_c , label='c', color = 'b')
# # ax.plot(w_lm_d ,ReK2_d , label='d', color = 'green')
# # ax.plot(w_lm_e ,ReK2_e , label='e', color = 'red')
# # ax.plot(w_lm_f ,ReK2_f , label='f', color = 'red')
# # ax.plot(w_lm_g ,ReK2_g , label='g', color = 'red')
# # ax.plot(w_lm_h ,ReK2_h , label='h', color = 'red')
# ax.set_ylabel(field+" ")
# ax.set_xlabel('Excitative frequency $\sigma$  ($s^{-1}$)')
# ax.set_xscale('log')
# ax.set_xlim(right=1E-2, left=1E-14)
# ax.set_ylim(top=2., bottom=0.)
# ax.legend(loc=0, prop={'size':8})
# #plt.setp(ax.get_xticklabels(), visible=False)

# #==============================================================================
# i=i+1
# ax = fig.add_subplot(ligne,colonne,i)
# field = 'Imaginary part $k_2$'
# # ax.plot(w_lm_b ,ImK2_b , label='b', color = 'red')
# # ax.plot(w_lm_c ,ImK2_c , label='c', color = 'b')
# ax.plot(w_lm_d ,ImK2_d , label='d', color = 'green')
# # ax.plot(w_lm_e ,ImK2_e , label='e', color = 'red')
# # ax.plot(w_lm_f ,ImK2_f , label='f', color = 'red')
# # ax.plot(w_lm_g ,ImK2_g , label='g', color = 'red')
# # ax.plot(w_lm_h ,ImK2_h , label='h', color = 'red')
# ax.set_ylabel(field+" ")
# ax.set_xlabel('Excitative frequency $\sigma$  ($s^{-1}$)')
# ax.set_xscale('log')
# ax.set_xlim(right=1E-2, left=1E-14)
# ax.set_ylim(top=0.5, bottom=0)
# ax.legend(loc=0, prop={'size':8})
# #plt.setp(ax.get_xticklabels(), visible=False)


# i=i+1
# ax = fig.add_subplot(ligne,colonne,i)
# field = 'Real part $k_2$'
# # ax.plot(w_lm_b ,ReK2_b , label='b', color = 'red')
# # ax.plot(w_lm_c ,ReK2_c , label='c', color = 'b')
# ax.plot(w_lm_d ,ReK2_d , label='d', color = 'green')
# # ax.plot(w_lm_e ,ReK2_e , label='e', color = 'red')
# # ax.plot(w_lm_f ,ReK2_f , label='f', color = 'red')
# # ax.plot(w_lm_g ,ReK2_g , label='g', color = 'red')
# # ax.plot(w_lm_h ,ReK2_h , label='h', color = 'red')
# ax.set_ylabel(field+" ")
# ax.set_xlabel('Excitative frequency $\sigma$  ($s^{-1}$)')
# ax.set_xscale('log')
# ax.set_xlim(right=1E-2, left=1E-14)
# ax.set_ylim(top=2., bottom=0.)
# ax.legend(loc=0, prop={'size':8})
# #plt.setp(ax.get_xticklabels(), visible=False)

# #==============================================================================
# i=i+1
# ax = fig.add_subplot(ligne,colonne,i)
# field = 'Imaginary part $k_2$'
# # ax.plot(w_lm_b ,ImK2_b , label='b', color = 'red')
# # ax.plot(w_lm_c ,ImK2_c , label='c', color = 'b')
# # ax.plot(w_lm_d ,ImK2_d , label='d', color = 'green')
# ax.plot(w_lm_e ,ImK2_e , label='e', color = 'red')
# # ax.plot(w_lm_f ,ImK2_f , label='f', color = 'red')
# # ax.plot(w_lm_g ,ImK2_g , label='g', color = 'red')
# # ax.plot(w_lm_h ,ImK2_h , label='h', color = 'red')
# ax.set_ylabel(field+" ")
# ax.set_xlabel('Excitative frequency $\sigma$  ($s^{-1}$)')
# ax.set_xscale('log')
# ax.set_xlim(right=1E-2, left=1E-14)
# ax.set_ylim(top=0.5, bottom=0)
# ax.legend(loc=0, prop={'size':8})
# #plt.setp(ax.get_xticklabels(), visible=False)


# i=i+1
# ax = fig.add_subplot(ligne,colonne,i)
# field = 'Real part $k_2$'
# # ax.plot(w_lm_b ,ReK2_b , label='b', color = 'red')
# # ax.plot(w_lm_c ,ReK2_c , label='c', color = 'b')
# # ax.plot(w_lm_d ,ReK2_d , label='d', color = 'green')
# ax.plot(w_lm_e ,ReK2_e , label='e', color = 'red')
# # ax.plot(w_lm_f ,ReK2_f , label='f', color = 'red')
# # ax.plot(w_lm_g ,ReK2_g , label='g', color = 'red')
# # ax.plot(w_lm_h ,ReK2_h , label='h', color = 'red')
# ax.set_ylabel(field+" ")
# ax.set_xlabel('Excitative frequency $\sigma$  ($s^{-1}$)')
# ax.set_xscale('log')
# ax.set_xlim(right=1E-2, left=1E-14)
# ax.set_ylim(top=2., bottom=0.)
# ax.legend(loc=0, prop={'size':8})
# #plt.setp(ax.get_xticklabels(), visible=False)

# #==============================================================================
# i=i+1
# ax = fig.add_subplot(ligne,colonne,i)
# field = 'Imaginary part $k_2$'
# # ax.plot(w_lm_b ,ImK2_b , label='b', color = 'red')
# # ax.plot(w_lm_c ,ImK2_c , label='c', color = 'b')
# # ax.plot(w_lm_d ,ImK2_d , label='d', color = 'green')
# # ax.plot(w_lm_e ,ImK2_e , label='e', color = 'red')
# ax.plot(w_lm_f ,ImK2_f , label='f', color = 'red')
# # ax.plot(w_lm_g ,ImK2_g , label='g', color = 'red')
# # ax.plot(w_lm_h ,ImK2_h , label='h', color = 'red')
# ax.set_ylabel(field+" ")
# ax.set_xlabel('Excitative frequency $\sigma$  ($s^{-1}$)')
# ax.set_xscale('log')
# ax.set_xlim(right=1E-2, left=1E-14)
# ax.set_ylim(top=0.5, bottom=0)
# ax.legend(loc=0, prop={'size':8})
# #plt.setp(ax.get_xticklabels(), visible=False)


# i=i+1
# ax = fig.add_subplot(ligne,colonne,i)
# field = 'Real part $k_2$'
# # ax.plot(w_lm_b ,ReK2_b , label='b', color = 'red')
# # ax.plot(w_lm_c ,ReK2_c , label='c', color = 'b')
# # ax.plot(w_lm_d ,ReK2_d , label='d', color = 'green')
# # ax.plot(w_lm_e ,ReK2_e , label='e', color = 'red')
# ax.plot(w_lm_f ,ReK2_f , label='f', color = 'red')
# # ax.plot(w_lm_g ,ReK2_g , label='g', color = 'red')
# # ax.plot(w_lm_h ,ReK2_h , label='h', color = 'red')
# ax.set_ylabel(field+" ")
# ax.set_xlabel('Excitative frequency $\sigma$  ($s^{-1}$)')
# ax.set_xscale('log')
# ax.set_xlim(right=1E-2, left=1E-14)
# ax.set_ylim(top=2., bottom=0.)
# ax.legend(loc=0, prop={'size':8})
# #plt.setp(ax.get_xticklabels(), visible=False)

# #==============================================================================
# i=i+1
# ax = fig.add_subplot(ligne,colonne,i)
# field = 'Imaginary part $k_2$'
# # ax.plot(w_lm_b ,ImK2_b , label='b', color = 'red')
# # ax.plot(w_lm_c ,ImK2_c , label='c', color = 'b')
# # ax.plot(w_lm_d ,ImK2_d , label='d', color = 'green')
# # ax.plot(w_lm_e ,ImK2_e , label='e', color = 'red')
# # ax.plot(w_lm_f ,ImK2_f , label='f', color = 'red')
# ax.plot(w_lm_g ,ImK2_g , label='g', color = 'red')
# # ax.plot(w_lm_h ,ImK2_h , label='h', color = 'red')
# ax.set_ylabel(field+" ")
# ax.set_xlabel('Excitative frequency $\sigma$  ($s^{-1}$)')
# ax.set_xscale('log')
# ax.set_xlim(right=1E-2, left=1E-14)
# ax.set_ylim(top=0.5, bottom=0)
# ax.legend(loc=0, prop={'size':8})
# #plt.setp(ax.get_xticklabels(), visible=False)


# i=i+1
# ax = fig.add_subplot(ligne,colonne,i)
# field = 'Real part $k_2$'
# # ax.plot(w_lm_b ,ReK2_b , label='b', color = 'red')
# # ax.plot(w_lm_c ,ReK2_c , label='c', color = 'b')
# # ax.plot(w_lm_d ,ReK2_d , label='d', color = 'green')
# # ax.plot(w_lm_e ,ReK2_e , label='e', color = 'red')
# # ax.plot(w_lm_f ,ReK2_f , label='f', color = 'red')
# ax.plot(w_lm_g ,ReK2_g , label='g', color = 'red')
# # ax.plot(w_lm_h ,ReK2_h , label='h', color = 'red')
# ax.set_ylabel(field+" ")
# ax.set_xlabel('Excitative frequency $\sigma$  ($s^{-1}$)')
# ax.set_xscale('log')
# ax.set_xlim(right=1E-2, left=1E-14)
# ax.set_ylim(top=2., bottom=0.)
# ax.legend(loc=0, prop={'size':8})
# #plt.setp(ax.get_xticklabels(), visible=False)

# #==============================================================================

# i=i+1
# ax = fig.add_subplot(ligne,colonne,i)
# field = 'Imaginary part $k_2$'
# # ax.plot(w_lm_b ,ImK2_b , label='b', color = 'red')
# # ax.plot(w_lm_c ,ImK2_c , label='c', color = 'b')
# # ax.plot(w_lm_d ,ImK2_d , label='d', color = 'green')
# # ax.plot(w_lm_e ,ImK2_e , label='e', color = 'red')
# # ax.plot(w_lm_f ,ImK2_f , label='f', color = 'red')
# # ax.plot(w_lm_g ,ImK2_g , label='g', color = 'red')
# ax.plot(w_lm_h ,ImK2_h , label='h', color = 'red')
# ax.set_ylabel(field+" ")
# ax.set_xlabel('Excitative frequency $\sigma$  ($s^{-1}$)')
# ax.set_xscale('log')
# ax.set_xlim(right=1E-2, left=1E-14)
# ax.set_ylim(top=0.5, bottom=0)
# ax.legend(loc=0, prop={'size':8})
# #plt.setp(ax.get_xticklabels(), visible=False)


# i=i+1
# ax = fig.add_subplot(ligne,colonne,i)
# field = 'Real part $k_2$'
# # ax.plot(w_lm_b ,ReK2_b , label='b', color = 'red')
# # ax.plot(w_lm_c ,ReK2_c , label='c', color = 'b')
# # ax.plot(w_lm_d ,ReK2_d , label='d', color = 'green')
# # ax.plot(w_lm_e ,ReK2_e , label='e', color = 'red')
# # ax.plot(w_lm_f ,ReK2_f , label='f', color = 'red')
# # ax.plot(w_lm_g ,ReK2_g , label='g', color = 'red')
# ax.plot(w_lm_h ,ReK2_h , label='h', color = 'red')
# ax.set_ylabel(field+" ")
# ax.set_xlabel('Excitative frequency $\sigma$  ($s^{-1}$)')
# ax.set_xscale('log')
# ax.set_xlim(right=1E-2, left=1E-14)
# ax.set_ylim(top=2., bottom=0.)
# ax.legend(loc=0, prop={'size':8})
# #plt.setp(ax.get_xticklabels(), visible=False)

# #==============================================================================

# plt.show()

# win.show_all()
# Gtk.main()
plt.savefig('plot_.png')