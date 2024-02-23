#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 13:57:31 2024

This mean-motion resonance (MMR) analysis script was created based on Deck (2013),
note that this script is valid only for first-order resonances (m=1).

@author: leonkwok
"""

import sys
import time
# print(sys.path)
# sys.path.append("/Users/leonkwok/opt/anaconda3/lib/python3.9/site-packages")
# This is a temporary solution for adding the path which the packages were installed in Py 3.11

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.ascii as ascii
import os
import math
from matplotlib.ticker import MultipleLocator
from matplotlib.animation import FuncAnimation
import pandas as pd

# %%

PI = np.pi
G = 1.
G_SI = 6.67428e-11  #.   m^3.kg^-1.s^-2 (S.I. units)
AU = 1.49597870700e11 #.   m

DAY = 24.*3600 #.   s

M_SUN =  1.9818e30 #.   Kg
M_EARTH = 5.972e24 #.   Kg
M_JUP = 1.898e27

R_SUN = 696340000 # meter

# G_alt = G_SI/(AU*AU*AU) * M_SUN * (DAY*DAY) #.   Gravitational constant with Mercury units

fm1_27 = [-1.19049, -2.02522, -2.84043 ,-3.64962 ,4.45614, -5.26125 ,-6.06552, -6.86925, -7.67261 , -8.47571 , -9.27861 , -10.0814, -10.884]
fm1_31 = [0.42839, 2.48401, 3.28326, 4.08371, 4.88471, 5.68601, 6.48749, 7.28909, 8.09077, 8.89251, 9.69429, 10.4961, 11.2979]

# %%
'''
M: Stelalr mass
m: Planetary mass
P: Orbital period
e: Eccentricity
l: mean longitude
w: longitude of periastron
'''
# %%

def MMR_1st_order(M, m1, m2, a1, a2, e1, e2, l1, l2, w1, w2, m=0):

    period_1 = 2*PI*np.sqrt(a1**3/(G*M))
    period_2 = 2*PI*np.sqrt(a2**3/(G*M))

    # A switch for which kind of first order resonance
    # simply set m = 0 and let the program decide itself    
    if m != 0:
        if round(period_1/(period_2-period_1)) != m:
            print(f"The periods do not corespond to the given resonance number: m(P)={round(period_1/(period_2-period_1))}, m={m}")
            return 1
    else:
        m = round(period_1/(period_2-period_1))
        
    if m == 0:
        return (np.nan,) * 3

    eps1 = m1/M
    eps2 = m2/M
    epsp = eps1 + eps2

    p1 = -w1
    p2 = -w2

    alpha = a1/a2

    # Eq.2
    Lambda1 = m1 * np.sqrt(G*M*a1)
    Lambda2 = m2 * np.sqrt(G*M*a2)
    
    P1 = m1 *np.sqrt(G*M*a1)*(1-np.sqrt(1-e1**2))
    P2 = m2 *np.sqrt(G*M*a2)*(1-np.sqrt(1-e2**2))

    # Eq.9
    theta = (m+1)*l2 - m*l1
    theta1 = l1 

    Theta = Lambda2/(m+1)
    Theta1 = m/(m+1)*Lambda2 + Lambda1

    x1 = np.sqrt(2*P1)*np.cos(p1)
    x2 = np.sqrt(2*P2)*np.cos(p2)
    y1 = np.sqrt(2*P1)*np.sin(p1)
    y2 = np.sqrt(2*P2)*np.sin(p2)
    
    x1_cap = x1/np.sqrt(Theta1)
    x2_cap = x2/np.sqrt(Theta1)
    y1_cap = y1/np.sqrt(Theta1)
    y2_cap = y2/np.sqrt(Theta1)


    alpha_res = (m/(m+1))**(2/3)
    xi = m1/m2
    
    Theta_cap = Theta/Theta1
    Theta_bar = alpha_res /(m*(alpha_res + xi))

    s = np.sqrt(alpha/alpha_res)
    delta_Theta = xi/m *(1 - s) * 1/(alpha_res + xi + s*xi*(1+xi/alpha_res))

    delta1 = 1/((m+1) * Theta_cap)**2 * fm1_27[m-1]/(np.sqrt(1 - m*Theta_cap))
    delta2 = 1/((m+1) * Theta_cap)**2 * fm1_31[m-1]/(np.sqrt((m+1) * Theta_cap))

    delta1_bar = 1/((m+1) * Theta_bar)**2 * fm1_27[m-1]/(np.sqrt(1 - m*Theta_bar))
    delta2_bar = 1/((m+1) * Theta_bar)**2 * fm1_31[m-1]/(np.sqrt((1 + m)*Theta_bar))
    delta_bar = np.sqrt(delta1_bar**2 + delta2_bar**2)

    r1 = (delta1_bar*x1_cap + delta2_bar*x2_cap)/delta_bar
    s1 = (delta1_bar*y1_cap + delta2_bar*y2_cap)/delta_bar
    r2 = (delta2_bar*x1_cap - delta1_bar*x2_cap)/delta_bar
    s2 = (delta2_bar*y1_cap - delta1_bar*y2_cap)/delta_bar

    # Psi2 is not needed, noted that r2 and s2 are conserved
    Psi = 1/2 * (r1**2 + s1**2)
    Psi2 = 1/2 * (r2**2 + s2**2)
    psi = np.arctan2(s1, r1)
    
    phi = theta + psi
    gamma = -theta

    Phi = Psi
    Gamma = Psi - delta_Theta

    R = abs(fm1_27[m-1]/fm1_31[m-1])
    Q = epsp**(2./3.) * xi * (fm1_31[m-1]**2. / (9.*(m+1.)*(1.+xi)**2.))**(1./3.) * (alpha_res**(5./6.) / m) * ((R**2 + xi*np.sqrt(alpha_res))/(alpha_res+xi)**5)**(1/3)

    # "rescale the variables", that's what she said
    Phi_p = Phi/Q 
    Gamma_p = Gamma/Q
    X = np.sqrt(2*Phi_p) * np.cos(phi)
    Y = np.sqrt(2*Phi_p) * np.sin(phi)

    return X, Y, Gamma_p

def Period_ratio(M, m1, m2, a1, a2):
    period_1 = 2*PI*np.sqrt(a1**3/(G*M))
    period_2 = 2*PI*np.sqrt(a2**3/(G*M))

    return period_2/period_1

def find_X_fix(Gamma):
    # a specific function for solving X (Eq. 38 from Deck 2013)
    roots = np.roots([1., 0., -2*Gamma, 2.])
    roots = roots[np.isreal(roots)].real
    if len(roots) == 1:
        return np.array([roots[0], np.nan, np.nan])
    else:
        return roots
    
def find_X_sol (Y, Gamma_p, H_p):
    roots = np.roots([1/4., 0., 0.5*(Y**2) - Gamma_p, 2, Y**4. + Gamma_p**2 + 2*H_p])
    roots = roots[np.isreal(roots)].real
    tem = np.full(4, np.NAN)
    tem[:len(roots)] = roots
    return tem

def Hamiltonian_prime (X, Y, Gamma_p):
    H_p = -0.5*(0.5*(X*X + Y*Y) - Gamma_p)**2 - X
    return H_p

def func_Hamiltonian_solving_X (X, *data):
    Y, Gamma_p, H_sep = data
    return Hamiltonian_prime(X, Y, Gamma_p) - H_sep

def print_progress(progress, i, length):
    # Construct the progress message
    message = f"Progress: {progress:.2f}%, Current frame: {i}/{length}"
    # Clear the current line
    sys.stdout.write('\r' + ' ' * (len(message) + 10))  # Overwrite with spaces
    sys.stdout.write('\r' + message)  # Carriage return and print
    sys.stdout.flush()


# %%: Data retreive, everything in manual

filename = "./data/r_21_a02_samem_1J1J_history.txt"
base_dir = os.path.dirname(filename)
filename_without_extension = os.path.splitext(os.path.basename(filename))[0]

raw_data = ascii.read(filename)
number = np.array(raw_data['planet'])
raw_t = np.array(raw_data['current_time'][np.where(number==1)])
raw_a = np.array(raw_data['semi-major_axis_AU'])
raw_e = np.array(raw_data['eccentricity'])
raw_w = np.array(raw_data['longitude_of_perihelion'])
raw_l = np.array(raw_data['mean_longitude'])
raw_Omega = 2*PI/np.array(raw_data['stellar_rotation_period'][np.where(number==1)]) # in DAYS^-1
raw_n = 2*PI/np.array(raw_data['orbital_period']) # in DAYS^-1
raw_delta_L = np.array(raw_data['conservation_of_angular_momentum'][np.where(number==1)])
# raw_M = np.array(raw_data['mean_anomaly'])

M = np.array(raw_data['stellar_mass'][0])
m1 = np.array(raw_data['planet_mass'][np.where(number==1)][0])
m2 = np.array(raw_data['planet_mass'][np.where(number==2)][0])

Omega_SI = raw_Omega/86400.
corot_r = pow((G_SI*(M*M_SUN)/(Omega_SI*Omega_SI)), 1./3.)/AU

a_1 = raw_a[np.where(number==1)]
e_1 = raw_e[np.where(number==1)]
l_1 = raw_l[np.where(number==1)]
w_1 = raw_w[np.where(number==1)]
r_Omega_n_1 = raw_Omega/raw_n[np.where(number==1)]

a_2 = raw_a[np.where(number==2)]
e_2 = raw_e[np.where(number==2)]
l_2 = raw_l[np.where(number==2)]
w_2 = raw_w[np.where(number==2)]
r_Omega_n_2 = raw_Omega/raw_n[np.where(number==2)]

# %%: Solve for Resonance

length = len(raw_t)

X = np.zeros(length)
Y = np.zeros(length)
Gamma_p = np.zeros(length)
P_ratio = np.zeros(length)
X_fix = np.full((length, 3), np.nan)
X_sol = np.zeros((length, 4))
X_sol_tem = np.zeros((length, 4))
X = np.zeros(length)
Res_state = np.full(length, np.nan)

separatrix_points = []

# We manually input M=1 for stellar mass
separatrix = False
for i in range(length):
    X[i], Y[i], Gamma_p[i] = MMR_1st_order(M, m1, m2, a_1[i], a_2[i], e_1[i], e_2[i], l_1[i], l_2[i], w_1[i], w_2[i])
    
    if np.isnan(Gamma_p[i]):
        # Resonance = False
        P_ratio[i] = Period_ratio(M, m1, m2, a_1[i], a_2[i])
    else:
        X_fix[i] = find_X_fix(Gamma_p[i])
        X_fix[i] = np.sort(X_fix[i])
        H_p = Hamiltonian_prime (X[i], Y[i], Gamma_p[i])
        X_sol[i] = find_X_sol(0., Gamma_p[i], H_p)
        X_sol[i] = np.sort(X_sol[i])
        P_ratio[i] = Period_ratio(M, m1, m2, a_1[i], a_2[i])
        
X3 = X_fix[:,2]
X_cross_1 = -X3 + 2./np.sqrt(X3)
X_cross_2 = -X3 - 2./np.sqrt(X3)

# %% Script for generating a video of resonance plots

def resonance_series(frame_index):
    # Calculate the actual index (i) to use, scaling by 10 to get every 10th iteration
    global iteration_counter  # Reference the global counter
    i = valid_indices[frame_index]
    if (resonance_plot_real == True and Gamma_p[i] >= 3./2.): 
        ax.clear()
        Gamma_p_tem = Gamma_p[i]  # Example value, adjust as necessary
        X3_tem = X_fix[i, 2]   
        X_sep_Y0 = -X3_tem - 2./np.sqrt(X3_tem)
        X_cross_1_tem = -X3_tem + 2./np.sqrt(X3_tem)
        X_cross_2_tem = -X3_tem - 2./np.sqrt(X3_tem)
        H_sep_value = Hamiltonian_prime(X_sep_Y0, 0., Gamma_p_tem)
        
        x = np.linspace(-math.ceil(abs(X_cross_2_tem))*1.05, abs(X_cross_2_tem), 400)  # Adjust range and density as necessary
        y = np.linspace(-math.ceil(abs(X_cross_2_tem))*1.05, abs(X_cross_2_tem), 400)  # Adjust range and density as necessary
        X_plot, Y_plot = np.meshgrid(x, y)
        
        H = Hamiltonian_prime(X_plot, Y_plot, Gamma_p_tem)
    
        contour_outside = plt.contour(X_plot, Y_plot, H, levels=np.linspace(H_sep_value - np.abs(H_sep_value)*0.5, H_sep_value, 2), colors='black', linestyles='solid')
        contour_inside = plt.contour(X_plot, Y_plot, H, levels=np.linspace(H_sep_value, np.max(H), 5), colors='black', linestyles='solid')
        contour = ax.contour(X_plot, Y_plot, H, levels=[H_sep_value], colors='red', linewidths=2, linestyles='solid')
    
        separatrix_paths = contour.collections[0].get_paths()
        if separatrix_paths:
            separatrix_path = separatrix_paths[0]
            Res_state[i] = separatrix_path.contains_point((X[i], Y[i]))
            
        for x_sol in X_sol[i]:  # Iterate through solutions for this gamma_prime
            if separatrix_path.contains_point((X[i], Y[i])) == separatrix_path.contains_point((x_sol, 0.)):
                X_Y0_final[i] = x_sol
                break
        
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))  # e.g., minor ticks at 0.2 intervals
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.tick_params(axis='both', which='major', length=7, width=2)  # Major ticks
        ax.tick_params(axis='both', which='minor', length=4, width=1)   
        
        zorder_top = 10
        
        ax.scatter(X_fix[i, 0], 0., c="black", s=30, marker='o', zorder=zorder_top)
        ax.text(X_fix[i, 0], 0., 'X1', fontsize=20, verticalalignment='bottom', horizontalalignment='left')
        
        ax.scatter(X_fix[i, 1], 0., c="black", s=30, marker='o', zorder=zorder_top)
        ax.text(X_fix[i, 1], 0., 'X2', fontsize=20, verticalalignment='bottom', horizontalalignment='left')
        
        ax.scatter(X_fix[i, 2], 0., c="black", s=30, marker='o', zorder=zorder_top)
        ax.text(X_fix[i, 2], 0., 'X3', fontsize=20, verticalalignment='bottom', horizontalalignment='left')
        
        x_coords = [X_cross_1_tem, X_cross_2_tem]
        y_coords = [0, 0]
        
        ax.plot(x_coords, y_coords, linestyle='-.', color='blue', linewidth=1.5, zorder=5)
              
        time_tem = raw_t[i]
        title_text = r"$\Gamma^\prime = {:.2f}$, Time: {:.2e}".format(Gamma_p_tem, time_tem)
        ax.set_title(title_text, fontsize=22)
        
        ax.scatter(X[i], Y[i], c="purple", s=200, marker='*', zorder=10)
            
        ax.set_xlabel('X', fontsize=20)
        ax.set_ylabel('Y', fontsize=20)
        ax.tick_params(labelsize=20, direction='in', top=True, right=True, which='both')
        
        X_min = -math.ceil(abs(X_cross_2_tem))
        X_max = math.ceil(abs(X_cross_2_tem))
        ax.set_xlim(X_min, X_max)
        ax.set_ylim(X_min, X_max)
               
        total_valid_frames = len(valid_indices)
        progress = iteration_counter / total_valid_frames * 100  # Calculate progress percentage
        print_progress(progress, i, length)  
        
        iteration_counter += 1
        
    elif (resonance_plot_real == True and Gamma_p[i] <= 3./2.): 
        Res_state[iteration_counter] = False
    else:
        print("This script is switched off.")
        return

resonance_plot_real = True

if resonance_plot_real:    
    new_base_name = filename_without_extension.replace('history', 'resonance')    
    res_gif = os.path.join(base_dir, new_base_name + ".gif")
    res_txt = os.path.join(base_dir, new_base_name + ".txt")
    
    start_time = time.time()
    iteration_counter = 0
    
    valid_indices = [i for i, gamma in enumerate(Gamma_p) if gamma >= 3./2.]
    # it means we take a snapshot each 10 valid frames
    valid_indices = [i for i in valid_indices if (valid_indices.index(i) % 10 == 0)] 
    num_frames = len(valid_indices)
    
    fig, ax = plt.subplots(figsize=(10,8), dpi = 120)
    # fig.tight_layout()
    X_Y0_final = np.full(length, np.nan)
    anim = FuncAnimation(fig, resonance_series, frames=range(num_frames), repeat=False)
    
    print("Saving animation, please wait...")
    anim.save(res_gif, writer='Pillow', fps=24)
    end_time = time.time()  # Record the end time
    duration = (end_time - start_time)/60  # Calculate the duration
    print(f"\nAnimation saved. Total generation time: {duration:.2f} minutes.")
    
    data = pd.DataFrame()
    data['Gamma_p'] = Gamma_p
    data['X'] = X
    data['Y'] = Y
    data['X1_Y0'] = X_sol[:, 0]
    data['X2_Y0'] = X_sol[:, 1]
    data['X3_Y0'] = X_sol[:, 2]
    data['X4_Y0'] = X_sol[:, 3]
    data['X_Y0_final'] = X_Y0_final
    data['X1_fix'] = X_fix[:, 0]
    data['X2_fix'] = X_fix[:, 1]
    data['X3_fix'] = X_fix[:, 2]
    data['X_cross_1'] = X_cross_1
    data['X_cross_2'] = X_cross_2
    data['res_state'] = Res_state
    
    data.to_csv(res_txt, sep="\t", index=False)

# %% Script for generating a single resonance plot, especially for checking resonance for some specific time frame.

resonance_plot_real = False

i = 25647
if (resonance_plot_real == True and Gamma_p[i] >= 3./2.): 
    fig, ax = plt.subplots(figsize=(10,8), dpi = 120, sharex=True)
    
    Gamma_p_tem = Gamma_p[i]  # Example value, adjust as necessary
    X3_tem = X_fix[i, 2]   
    X_sep_Y0 = -X3_tem - 2./np.sqrt(X3_tem)
    X_cross_1_tem = -X3_tem + 2./np.sqrt(X3_tem)
    X_cross_2_tem = -X3_tem - 2./np.sqrt(X3_tem)
    H_sep_value = Hamiltonian_prime(X_sep_Y0, 0., Gamma_p_tem)
    
    x = np.linspace(-math.ceil(abs(X_cross_2_tem))*1.05, abs(X_cross_2_tem), 400)  # Adjust range and density as necessary
    y = np.linspace(-math.ceil(abs(X_cross_2_tem))*1.05, abs(X_cross_2_tem), 400)  # Adjust range and density as necessary
    X_plot, Y_plot = np.meshgrid(x, y)
    
    H = Hamiltonian_prime(X_plot, Y_plot, Gamma_p_tem)

    contour_outside = ax.contour(X_plot, Y_plot, H, levels=np.linspace(H_sep_value - np.abs(H_sep_value)*0.5, H_sep_value, 2), colors='black', linestyles='solid')
    contour_inside = ax.contour(X_plot, Y_plot, H, levels=np.linspace(H_sep_value, np.max(H), 5), colors='black', linestyles='solid')
        
    contour = ax.contour(X_plot, Y_plot, H, levels=[H_sep_value], colors='red', linewidths=2, linestyles='solid')

    separatrix_paths = contour.collections[0].get_paths()
    if separatrix_paths:
        separatrix_path = separatrix_paths[0]
        
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))  # e.g., minor ticks at 0.2 intervals
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.tick_params(axis='both', which='major', length=7, width=2)  # Major ticks
    ax.tick_params(axis='both', which='minor', length=4, width=1)   
    
    ax.scatter(X_sol[i,0], 0., color = "orange", s = 150, marker = "x")
    ax.scatter(X_sol[i,1], 0., color = "orange", s = 150, marker = "x")
    ax.scatter(X_sol[i,2], 0., color = "orange", s = 150, marker = "x")
    ax.scatter(X_sol[i,3], 0., color = "orange", s = 150, marker = "x")
    
    point = (X_sol[i,0], 0.)  # The point you want to test
    inside = separatrix_path.contains_point(point)
    
    zorder_top = 10
    
    ax.scatter(X_fix[i, 0], 0., c="black", s=30, marker='o', zorder=zorder_top)
    ax.text(X_fix[i, 0], 0., 'X1', fontsize=20, verticalalignment='bottom', horizontalalignment='left')
    
    ax.scatter(X_fix[i, 1], 0., c="black", s=30, marker='o', zorder=zorder_top)
    ax.text(X_fix[i, 1], 0., 'X2', fontsize=20, verticalalignment='bottom', horizontalalignment='left')
    
    ax.scatter(X_fix[i, 2], 0., c="black", s=30, marker='o', zorder=zorder_top)
    ax.text(X_fix[i, 2], 0., 'X3', fontsize=20, verticalalignment='bottom', horizontalalignment='left')
    
    x_coords = [X_cross_1_tem, X_cross_2_tem]
    y_coords = [0, 0]
    
    ax.plot(x_coords, y_coords, linestyle='-.', color='blue', linewidth=1.5, zorder=5)
        
    time_tem = raw_t[i]
    
    title_text = r"$\Gamma^\prime = {:.2f}$, Time: {:.2e}".format(Gamma_p_tem, time_tem)
    ax.set_title(title_text, fontsize=22)
    
    ax.scatter(X[i], Y[i], c="purple", s=200, marker='*', zorder=10)
    print(i, X[i], Y[i])
    print(separatrix_path.contains_point((X[i], Y[i])))
        
    ax.set_xlabel('X', fontsize=20)
    ax.set_ylabel('Y', fontsize=20)
    ax.tick_params(labelsize=20, direction='in', top=True, right=True, which='both')
    
    X_min = -math.ceil(abs(X_cross_2_tem))
    X_max = math.ceil(abs(X_cross_2_tem))
    ax.set_xlim(X_min, X_max)
    ax.set_ylim(X_min, X_max)
    
    # ax.set_xlim(-math.ceil(abs(X_cross_2_tem))*1.1, abs(X_cross_2_tem))
    # ax.set_ylim(-math.ceil(abs(X_cross_2_tem))*1.1, math.ceil(abs(X_cross_2_tem)))

    plt.show()
    
elif resonance_plot_real == False:
    print("This script is switched off.")

else:
    print("Not in resonance.")
    # plt.ioff()

# %%
 
fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize=(14, 14), sharex=True)

ax1.plot(raw_t, a_1, color='red', label='1')
ax1.plot(raw_t, a_1*(1 + e_1), color='red')
ax1.plot(raw_t, a_1*(1 - e_1), color='red')

ax1.plot(raw_t, a_2, color='mediumblue', label='2')
ax1.plot(raw_t, a_2*(1 + e_2), color='mediumblue')
ax1.plot(raw_t, a_2*(1 - e_2), color='mediumblue')

ax1.plot(raw_t, corot_r, color='grey', linestyle = '--', label='CR')
ax1.set_ylabel("a (AU)", fontsize=25)
ax1.legend(fontsize=20)
ax1.tick_params(labelsize=25, direction='in', top=True, right=True, which='both')
ax1.set_xlim(1e2, )
ax1.set_xscale('log')
plt.setp(ax1.get_xticklabels(), visible=False)  # Hide x-axis labels for ax1

ax2.plot(raw_t, e_1, color='red')
ax2.plot(raw_t, e_2, color='mediumblue')
ax2.set_ylabel("e", fontsize=25)
ax2.tick_params(labelsize=25, direction='in', top=True, right=True, which='both')
ax2.set_xlim(1e2, )
ax2.set_ylim(1e-4, )
ax2.set_xscale('log')
ax2.set_yscale('log')
plt.setp(ax2.get_xticklabels(), visible=False)  # Hide x-axis labels for ax2

ax3.plot(raw_t, r_Omega_n_1, color='red')
ax3.plot(raw_t, r_Omega_n_2, color='mediumblue')
ax3.axhline(y=18/11, color='grey', linestyle='--', label=r'$\dot e = 0$')
ax3.legend(fontsize=20)
ax3.set_ylabel(r"$\Omega/n$", fontsize=25)
ax3.set_xlabel("Time (yrs)", fontsize=25)  # Set x-axis label for ax3 only
ax3.tick_params(labelsize=25, direction='in', top=True, right=True, which='both')
ax3.set_xlim(1e2, )
ax3.set_xscale('log')

plt.tight_layout()
plt.subplots_adjust(hspace = 0)  
plt.savefig(filename_without_extension + "_evolution.jpg")
plt.show()
# plt.ioff()

# %%
 
fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize=(14, 14), sharex=True)

starting_time = int(1e5/100) # data was saved by each 100 years (adjustable)

ax1.plot(raw_t[starting_time:], a_1[starting_time:], color='red', label='1')
ax1.plot(raw_t[starting_time:], (a_1*(1 + e_1))[starting_time:], color='red')
ax1.plot(raw_t[starting_time:], (a_1*(1 - e_1))[starting_time:], color='red')

ax1.plot(raw_t[starting_time:], a_2[starting_time:], color='mediumblue', label='2')
ax1.plot(raw_t[starting_time:], (a_2*(1 + e_2))[starting_time:], color='mediumblue')
ax1.plot(raw_t[starting_time:], (a_2*(1 - e_2))[starting_time:], color='mediumblue')

ax1.plot(raw_t[starting_time:], corot_r[starting_time:], color='grey', linestyle = '--', label='CR')
ax1.set_ylabel("a (AU)", fontsize=30)
ax1.legend(fontsize=20)
ax1.tick_params(labelsize=25, direction='in', top=True, right=True, which='both')
ax1.set_xlim(starting_time*100, )
ax1.set_xscale('log')
plt.setp(ax1.get_xticklabels(), visible=False)  # Hide x-axis labels for ax1

ax2.plot(raw_t[starting_time:], e_1[starting_time:], color='red')
ax2.plot(raw_t[starting_time:], e_2[starting_time:], color='mediumblue')
ax2.set_ylabel("e", fontsize=30)
ax2.tick_params(labelsize=25, direction='in', top=True, right=True, which='both')
ax2.set_xlim(starting_time*100, )
ax2.set_ylim(1e-4, )
ax2.set_xscale('log')
ax2.set_yscale('log')
plt.setp(ax2.get_xticklabels(), visible=False)  # Hide x-axis labels for ax2

ax3.plot(raw_t[starting_time:], r_Omega_n_1[starting_time:], color='red')
ax3.plot(raw_t[starting_time:], r_Omega_n_2[starting_time:], color='mediumblue')
ax3.axhline(y=18/11, color='grey', linestyle='--', label=r'$\dot e = 0$')
ax3.legend(fontsize=20)
ax3.set_ylabel(r"$\Omega/n$", fontsize=30)
ax3.set_xlabel("Time (yrs)", fontsize=30)  # Set x-axis label for ax3 only
ax3.tick_params(labelsize=25, direction='in', top=True, right=True, which='both')
ax3.set_xlim(starting_time*100, )
ax3.set_xscale('log')

plt.tight_layout()
plt.subplots_adjust(hspace = 0)   
plt.savefig(filename_without_extension + "_evolution_enlarged.jpg")
plt.show()
# plt.ioff()

# %%

for i in range(length):
    if Gamma_p[i] >= 3./2:
        pass
    else:
        X_Y0_final[i] = np.min(X_sol[i][~np.isnan(X_sol[i])]) 

fig, ax = plt.subplots(figsize=(10,8), dpi = 80, sharex=True)
colors = np.linspace(raw_t[0], raw_t[-1], length)

sc = ax.scatter(Gamma_p, X_Y0_final, c=colors, cmap='inferno', s = 4, marker = "o")
# sc = ax.scatter(Gamma_p, X_sol[:,0], c=colors, cmap='inferno', s = 4, marker = "o")
# ax.scatter(Gamma_p, X_sol[:,1], c=colors, cmap='inferno', s = 4, marker = "v")
# ax.scatter(Gamma_p, X_sol[:,2], c=colors, cmap='inferno', s = 4, marker = "^")
# ax.scatter(Gamma_p, X_sol[:,3], c=colors, cmap='inferno', s = 4, marker = "*")

cbar = fig.colorbar(sc, orientation='horizontal', pad=0.1)
cbar.set_label('Time', fontsize=25)
cbar.ax.tick_params(labelsize=20)
sc.set_clim(raw_t[0], raw_t[-1])

ax.scatter(Gamma_p, X_fix[:,0], c="dimgrey", s = 1, marker = ",", label = 'X1')
if not np.all(np.isnan(X_fix[:,1])):
    ax.scatter(Gamma_p, X_fix[:,1], c="darkgrey", s = 1, marker = ",", label = 'X2')
    ax.scatter(Gamma_p, X_fix[:,2], c="lightgrey", s = 1, marker = ",", label = 'X3')
    ax.scatter(Gamma_p, X_cross_1, c="red", s = 1, label = 'Resonance boundary')
    ax.scatter(Gamma_p, X_cross_2, c="red", s = 1)

ax.set_xlabel(r"$\Gamma^{\prime}$", fontsize=25)
ax.set_ylabel("X", fontsize=25)
ax.legend(fontsize=14)
ax.tick_params(labelsize=20, direction='in', top=True, right=True, which='both')

plt.tight_layout()
plt.savefig(filename_without_extension + "_resonance.jpg")
plt.show()
# plt.ioff()

# %%

fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize=(14, 14), sharex=True)

ax1.scatter(raw_t, Gamma_p, color='green', s=1)
ax1.set_ylabel(r"$\Gamma^{\prime}$", fontsize=30)
ax1.tick_params(labelsize=25, direction='in', top=True, right=True, which='both')
ax1.set_xlim(1e2, )
ax1.set_xscale('log')
plt.setp(ax1.get_xticklabels(), visible=False)  

ax2.scatter(raw_t, P_ratio, color='purple', s=1)
ax2.set_xlabel("Time (yrs)", fontsize=30)
ax2.set_ylabel("Orbital Period Ratio", fontsize=30)
ax2.tick_params(labelsize=25, direction='in', top=True, right=True, which='both')
ax2.set_xlim(1e2, )
ax2.set_ylim(1, 5)
ax2.set_xscale('log')

ax3.plot(raw_t, raw_delta_L, color='orange')
ax3.set_ylabel(r"$\Delta L$", fontsize=30)
ax3.set_xlabel("Time (yrs)", fontsize=30)  # Set x-axis label for ax3 only
ax3.tick_params(labelsize=25, direction='in', top=True, right=True, which='both')
ax3.set_xlim(1e2, )
ax3.set_xscale('log')
ax3.set_yscale('log')

plt.tight_layout()
plt.subplots_adjust(hspace = 0)   
plt.savefig(filename_without_extension + "_analysis.jpg")
plt.show()
# plt.ioff()

























