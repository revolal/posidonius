
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

PI = np.pi
G = 1.
G_SI = 6.67428e-11  #.   m^3.kg^-1.s^-2 (S.I. units)
AU = 1.49597870700e11 #.   m

DAY = 24.*3600 #.   s

M_SUN =  1.9818e30 #.   Kg
M_EARTH = 5.972e24 #.   Kg
M_JUP = 1.898e27

R_SUN = 696340000 # meter


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

    fm1_27 = [-1.19049, -2.02522, -2.84043 ,-3.64962 ,4.45614, -5.26125 ,-6.06552, -6.86925, -7.67261 , -8.47571 , -9.27861 , -10.0814, -10.884]
    fm1_31 = [0.42839, 2.48401, 3.28326, 4.08371, 4.88471, 5.68601, 6.48749, 7.28909, 8.09077, 8.89251, 9.69429, 10.4961, 11.2979]

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

