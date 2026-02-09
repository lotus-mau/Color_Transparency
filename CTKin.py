# -*- coding: utf-8 -*-
"""
Kinematic Phase Space and Histogram script

@author: Lotus
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import CTHelp as cth

plt.close('all')

# CONSTANTS

m_p  = 0.9382720813     # proton mass
m_n = 0.93956563        # neutron mass
m_pi = 0.139570611      # charged pion mass

# INPUTS

A = 12; Z = 6           # nucleon and atomic number : Carbon-12
m_A = Z*m_p + (A-Z)*m_n # target mass in GeV/c^2
m_A_1 = (A-1)*m_p       #

# TEST INPUT

# E_beam = np.arange(4, 6, 0.1)
# E_beam = np.array([6.0])
# Q2 = np.arange(1, 5, 0.02)
# E_prime = np.arange(0.5, E_beam[-1], 0.01)

# VARIABLE INPUT

# E_beam = np.arange(10.5, 11.1, 0.1)
E_beam = np.array([11.0])
Q2 = np.arange(4.0, 5.0, 0.01)
E_prime = np.arange(0.3, E_beam[-1], 0.01)

# DEFINITIONS

def Theta_e(q2, E_beam, Eep):
    
    denom = 4.0 * E_beam * Eep
    
    # if denom <= 0.0: return None
    
    sin2half = q2 / denom
    
    # if sin2half < 1.0 or sin2half > 0.0: return None
    
    return 2.0 * np.arcsin(np.sqrt(sin2half))

def Calc_kin(Ebeam, Eprime, q2, theta_e, forward=True):
    """
    Compute the hadronic momentum transfer t for e + p -> e' + π+ + n
    using two-body exclusive kinematics.
    forward=True gives t_min (pion emitted along q-vector).

    """
    # Energy and momentum transfer
    omega = Ebeam - Eprime
    q_abs = np.sqrt(q2 + omega**2)
    
    q_par = Ebeam - Eprime*np.cos(theta_e)
    q_perp = -1*Eprime*np.sin(theta_e)
    q_abs_comp = np.sqrt(q_par**2 + q_perp**2)

    q_vec = np.array([q_perp, 0, q_par])
    
    # Invariant mass of the hadronic system
    W2 = m_p**2 + 2*m_p*omega - q2
    if W2 <= (m_p + m_pi)**2:        
        raise ValueError("Unphysical kinematics: W too small")
    W = np.sqrt(W2)
    
    # CM pion momentum and energy
    ppi_star = np.sqrt((W2 - (m_p + m_pi)**2)*(W2 - (m_p - m_pi)**2)) / (2*W)
    Epi_star = np.sqrt(ppi_star**2 + m_pi**2)
    
    # Boost to lab frame
    beta = q_abs / (m_p + omega)
    gamma = (m_p + omega) / W
    
    # For forward pion emission (θ* = 0) → minimal |t|
    cos_theta_pq = 1.0 if forward else -1.0
    Epi_lab = gamma * (Epi_star + beta * ppi_star * cos_theta_pq)
    ppi_lab = np.sqrt(Epi_lab**2 - m_pi**2)
    
    #Pion angle and vector
    cos_q = q_par / q_abs_comp if q_abs_comp != 0 else 1.0
    theta_pi = (np.arccos(cos_q) + np.arccos(cos_theta_pq)) * 180 / np.pi
    p_pi_vec = np.array([ppi_lab * np.sin(np.deg2rad(theta_pi)), 0, ppi_lab * np.cos(np.deg2rad(theta_pi))])
    
    k_pi = np.sqrt(max(0.0, ppi_lab**2 + q_abs_comp**2 - 2.0 * ppi_lab * q_abs * cos_theta_pq))

    # Missing energy, momentum, and mass calculations of the nucleon
    E_x = omega - Epi_lab + m_A
    P_x_vec = q_vec - p_pi_vec; P_x2 = np.dot(P_x_vec, P_x_vec); P_x = np.sqrt(P_x2)
    M_x2 = E_x**2 - P_x2; M_x = np.sqrt(max(0.0, M_x2))
    
    # Compute t (minimal value)
    t = -q2 + m_pi**2 - 2 * (omega * Epi_lab - q_abs * ppi_lab * cos_theta_pq)
    
    return t, ppi_lab, theta_pi, k_pi, W, E_x, P_x, M_x

def Q_squared(Ee, Eep, theta): # theta is in radians
    
    return 4*Ee*Eep*np.sin(theta/2)**2

def XBjorken(Ee, Eep, theta): # elastic (jlab)
    
    return Q_squared(Ee, Eep, theta)/(2*m_p*(Ee-Eep))

# RESULTS

results = []

q2_results = []
xb_results = []
theta_e_results = []
p_pi_results = []
theta_pi_results = []
k_pi_results = []
t_results = []
W_results = []
Ex_results = []
Px_results = []
Mx_results = []

# INPUT 

Ebeam_results = []
Eprime_results = []


for Ebeam in E_beam:
    #print(Ebeam)
    for Q2_val in Q2:
        #print(Q2_val)
        for Eprime in E_prime:
            #print(Eprime)
            try:
                #print(Eprime)
                
                theta_e = Theta_e(Q2_val, Ebeam, Eprime)
                
                q2_val = Q_squared(Ebeam, Eprime, theta_e)
                
                xb_val = XBjorken(Ebeam, Eprime, theta_e)
                
                theta_e *= 180/np.pi # radians conversion
                
                t, p_pi, theta_pi, k_pi, W, E_x, P_x, M_x = Calc_kin(Ebeam, Eprime, Q2_val, theta_e, forward=True)

                t *= -1
                
                if np.isnan(t): continue
                
                if np.abs(t) > 1.0: continue
                
                #results
                Ebeam_results.append(Ebeam); Eprime_results.append(Eprime)
                q2_results.append(q2_val); xb_results.append(xb_val)
                theta_e_results.append(theta_e); theta_pi_results.append(theta_pi)
                t_results.append(t); W_results.append(W)
                p_pi_results.append(p_pi); k_pi_results.append(k_pi)
                Ex_results.append(E_x); Px_results.append(P_x)
                Mx_results.append(M_x); 

                results.append([Q2_val, xb_val, t, Ebeam, theta_e, Eprime, theta_pi, p_pi, k_pi, W, E_x, P_x, M_x])
                
            except Exception:
                continue

df = pd.DataFrame(results, columns=["Q2", "Bjorken X", "t", "Ebeam", "Theta_e", "Eprime", "Theta_p", "p_pi", "k_pi", "W", "Em", "Pm", "Mm"])
df.to_csv("t_scan_results.csv", index=False, float_format="%.6f")
print("Saved results to t_scan_results.csv with", len(df), "rows.")

Ebeam_results = np.array(Ebeam_results)
Eprime_results = np.array(Eprime_results)

nu_results = Ebeam_results - Eprime_results

q2_results = np.array(q2_results)
xb_results = np.array(xb_results)
theta_e_results = np.array(theta_e_results)
p_pi_results = np.array(p_pi_results)
theta_pi_results = np.array(theta_pi_results)
k_pi_results = np.array(k_pi_results)
t_results = np.array(t_results)
W_results = np.array(W_results)
Em_results = np.array(Ex_results)
Pm_results = np.array(Px_results)
Mm_results = np.array(Mx_results)

# PLOTTING with constraints

mask = (xb_results > 0) & (xb_results < 0.8) & (theta_e_results > 12) & (theta_e_results < 90)
massk = (Mm_results >= 1e-6) 

# fixed constraints
xb_fixed = (xb_results > 0.45) & (xb_results < 0.55) & mask     # x ~ 0.5
t_fixed = (t_results > 0.4) & (t_results < 0.45) & mask         # t ~ 0.4 - 0.45

# Labelling Ebeam, checking if range or fixed.

Ebeam_name = E_beam[0]

if len(E_beam) != 1:
    
    Ebeam_name = [E_beam[0], E_beam[-1]]
    
label_xb = (r'$x_b = 0.5$')

label_xb_t = (r'$x_b = 0.5$' '\n' 
              r'$t = -0.4$ (GeV/c)$^2$')

def Label(label):
    plt.text(0.98, 0.98, label,
             transform=plt.gca().transAxes,
             va='top', ha='right',
             bbox=dict(facecolor='white',
                       alpha=0.8,edgecolor='black'))

# Variable dictionaries

results = {"nu": nu_results,
           "Q2": q2_results,
           "W": W_results,
           "xb": xb_results,
           "Es": Eprime_results,
           "theta_e": theta_e_results,
           "t": t_results,
           "p_pi": p_pi_results,
           "theta_pi": theta_pi_results,
           "k_pi": k_pi_results,
           "Em": Em_results,
           "Pm": Pm_results,
           "Mm": Mm_results}

# dictionary for plotting phase space graphs
plotPS = [("xb", "Q2", "theta_e", mask),
         ("theta_pi", "p_pi", "Es", mask),
         ("theta_e", "Q2", "Es", xb_fixed),
         ("t", "p_pi", "Es", xb_fixed),
         ("t", "theta_pi", "Es", xb_fixed),
         ("theta_pi", "k_pi", "Es", xb_fixed),
         ("theta_e", "Q2", "t", xb_fixed),
         ("theta_pi", "Q2", "t", xb_fixed),
         ("W", "Q2", "t", xb_fixed),
         ("Mm", "Pm", "Es", mask & massk),
         ("Mm", "Em", "Es", mask & massk),
         ("Em", "Pm", "Es", mask & massk),]

plotH = [("Mm", 100, mask & massk)]

plot2H = [("Q2", "W", 100, mask),
         ("Mm", "Pm", 100, mask)]

# RESULTS

for xkey, ykey, zkey, mask in plotPS:

    cth.scatter(results[xkey], results[ykey], results[zkey], mask)
    cth.format(cth.labels[xkey], cth.labels[ykey], cth.labels[zkey],
                fr'Phase Space for $E_b=$ {Ebeam_name} GeV')

for key, binsize, mask in plotH:

    cth.hist(results[key], bins=binsize, weights=None, mask=mask)
    cth.format(cth.labels[key], ylabel='Counts', colorbar=None, title=
               fr'Counts Graphs for $E_b=$ {Ebeam_name} GeV')

for xkey, ykey, binsize, mask in plot2H:

    cth.hist2D(results[xkey], results[ykey], bins=binsize, weights=None, mask=mask)
    cth.format(cth.labels[xkey], cth.labels[ykey], colorbar=None, title=
                fr'Counts Graphs for $E_b=$ {Ebeam_name} GeV')

plt.show()