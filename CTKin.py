# -*- coding: utf-8 -*-
"""
Created on Sun Oct 26 19:31:09 2025

@author: Lotus
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import griddata

plt.close('all')

# CONSTANTS

m_p  = 0.9382720813     # proton mass
m_n = 0.93956563        # neutron mass
m_pi = 0.139570611      # charged pion mass

# INPUTS

A = 12; Z = 6           # nucleon and atomic number : Carbon-12
m_A = Z*m_p + (A-Z)*m_n # target mass in GeV/c^2

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
    if W2 < (m_n + m_pi)**2:
        raise ValueError("Unphysical kinematics: W too small")
    W = np.sqrt(W2)
    
    # CM pion momentum and energy
    ppi_star = np.sqrt((W2 - (m_n + m_pi)**2)*(W2 - (m_n - m_pi)**2)) / (2*W)
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

    # Missing energy, momentum, and mass calculations
    E_x = Ebeam - Eprime + m_p - Epi_lab
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
                Mx_results.append(M_x)

                results.append([Q2_val, xb_val, t, Ebeam, theta_e, Eprime, theta_pi, p_pi, k_pi, W, E_x, P_x, M_x])
                
            except Exception:
                continue

df = pd.DataFrame(results, columns=["Q2", "Bjorken X", "t", "Ebeam", "Theta_e", "Eprime", "Theta_p", "p_pi", "k_pi", "W", "Em", "Pm", "Mm"])
df.to_csv("t_scan_results.csv", index=False, float_format="%.6f")
print("Saved results to t_scan_results.csv with", len(df), "rows.")

Ebeam_results = np.array(Ebeam_results)
Eprime_results = np.array(Eprime_results)

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

# PLOTTING HELPERS

def plot_contour(x, y, z, mask):
    
    x_masked = x[mask]; y_masked = y[mask]; z_masked = z[mask]
    
    x_grid = np.linspace(x_masked.min(), x_masked.max(), 200)
    y_grid = np.linspace(y_masked.min(), y_masked.max(), 200)
    X, Y = np.meshgrid(x_grid, y_grid)
    
    z_grid = griddata((x_masked, y_masked), z_masked, (X, Y), method='linear')
    
    # Plot
    plt.figure()
    plt.contourf(X, Y, z_grid, levels=50)
    
def plot_scatter(x, y, z, mask):
    
    x_masked = x[mask]; y_masked = y[mask]; z_masked = z[mask]
    
    # Plot
    #plt.figure()
    #plt.scatter(x_masked, y_masked, c=z_masked)
    
    fig, ax = plt.subplots()
    plt.scatter(x_masked, y_masked, c=z_masked, marker='s', 
                s=(450./fig.dpi)**2, edgecolors="None", cmap='bone')
    
def plot_format(xlabel, ylabel, colorbar, title):
    
    plt.xlabel(xlabel); plt.ylabel(ylabel); plt.colorbar(label=colorbar); plt.title(title)

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


# plot_contour(xb_results, q2_results, theta_e_results, xb_mask)

plot_scatter(xb_results, q2_results, theta_e_results, mask)
plot_format(r'$x_b$', r'$Q^2$ (GeV/c)$^2$', r'$\theta_e$ (deg)', 
            fr'($Q^2$, $x_b$, $\theta_e$) Phase Space for $E_b=$ {Ebeam_name} GeV')

plot_scatter(theta_pi_results, p_pi_results, Eprime_results, mask)
plot_format(r'$\theta_\pi$ (deg)', r'$p_\pi$ (GeV/c)', r'$E_s$ (GeV)', 
            fr'($p_\pi$, $\theta_\pi$, $E_s$) Phase Space for $E_b=$ {Ebeam_name} GeV')

plot_scatter(theta_e_results, q2_results, Eprime_results, xb_fixed)
plot_format(r'$\theta_e$ (deg)', r'$Q^2$ (GeV/c)$^2$', r'$E_s$ (GeV)', 
            fr'($Q^2$, $\theta_e$, $E_s$) Phase Space for $E_b=$ {Ebeam_name} GeV')
Label(label_xb)

plot_scatter(t_results, p_pi_results, Eprime_results, xb_fixed)
plot_format(r'$-t$ (GeV/c)$^2$', r'$p_\pi$ (GeV/c)', r'$E_s$ (GeV)', 
            fr'($p_\pi$, $-t$, $E_s$) Phase Space for $E_b=$ {Ebeam_name} GeV')
Label(label_xb)

plot_scatter(t_results, theta_pi_results, Eprime_results, xb_fixed)
plot_format(r'$-t$ (GeV/c)$^2$', r'$\theta_\pi$ (deg)', r'$E_s$ (GeV)', 
            fr'($\theta_\pi$, $-t$, $E_s$) Phase Space for $E_b=$ {Ebeam_name} GeV')
Label(label_xb)

plot_scatter(theta_pi_results, k_pi_results, Eprime_results, xb_fixed)
plot_format(r'$\theta_\pi$ (deg)', r'$k_\pi$ (GeV/c)', r'$E_s$ (GeV)',
            fr'($k_\pi$, $\theta_\pi$, $E_s$) Phase Space for $E_b=$ {Ebeam_name} GeV')
Label(label_xb)

plot_scatter(theta_pi_results, p_pi_results, Eprime_results, xb_fixed & t_fixed)
plot_format(r'$\theta_\pi$ (deg)', r'$p_\pi$ (GeV/c)', r'$E_s$ (GeV)', 
            fr'($p_\pi$, $\theta_\pi$, $E_s$) Phase Space for $E_b=$ {Ebeam_name} GeV')
Label(label_xb_t)

plot_scatter(theta_e_results, q2_results, Eprime_results, xb_fixed & t_fixed)
plot_format(r'$\theta_e$ (deg)', r'$Q^2$ (GeV/c)$^2$', r'$E_s$ (GeV)', 
            fr'($Q^2$, $\theta_e$, $E_s$) Phase Space for $E_b=$ {Ebeam_name} GeV')
Label(label_xb_t)

plot_scatter(theta_pi_results, q2_results, Eprime_results, xb_fixed & t_fixed)
plot_format(r'$\theta_\pi$ (deg)', r'$Q^2$ (GeV/c)$^2$', r'$E_s$ (GeV)', 
            fr'($Q^2$, $\theta_\pi$, $E_s$) Phase Space for $E_b=$ {Ebeam_name} GeV')
Label(label_xb_t)

plot_scatter(W_results, q2_results, Eprime_results, mask)
plot_format(r'$W$ (GeV/c)$^2$', r'$Q^2$ (GeV/c)$^2$', r'$E_s$ (GeV)',
            fr'($Q^2$, $W$, $E_s) Phase Space for $E_b=$ {Ebeam_name} GeV')

plot_scatter(Mm_results, Pm_results, Eprime_results, mask & massk)
plot_format(r'$M_m$ (GeV/c$^2$)', r'$P_m$ (GeV/c)', r'$E_s$ (GeV)',
            fr'($P_m$, $M_m$, $E_s$) Phase Space for $E_b=$ {Ebeam_name} GeV')

plot_scatter(Mm_results, Em_results, Eprime_results, mask & massk)

plot_scatter(Em_results, Pm_results, Eprime_results, mask & massk)

Mm_results = Mm_results[Mm_results >= 1e-6]
plt.figure()
plt.hist(Mm_results, bins=100)

plt.show()