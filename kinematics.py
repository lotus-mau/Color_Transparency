# -*- coding: utf-8 -*-
"""
Created on Fri Oct  3 17:04:40 2025

@author: Lotus
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

plt.close('all')

mass_p = 0.938

def bjorken(Ee, Eep, theta): # elastic (jlab)
    
    return q2(Ee, Eep, theta)/(2*mass_p*(Ee-Eep))

def q2(Ee, Eep, theta):
    
    return 4*Ee*Eep*np.sin(theta/2)**2

def s(Ee, Ep):
    
    return 4 * Ee * Ep

def y(Ee, Eep):
    
    return 1 - Eep/Ee

# Jlab

Ee_j = 11       # 6, 11

Eep_j = np.arange(0.5, Ee_j, 0.01)

theta_j = np.radians(np.arange(1, 180, 1))

def Ee_scat(E, theta): # elastic (proton at rest)
    
    return E / (1 + (2*E/mass_p)*np.sin(theta/2)**2)

# Eep_j = Ee_scat(Ee_j, theta_vals)

Eep_j, theta_j = np.meshgrid(Eep_j, theta_j)

q2_j = q2(Ee_j, Eep_j, theta_j)

xb_j = bjorken(Ee_j, Eep_j, theta_j)

# PLOT

# Q^2
plt.figure()
cp1_j = plt.contourf(theta_j*180/np.pi, Eep_j, q2_j, levels=50)
plt.colorbar(cp1_j, label=r'$Q^2\ \mathrm{(GeV^2)}$')
plt.ylabel(r"Scattered electron energy $E'$ [GeV]")
plt.xlabel(r"Scattering angle $\theta_e$ [deg]")
plt.title(r"$Q^2$ contours for $E_e = {Ee_j}$ GeV")
plt.show()

# Bjorken X
plt.figure()
cp2_j = plt.contourf(theta_j*180/np.pi, Eep_j, xb_j, levels=np.linspace(0, 1, 50))
plt.xscale('linear')
plt.yscale('linear')
plt.colorbar(cp2_j, label=r"$x_{Bj}$")
plt.ylabel(r"Scattered electron energy $E'$ [GeV]")
plt.xlabel(r"Scattering angle $\theta_e$ [deg]")
plt.title(r"$x_{Bj}$ contours for $E_e = {Ee_j}$ GeV")

# EIC

Ee_eic = 18     # 5, 18
Ep_eic = 275    # 41, 275

s_eic = s(Ee_eic, Ep_eic)

Eep_eic = np.arange(0.5, Ee_eic, 0.01)

theta_eic = np.radians(np.arange(1, 180, 1))

Eep_eic, theta_eic = np.meshgrid(Eep_eic, theta_eic)

y_eic = y(Ee_eic, Eep_eic)

q2_eic = q2(Ee_eic, Eep_eic, theta_eic)

xb_eic = q2_eic / (s_eic * y_eic)

# PLOT

# Q^2
plt.figure()
cp1_eic = plt.contourf(theta_eic*180/np.pi, Eep_eic, q2_eic, levels=50)
plt.colorbar(cp1_eic, label=r'$Q^2\ \mathrm{(GeV^2)}$')
plt.ylabel(r"Scattered electron energy $E'$ [GeV]")
plt.xlabel(r"Scattering angle $\theta_e$ [deg]")
plt.title(r"$Q^2$ contours for $E_e = {Ee_{eic}}$ GeV") 

# Bjorken X
plt.figure()
cp2_eic = plt.contourf(theta_eic*180/np.pi, Eep_eic, xb_eic, levels=np.linspace(0, 1, 50))
plt.xscale('linear')
plt.yscale('linear')
plt.colorbar(cp2_eic, label=r"$x_{Bj}$")
plt.ylabel(r"Scattered electron energy $E'$ [GeV]")
plt.xlabel(r"Scattering angle $\theta_e$ [deg]")
plt.title(r"$x_{Bj}$ contours for $E_e = {Ee_{eic}}$ GeV")

# PLOT 

def plot(Ee, Eep, theta, title): # Q2 vs. XB
    
    xb = bjorken(Ee, Eep, theta); Q2 = q2(Ee, Eep, theta)
    
    mask = (xb > 0) & (xb < 1)
    
    xb_masked = xb[mask]; q2_masked = Q2[mask]; theta_masked = theta[mask]*180/np.pi
    
    xb_grid = np.linspace(xb_masked.min(), xb_masked.max(), 200)
    q2_grid = np.linspace(q2_masked.min(), q2_masked.max(), 200)
    XB, Q2 = np.meshgrid(xb_grid, q2_grid)
    
    Theta_grid = griddata((xb_masked, q2_masked), theta_masked, (XB, Q2), method='linear')
    
    # Plot
    plt.figure()
    plt.contourf(XB, Q2, Theta_grid, levels=50)
    plt.colorbar(label=r'$\theta_e$')
    plt.xlabel(r'$x_B$')
    plt.ylabel(r'$Q^2$')
    plt.title(fr'Phase Space for {title}')
    
    # Plot scatter
    plt.figure()
    plt.scatter(xb_masked, q2_masked, c=theta_masked, marker='.')
    plt.colorbar(label=r'$\theta_e$')
    plt.xlabel(r'$x_B$'); plt.ylabel(r'$Q^2$')
    plt.title(fr'Phase Space for {title}')
    
plot(Ee_j, Eep_j, theta_j, 'JLab')

plot(Ee_eic, Eep_eic, theta_eic, 'EIC')

def plot_J(Ee, Eep, theta): # Q2 and XB -> E' vs. theta
    
    Eep, theta = np.meshgrid(Eep, theta)

    Q2 = q2(Ee, Eep, theta)

    XB = bjorken(Ee, Eep, theta)
    
    mask = (XB > 0) & (XB < 1)
    
    xb_masked = XB[mask]; q2_masked = Q2[mask]; theta_masked = theta[mask]*180/np.pi