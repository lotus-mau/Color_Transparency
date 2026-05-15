"""
Helper file used to import certain redundant functions 
that are helpful when conducting research for CT

primarily uses matplotlib.pyplot for plotting functions
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata

## CONSTANTS

q_e = 1.602e-19 * 1000  # C -> mC           elementary charge
N_A = 6.022e23          # 1/mol             Avogadro's number

m_p  = 0.9382720813     # proton mass
m_n = 0.93956563        # neutron mass
m_pi = 0.139570611      # charged pion mass

A = 12
Z = 6

q2tags = {
    5.0: 'q5',
    6.5: 'q6p5',
    7.5: 'q7p5',
    8.5: 'q8p5'
}

# SPECS (figure out later how to extract from hist file)

ngen = 10000        # number of generated events

normfac = { # from hist file
            "q5":   {"1pi": 0.144897E+07,
                     "2pi": 0.144897E+07,
                     "norad": 0.144393E+08},
            "q6p5": {"1pi": 0.117848E+08,
                     "2pi": 0.117848E+08,
                     "norad": 1},
            "q7p5": {"1pi": 0.966871E+07,
                     "2pi": 0.966871E+07,
                     "norad": 1},
            "q8p5": {"1pi": 0.799880E+07,
                     "2pi": 0.799880E+07,
                     "norad": 1}
            }

length = {  # cm
            "q5":   {"1pi": 0.995858E+00,
                     "2pi": 0.995858E+00,
                     "norad": 0.995858E+00},
            "q6p5": {"1pi": 0.995858E+00,
                     "2pi": 0.995858E+00,
                     "norad": 0.995858E+00},
            "q7p5": {"1pi": 0.995858E+00,
                     "2pi": 0.995858E+00,
                     "norad": 0.995858E+00},
            "q8p5": {"1pi": 0.995858E+00,
                     "2pi": 0.995858E+00,
                     "norad": 0.995858E+00}
            }

density = { # g/cm^3
            "q5":   {"1pi": 0.169000E+01,
                     "2pi": 0.169000E+01,
                     "norad": 0.169000E+01},
            "q6p5": {"1pi": 0.169000E+01,
                     "2pi": 0.169000E+01,
                     "norad": 0.169000E+01},
            "q7p5": {"1pi": 0.169000E+01,
                     "2pi": 0.169000E+01,
                     "norad": 0.169000E+01},
            "q8p5": {"1pi": 0.169000E+01,
                     "2pi": 0.169000E+01,
                     "norad": 0.169000E+01}
            }


current = 40 / 1000     # muA -> mA : mC/s

## HELPERS

def format(xlabel, ylabel, colorbar, title):

    plt.xlabel(xlabel); plt.ylabel(ylabel); plt.title(title)

    if colorbar is not None:
        plt.colorbar(label=colorbar)

def hist(x, bins, weights, mask, type):
    
    if mask is not None:
        
        x = x[mask]; 
        
        if weights is not None:
            
            weights = weights[mask]

    plt.hist(x, bins=bins, weights=weights, histtype=type)

def hist2D(x, y, bins, weights, mask):

    # x = np.asarray(x)
    # y = np.asarray(y)
    # weights = np.asarray(weights)

    if mask is not None:
        
        x = x[mask]; y = y[mask]; 
        
        if weights is not None:
            
            weights=weights[mask]

    plt.hist2d(x, y, bins=bins, weights=weights 
               , norm=LogNorm()
               )

def scatter(x, y, z, mask, fig):

    if mask is not None:

        x = x[mask]; y = y[mask]; z = z[mask]

    plt.scatter(x, y, c=z, marker='s', 
                s=(450./fig.dpi)**2, edgecolors="None", cmap='viridis')
    
def contour(x, y, z, mask):
    
    if mask is not None:

        x = x[mask]; y = y[mask]; z = z[mask]
    
    x_grid = np.linspace(x.min(), x.max(), 200)
    y_grid = np.linspace(y.min(), y.max(), 200)
    X, Y = np.meshgrid(x_grid, y_grid)
    
    z_grid = griddata((x, y), z, (X, Y), method='linear')
    
    plt.contourf(X, Y, z_grid, levels=50)

def label(label):
    plt.text(0.98, 0.98, label,
             transform=plt.gca().transAxes,
             va='top', ha='right',
             bbox=dict(facecolor='white',
                       alpha=0.8,edgecolor='black'))
    
def savefig(folder, name):

    os.makedirs(folder, exist_ok=True)

    path = os.path.join(folder, name)

    plt.savefig(path, dpi=300, bbox_inches='tight')

# DEFINITIONS

def calclum(current, density, length, A):

    lum_B = current / q_e                     # number of electrons per second
    lum_T = density * length / A * N_A        # number of particles (nucleon per nucleus) in unit area

    return lum_B * lum_T

def luminosity(q2):
    return {
            "1pi": calclum(current, density[q2tags[q2]]['1pi'], length[q2tags[q2]]['1pi'], A),
            "2pi": calclum(current, density[q2tags[q2]]['2pi'], length[q2tags[q2]]['2pi'], A),
            "norad": calclum(current, density[q2tags[q2]]['norad'], length[q2tags[q2]]['norad'], A)
            }

## USEFUL VARIABLES

vars = [
    "q",        # magnitude of q vector 
    "nu",       # Energy transfer := E - E'
    "Q2",       # virtual photon momentum transfer
    "W",        # invariant rest mass of system
    "epsilon",  # ?
    "Eprime",   # scattered electron energy
    "theta_e",  # scattered electron angle
    "Em",       # reconstructed missing energy
    "Pm",       # reconstructed missing momentum
    "k_pi",     # pion 3-momentum
    "p_pi",     # pion 4-momentum
    "theta_pi", # pion angle
    "thetapq",  # angle between pion and q
    "phipq",    # angle between pion and q atomic planes
    "mmnuc",    # reconstructed missing nuclear mass
    "phad",     # momentum of hadronic system ?
    "t",        # mandelstam t, momentum transfer of hadronic system
    "pmpar",    # ?
    "pmper",    # ?
    "pmoop",    # ?
    "radphot",  # ?
    "pfermi"    # fermi momentum ?
]

labels = {
    "q":        r"$|\vec{q}|\ \mathrm{(GeV/c)}$",
    "nu":       r"$\nu\ \mathrm{(GeV/c)}$",
    "Q2":       r"$Q^2\ \mathrm{(GeV/c)^2}$",
    "W":        r"$W\ \mathrm{(GeV)}$",
    "xb":       r"$x_b$",
    "epsilon":  r"$\epsilon$",
    "Eprime":   r"$E'\ \mathrm{(GeV)}$",
    "theta_e":  r"$\theta_{\mathrm{e}}\ \mathrm{(\deg)}$",
    "Em":       r"$E_m\ \mathrm{(GeV)}$",
    "Pm":       r"$P_m\ \mathrm{(GeV/c)}$",
    "k_pi":     r"$k_\pi\ \mathrm{(GeV/c)}$",
    "p_pi":     r"$p_\pi\ \mathrm{(GeV/c)}$",
    "theta_pi": r"$\theta_\pi\ \mathrm{(\deg)}$",
    "thetapq":  r"$\theta_{pq}\ \mathrm{(rad)}$",
    "phipq":    r"$\phi_{pq}\ \mathrm{(rad)}$",
    "Mm":       r"$M_{\mathrm{m}}\ \mathrm{(GeV)}$",
    "mmnuc":    r"$M_{\mathrm{m}}^{\mathrm{nuc}}\ \mathrm{(GeV)}$",
    "phad":     r"$p_{\mathrm{had}}\ \mathrm{(GeV/c)}$",
    "t":        r"$t\ \mathrm{(GeV^2)}$",
    "pmpar":    r"$p_{m}^{\parallel}\ \mathrm{(GeV/c)}$",
    "pmper":    r"$p_{m}^{\perp}\ \mathrm{(GeV/c)}$",
    "pmoop":    r"$p_{m,\mathrm{oop}}\ \mathrm{(GeV/c)}$",
    "radphot":  r"",
    "pfermi":   r"$p_{\mathrm{f}}\ \mathrm{(GeV/c)}$",
    "Counts_mC": r"$Counts/mC$",
    "Counts_s":  r"$Counts/s$"
}