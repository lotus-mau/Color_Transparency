"""
Helper file used to import certain redundant functions 
that are helpful when conducting research for CT

primarily uses matplotlib.pyplot for plotting functions
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

## CONSTANTS

# constants
q_e = 1.602e-19 * 1000  # C -> mC           elementary charge
N_A = 6.022e23          # 1/mol             Avogadro's number

## DEFINITIONS

def format(xlabel, ylabel, colorbar, title):

    plt.xlabel(xlabel); plt.ylabel(ylabel); plt.title(title)

    if colorbar != None:
        plt.colorbar(label=colorbar)

def hist(x, bins, weights, mask):

    x_masked = x[mask]

    plt.figure()
    plt.hist(x_masked, bins=bins, weights=weights)

def hist2D(x, y, bins, weights, mask):

    x_masked = x[mask] ; y_masked = y[mask]

    plt.figure()
    plt.hist2d(x_masked, y_masked, bins=bins, weights=weights, norm=LogNorm())

def scatter(x, y, z, mask):

    x_masked = x[mask]; y_masked = y[mask]; z_masked = z[mask]
    
    fig, _ = plt.subplots()
    plt.scatter(x_masked, y_masked, c=z_masked, marker='s', 
                s=(450./fig.dpi)**2, edgecolors="None", cmap='viridis')

def calclum(current, density, length, A):

    L_B = current / q_e                     # number of electrons per second
    L_T = density * length / A * N_A        # number of particles (nucleon per nucleus) in unit area

    return L_B * L_T

## USEFUL VARIABLES

vars = [
    "q",        # magnitude of q vector 
    "nu",       # Energy transfer := E - E'
    "Q2",       # virtual photon momentum transfer
    "W",        # invariant rest mass of system
    "epsilon",  # ?
    "Es",       # scattered electron energy
    "theta_e",  # scattered electron angle
    "Em",       # reconstructed missing energy
    "Pm",       # reconstructed missing momentum
    "p_pi",     # pion momentum
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
    "Es":       r"$E_s\ \mathrm{(GeV)}$",
    "theta_e":  r"$\theta_{\mathrm{e}}\ \mathrm{(\deg)}$",
    "Em":       r"$E_m\ \mathrm{(GeV)}$",
    "Pm":       r"$P_m\ \mathrm{(GeV/c)}$",
    "k_pi":     r"$k_\pi\ \mathrm{(GeV/c)}$",
    "p_pi":     r"$p_\pi\ \mathrm{(GeV/c)}$",
    "theta_pi": r"$\theta_\pi\ \mathrm{(\deg)}$",
    "thetapq":  r"$\theta_{pq}\ \mathrm{(rad)}$",
    "phipq":    r"$\phi_{pq}\ \mathrm{(rad)}$",
    "Mm":       r"$M_{\mathrm{m}}\ \mathrm{(GeV)}$",
    "Mmnuc":    r"$M_{\mathrm{m}}\ \mathrm{(GeV)}$",
    "phad":     r"$p_{\mathrm{had}}\ \mathrm{(GeV/c)}$",
    "t":        r"$t\ \mathrm{(GeV^2)}$",
    "pmpar":    r"$p_{m}^{\parallel}\ \mathrm{(GeV/c)}$",
    "pmper":    r"$p_{m}^{\perp}\ \mathrm{(GeV/c)}$",
    "pmoop":    r"$p_{m,\mathrm{oop}}\ \mathrm{(GeV/c)}$",
    "radphot":  r"",
    "pfermi":   r"$p_{\mathrm{f}}\ \mathrm{(GeV/c)}$",
}