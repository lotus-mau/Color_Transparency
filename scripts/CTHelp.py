"""
Helper file used to import certain redundant functions 
that are helpful when conducting research for CT

primarily uses matplotlib.pyplot for plotting functions
"""

import re
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

q2tags = ['q5', 'q6p5', 'q7p5', 'q8p5']
q2val_tags = {5.0: q2tags[0], 6.5: q2tags[1], 7.5: q2tags[2], 8.5: q2tags[3]}

settings = ['1pi', '2pi']

targets = ['C', 'H']

ngen = 10000            # number of generated events

current = 40 / 1000     # muA -> mA : mC/s

def parse_hist(q2tag, key, target):

    if key == '1pi':

        setting = ''

    elif key == '2pi':

        setting = '_multipi'

    elif key == 'norad':

        setting = '_norad'

    else:
        
        print('\nNot a valid setting.\n')

    filename = f'pionCT-simc/pion_{q2tag}_{target}{setting}.hist'

    with open(filename) as f:
        text = f.read()

    def get_value(variable):
        match = re.search(fr"{variable}\s*=\s*([-\d.E+]+)", text)
        return float(match.group(1)) if match else None
    
    variables = ['Ntried',          # number of events tried before reaching ngen successes
                 'normfac',         # normalization factor calculated using luminosity/ntried*ngen
                 'length',          # length of target
                 'rho',             # density of target
                 'mass'             # target mass in GeV
                 ]
    
    parsed = {variable: get_value(variable) for variable in variables}

    match1 = re.search(
        r"mass\s*=\s*([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s*MeV",
        text
    )

    if match1:
        parsed['mass'] = float(match1.group(1)) / 1000.0            # MeV to GeV : /1000.
    
    return parsed

specs = {tag: {key: {target: parse_hist(tag, key, target)
                     for target in targets}
               for key in settings}
         for tag in q2tags}

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

def getTarget(target):

    if target == 'C':

        A = 12; Z = 6
    
    elif target == 'H':
        
        A = 1; Z = 1

    else:

        print("\nNot a valid target. \n")

    return A, Z

# DEFINITIONS

def calclum(current, density, length, A):

    lum_B = current / q_e                     # number of electrons per second
    lum_T = density * length / A * N_A        # number of particles (nucleon per nucleus) in unit area

    return lum_B * lum_T

def luminosity(q2, target):

    A, _ = getTarget(target)

    return {key: calclum(current, 
                         specs[q2val_tags[q2]][key][target]['rho'], 
                         specs[q2val_tags[q2]][key][target]['length'], 
                         A) 
            for key in settings}

## USEFUL VARIABLES

vars = [
        #'h10',          # ?
        'hsdelta',      # ?
        'hsyptar',      # ?
        'hsxptar',      # ?
        'hsytar',       # ?
        'hsxfp',        # ?
        'hsxpfp',       # ?
        'hsyfp',        # ?
        'hsypfp',       # ?
        'hsdeltai',     # ?
        'hsyptari',     # ?
        'hsxptari',     # ?
        'hsytari',      # ?
        'ssdelta',      # ?
        'ssyptar',      # ?
        'ssxptar',      # ?
        'ssytar',       # ?
        'ssxfp',        # ?
        'ssxpfp',       # ?
        'ssyfp',        # ?
        'ssypfp',       # ?
        'ssdeltai',     # ?
        'ssyptari',     # ?
        'ssxptari',     # ?
        'ssytari',      # ?
        'q',            # magnitude of q vector 
        'nu',           # Energy transfer := E - E'
        'Q2',           # virtual photon momentum transfer
        'W',            # invariant rest mass of system
        'epsilon',      # ?
        'Em',           # reconstructed missing energy
        'Pm',           # reconstructed missing momentum
        'thetapq',      # angle between pion and q
        'phipq',        # angle between pion and q atomic planes
        'missmass',     # missing mass (using M = M_p in mm equation)
        'mmnuc',        # reconstructed missing nuclear mass (using M = M_A)
        'phad',         # momentum of hadron (maybe total hadron momentum?)
        't',            # mandelstam t, momentum transfer of hadronic system
        'pmpar',        # ?
        'pmper',        # ?
        'pmoop',        # ?
        'fry',          # ?
        'radphot',      # ?
        'pfermi',       # fermi momentum
        'siglab',       # ?
        'sigcm',        # ?
        'Weight',       # Monte Carlo event weight => (cross section x acceptance)
        'decdist',      # ?
        'Mhadron',      # ?
        'pdotqhat',     # ?
        'Q2i',          # ?
        'Wi',           # ?
        'ti',           # ?
        'phipqi',       # ?
        'dpimm'         # missing mass of second pion and nucleus
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
    "MMa":      r"$M_{\mathrm{m}}^{\mathrm{A}}\ \mathrm{(GeV)}$",
    "mmnuc":    r"$M_{\mathrm{m}}^{\mathrm{nuc}}\ \mathrm{(GeV)}$",
    "dpimm":    r"$M_{\mathrm{m}}^{\mathrm{nuc}}\ \mathrm{(GeV)}$",
    "missmass": r"$M_{\mathrm{m}}\ \mathrm{(GeV)}$",
    "Mhadron":  r"$M_{\mathrm{had}}\ \mathrm{(GeV)}$",
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