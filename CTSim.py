import uproot
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

plt.close('all')

# Open the ROOT file
file = uproot.open("pion_q5_og.root")

# Get the tree
tree = file["h10"]

# Read the branch data

# HMS variables you can read: 
# hsdelta, hsyptar, hsxptar, hsytar, hsxfp, hsxpfp, hsyfp, hsypfp

# SHMS variables you can read: 
# ssdelta, ssyptar, ssxptar, ssytar, ssxfp, ssxpfp, ssyfp, ssypfp

# Other variables for cross section: Weight

# Kinematic variables: 
# q, nu, Q2, W, epsilon*, Em, Pm, thetapq*, phipq*, 
# mmnuc, phad, t, pmpar, pmper, pmoop, radphot, pfermi

# Note that the weight will give you units of counts/mC

# you can plot any of these variables by doing something like the following:

# Monte Carlo event weight (cross section x acceptance)
weight = tree["Weight"].array()     

# Generator information, from .hist file
ngen = 10000        # number of generated events
normfac = 5383780   # total accumulated charge (mC)

vars = ["q",        # magnitude of q vector 
        "nu",       # Energy transfer := E - E'
        "Q2",       # virtual photon momentum transfer
        "W",        # invariant rest mass of system
        "epsilon",  # ?
        "Em",       # reconstructed missing energy
        "Pm",       # reconstructed missing momentum
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
kin = {v: tree[v].array() for v in vars}

labels = {
    "q":        r"$|\vec{q}|\ \mathrm{(GeV/c)}$",
    "nu":       r"$\nu\ \mathrm{(GeV/c)}$",
    "Q2":       r"$Q^2\ \mathrm{(GeV/c)^2}$",
    "W":        r"$W\ \mathrm{(GeV)}$",
    "epsilon":  r"$\epsilon$",
    "Em":       r"$E_m\ \mathrm{(GeV)}$",
    "Pm":       r"$P_m\ \mathrm{(GeV/c)}$",
    "thetapq":  r"$\theta_{pq}\ \mathrm{(rad)}$",
    "phipq":    r"$\phi_{pq}\ \mathrm{(rad)}$",
    "mmnuc":    r"$M_{\mathrm{m}}\ \mathrm{(GeV)}$",
    "phad":     r"$p_{\mathrm{had}}\ \mathrm{(GeV/c)}$",
    "t":        r"$t\ \mathrm{(GeV^2)}$",
    "pmpar":    r"$p_{m}^{\parallel}\ \mathrm{(GeV/c)}$",
    "pmper":    r"$p_{m}^{\perp}\ \mathrm{(GeV/c)}$",
    "pmoop":    r"$p_{m,\mathrm{oop}}\ \mathrm{(GeV/c)}$",
    "radphot":  r"",
    "pfermi":   r"$p_{\mathrm{f}}\ \mathrm{(GeV/c)}$",
}

def pltformat(xlabel, ylabel, title):
    plt.xlabel(xlabel); plt.ylabel(ylabel); plt.title(title)

#storing histograms
h_figs = []

plot1D = ['q', 'nu', 'Q2', 
          'W', 'epsilon', 'Em', 
          'thetapq', 'mmnuc', 'phad', 
          't']

# Plotting each graph
def pltallhist(vars, binsize, weights, 
               ylabel, title):
    for key in vars:
        
        plt.figure()
        
        # Current is in mA (mC/s)
        h = plt.hist(kin[key], bins=binsize, weights=weights)
        pltformat(labels.get(key, key), ylabel, title)

        h_figs.append(h)

# constants
q_e = 1.602e-19 * 1000  # C -> mC           elementary charge
N_A = 6.022e23          # 1/mol             Avogadro's number

# specs
current = 40 / 1000     # A -> muA = mC/s   beam current
length = 0.294089       # cm                beam length
density = 0.226700e01   # g/cm^3            beam density  
A = 12                  #                   nucleon number  

weight1 = np.ones_like(weight)

# Plotting all Counts/mC graphs
#pltallhist(plot1D, 100, weight * normfac/ngen, 
#           'Counts/mC', '1mC of Current on Carbon target')

# Plotting all Counts/s graphs 
#pltallhist(plot1D, 100, weight * normfac / ngen * current, 
#           'Counts/s', 'Count Rates on Carbon target')

# Calculating luminosity:
def calclum(current, density, length, A):

    L_B = current / q_e                     # number of electrons per second
    L_T = density * length / A * N_A        # number of particles (nucleon per nucleus) in unit area

    return L_B * L_T

luminosity = calclum(current, density, length, A)

# Plotting 2D Histograms

plot2D = [('Q2','thetapq'),
          ('Q2','t'),
          ('Q2','epsilon'),
          ('phad','thetapq')]

def pltall2D(vars, binsize, weights, title):

    for ykey, xkey in vars:

        plt.figure()
        h = plt.hist2d(np.asarray(kin[xkey]), np.asarray(kin[ykey]), 
                       bins=binsize, weights=np.asarray(weights), 
                       norm=LogNorm())
        plt.colorbar(label=r'Counts/s')
        pltformat(labels.get(xkey, xkey), labels.get(ykey, ykey), title)

pltall2D(plot2D, 100, weight * normfac / ngen * current, 
         'Count Rates on Carbon target, comparison')


plt.show()