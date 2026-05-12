import uproot
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import CTHelp as cth

def main(file_root):
    plt.close('all')

    # Open the ROOT file
    file = uproot.open(file_root)

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

    vars = ["q",        # magnitude of q vector 
            "nu",       # Energy transfer := E - E'
            "Q2",       # virtual photon momentum transfer
            "W",        # invariant rest mass of system
            "epsilon",  # ?
            "Em",       # reconstructed missing energy
            "Pm",       # reconstructed missing momentum
            "thetapq",  # angle between pion and q
            "phipq",    # angle between pion and q atomic planes
            "missmass", # missing mass
            "Mhadron",  # hadron mass
            "mmnuc",    # reconstructed missing nuclear mass
            "phad",     # momentum of hadronic system ?
            "t",        # mandelstam t, momentum transfer of hadronic system
            "pmpar",    # ?
            "pmper",    # ?
            "pmoop",    # ?
            "radphot",  # ?
            "pfermi",   # fermi momentum ?
            "Weight"    # Monte Carlo event weight => (cross section x acceptance)
            ]
    kin = {v: tree[v].array() for v in vars}

    print("Simulation Process Finished\n")

    return kin

if __name__ == "__main__":

    # file name(s)
    file = "pionCT-sim/pion_q5_og.root"

    main(file)