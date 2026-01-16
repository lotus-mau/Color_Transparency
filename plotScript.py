import uproot
import matplotlib.pyplot as plt

plt.close('all')

# Open the ROOT file
file = uproot.open("pion_q5_og.root")

# Get the tree
tree = file["h10"]

# Generator information, from .hist file
ngen = 10000
normfac = 5383780

# Read the branch data
# HMS variables you can read: hsdelta, hsyptar, hsxptar, hsytar, hsxfp, hsxpfp, hsyfp, hsypfp
# SHMS variables you can read: ssdelta, ssyptar, ssxptar, ssytar, ssxfp, ssxpfp, ssyfp, ssypfp
# Other variables for cross section: Weight
# Kinematic variables: q, nu, Q2, W, epsilon*, Em, Pm, thetapq*, phipq*, mmnuc, phad, t, pmpar, pmper, pmoop, radphot, pfermi
# Note that the weight will give you units of counts/mC

# you can plot any of these variables by doing something like the following:
mmass = tree["mmnuc"].array()
weight = tree["Weight"].array()

plt.figure()
# Create a histogram
plt.hist(mmass, bins=100, weights=weight*normfac/ngen)
# Set labels
plt.xlabel("Missing mass [GeV/c]")
plt.ylabel("Counts/mC")
plt.title("1mC of beam on carbon target")
# Show the plot
plt.show()

# Testing different plots.

# momentum of pion
p_pi = tree["phad"].array()

plt.figure()
plt.hist(p_pi, bins=100, weights=weight*normfac/ngen)
plt.xlabel(r"$p_\pi$ [GeV/c]"); plt.ylabel(r"Counts/mC")
plt.title("1mC of beam on carbon target")

# t, hadronic momentum transfer
t = tree["t"].array()

plt.figure()
plt.hist(t, bins=100, weights=weight*normfac/ngen)
plt.xlabel(r"t [GeV/c]^2"); plt.ylabel(r"Counts/mC")
plt.title("1mC of beam on carbon target")

# Q2, virtual photon momentum transfer
Q2 = tree["Q2"].array()

plt.figure()
plt.hist(Q2, bins=100, weights=weight*normfac/ngen)
plt.xlabel(r"Q2 [GeV/c]^2"); plt.ylabel(r"Counts/mC")
plt.title("1mC of beam on carbon target")

# Em, reconstructed energy of mass

Em = tree["Em"].array()

plt.figure()
plt.hist(Em, bins=100, weights=weight*normfac/ngen)
plt.xlabel(r"Em [GeV]"); plt.ylabel(r"Counts/mC")
plt.title("1mC of beam on carbon target")

# pfermi, Fermi momentum 

pfermi = tree["pfermi"].array()

plt.figure()
plt.hist(pfermi, bins=100, weights=weight*normfac/ngen)
plt.xlabel(r"$p_f$ [GeV/c]"); plt.ylabel(r"Counts/mC")
plt.title("1mC of beam on carbon target")