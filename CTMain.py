"""
File to compile calculations.
"""

import numpy as np
import matplotlib.pyplot as plt

import CTHelp as cth
import CTSim as cts
import CTKin as ctk
import CTKinTable as ctkt

# INPUTS

inputKT = [5.0,       # Q2
           11.0,      # Ebeam
           -0.527     # t_target
           ]

inputK = [np.array([11.0]),              # E_beam       
          np.arange(4.0, 5.0, 0.01),     # Q2
          12,                            # A
          6                              # 6
          ]

# Labelling Ebeam, checking if range or fixed.

Ebeam = inputK[0]; Ebeam_name = Ebeam[0]

if len(Ebeam) != 1:
    
    Ebeam_name = [Ebeam[0], Ebeam[-1]]
    
label_xb = (r'$x_b = 0.5$')

label_xb_t = (r'$x_b = 0.5$' '\n' 
            r'$t = -0.4$ (GeV/c)$^2$')

# RESULTS DICTIONARY

resultsSIM = {"1pi": cts.main("pionCT-sim/pion_q5_og.root"),
              "2pi": cts.main("pionCT-sim/pion_q5_2pi.root"),
              "norad": cts.main("pionCT-sim/pion_q5_norad.root")}

# resultsKT           = ctkt.main(inputKT)
# resultsK, masksK    = ctk.main(inputK)

# KIN PLOTTING

# (xkey, ykey, zkey, mask)
plotPS = [("xb", "Q2", "theta_e", "detection"),
            ("theta_pi", "p_pi", "Eprime", "detection"),
            ("theta_e", "Q2", "Eprime", "fix_xb"),
            ("t", "p_pi", "Eprime", "fix_xb"),
            ("t", "theta_pi", "Eprime", "fix_xb"),
            ("theta_pi", "k_pi", "Eprime", "fix_xb"),
            ("theta_e", "Q2", "t", "fix_xb"),
            ("theta_pi", "Q2", "t", "fix_xb"),
            ("W", "Q2", "t", "fix_xb"),
            ("Mm", "Pm", "Eprime", "detection"),
            ("Mm", "Em", "Eprime", "detection"),
            ("Em", "Pm", "Eprime", "detection"),
            ("Mm", "Q2", "xb", "detection")]

# (xkey, binsize, mask)
plotH = [("Mm", 100, "detection")]

# (xkey, ykey, binsize, mask)
plot2H = [("Q2", "W", 100, "detection"),
          ("Mm", "Pm", 100, "detection"),
          ("Mm", "Q2", 100, "detection")]

# for xkey, ykey, zkey, mask in plotPS:

#     fig, _ = plt.subplots()
#     cth.scatter(resultsK[xkey], resultsK[ykey], resultsK[zkey], 
#                 masksK[mask], fig)
#     cth.format(cth.labels[xkey], cth.labels[ykey], cth.labels[zkey],
#                 fr'Phase Space for $E_b=$ {Ebeam_name} GeV')

# for key, binsize, mask in plotH:

#     plt.figure()
#     cth.hist(resultsK[key], bins=binsize, 
#              weights=None, mask=masksK[mask])
#     cth.format(cth.labels[key], ylabel='Counts', colorbar=None, title=
#             fr'Counts Graphs for $E_b=$ {Ebeam_name} GeV')

# for xkey, ykey, binsize, mask in plot2H:

#     plt.figure()
#     cth.hist2D(resultsK[xkey], resultsK[ykey], binsize, 
#                weights=None, mask=masksK[mask])
#     cth.format(cth.labels[xkey], cth.labels[ykey], colorbar=None, title=
#                 fr'Counts Graphs for $E_b=$ {Ebeam_name} GeV')
    
# SIM PLOTTING

weights = {"1pi": resultsSIM['1pi']['Weight'] * cth.normfac['1pi'] / cth.ngen,
           "2pi": resultsSIM['2pi']['Weight'] * cth.normfac['2pi'] / cth.ngen,
           "norad": resultsSIM['norad']['Weight'] * cth.normfac['norad'] / cth.ngen}

weights_rate = {"1pi": weights['1pi'] * cth.current,
                "2pi": weights['2pi'] * cth.current,
                "norad": weights['norad'] * cth.current}

weights_lum = {"1pi": weights_rate['1pi'] * cth.luminosity['1pi'],
               "2pi": weights_rate['2pi'] * cth.luminosity['2pi'],
               "norad": weights_rate['norad'] * cth.luminosity['norad']}


print('Weights Evaluated\n')

plt.figure()
cth.hist(resultsSIM['1pi']['mmnuc'], 100, weights_rate['1pi'], mask=None, type='step')
cth.hist(resultsSIM['2pi']['mmnuc'], 100, weights_rate['2pi'], mask=None, type='step')
cth.hist(resultsSIM['norad']['mmnuc'], 100, weights_rate['norad'], mask=None, type='step')
cth.format(cth.labels['mmnuc'], cth.labels['Counts_s'], colorbar=None, title=
           fr'Missing Mass Rates Graphs for 1pi, 2pi, and norad Processes')
plt.show()

# above plot is... interesting. ask Holly

print("Plots Created\n")

# plot:

# Mm spectrum for og (1pi) and 2pi, 
# so that 2pi tail corresponds to 1pi, 
# then no rad histo, which should be thinner
# then rho histo. which should be wider and lower
# All on one figure to see. 

# then, cut Mm below a certain threshold.