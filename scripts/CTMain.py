"""
File to compile calculations.
"""

import sys
import time
import threading 

import numpy as np
import matplotlib.pyplot as plt

import CTHelp as cth
import CTSim as cts
import CTKin as ctk
import CTKinTable as ctkt

# sanity checks

_running = True

def heartbeat():
    while _running:
        print(".", end="", flush=True)
        time.sleep(1)

hb_thread = threading.Thread(target=heartbeat, daemon=True)
hb_thread.start()

t0 = time.time()

# INPUTS

if len(sys.argv) < 2:
    print("Usage: python CTMain.py <Q2:float>")
    sys.exit(1)

Q2 = float(sys.argv[1])
beam = 10.7

q2tags = {
    5.0: "q5",
    6.5: "q6p5",
    7.5: "q7p5",
    8.5: "q8p5",
}

inputKT = [Q2,       
           beam,        # Ebeam
           -0.527       # t_target
           ]

inputK = [np.array([beam]),                     # E_beam       
          np.arange(Q2-0.5, Q2+0.5, 0.01),  # Q2 range
          cth.A,                                # A
          cth.Z                                 # Z
          ]

# Labelling Ebeam, checking if range or fixed.

Ebeam = inputK[0]; Ebeam_name = Ebeam[0]

if len(Ebeam) != 1:
    
    Ebeam_name = [Ebeam[0], Ebeam[-1]]
    
label_xb = (r'$x_b = 0.5$')

label_xb_t = (r'$x_b = 0.5$' '\n' 
              r'$t = -0.4$ (GeV/c)$^2$')

# RESULTS DICTIONARY

resultsSIM = {
    "1pi": cts.main(f"pionCT-sim/pion_{q2tags[Q2]}_og.root"),
    "2pi": cts.main(f"pionCT-sim/pion_{q2tags[Q2]}_2pi.root"),
    "norad": cts.main(f"pionCT-sim/pion_{q2tags[Q2]}_norad.root"),
}

#%% Kinematic Table Calculations 
tKT = time.time()
resultsKT = ctkt.main(inputKT)
print(f"CTKinTable finished: {time.time() - tKT:.2f} s")

#%% Kinematic Phase Space Calculations
tK = time.time()
resultsK, masksK = ctk.main(inputK)
print(f"CTKin finished: {time.time() - tK:.2f} s")

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
            ("Mm", "Pm", "Eprime", "physical"),
            ("Mm", "Em", "Eprime", "detection"),
            ("Em", "Pm", "Eprime", "detection"),
            ("Mm", "Q2", "Eprime", "physical")
            ]

# (xkey, binsize, mask)
plotH = [("Mm", 100, "physical")]

# (xkey, ykey, binsize, mask)
plot2H = [("Mm", "Pm", 100, "physical"),
          ("Mm", "Q2", 100, "physical")]

tPS = time.time()
for xkey, ykey, zkey, mask in plotPS:

    fig, _ = plt.subplots()
    cth.scatter(resultsK[xkey], resultsK[ykey], resultsK[zkey], 
                masksK[mask], fig)
    cth.format(cth.labels[xkey], cth.labels[ykey], cth.labels[zkey],
                fr'Phase Space for $E_b=$ {Ebeam_name} GeV')
    cth.savefig(f'figures/Q2={Q2}', f'PS_{xkey}_{ykey}_{mask}')
print(f"PS plots created: {time.time() - tPS:.2f} s")

tH = time.time()
for key, binsize, mask in plotH:

    plt.figure()
    cth.hist(x=resultsK[key], bins=binsize, 
             weights=None, mask=masksK[mask], type='bar')
    cth.format(cth.labels[key], ylabel='Counts', colorbar=None, title=
            fr'Counts Graphs for $E_b=$ {Ebeam_name} GeV')
    cth.savefig(f'figures/Q2={Q2}', f'histo_{key}_{mask}')
print(f"Histograms created: {time.time() - tH:.2f} s")

t2H = time.time()
for xkey, ykey, binsize, mask in plot2H:

    plt.figure()
    cth.hist2D(resultsK[xkey], resultsK[ykey], binsize, 
               weights=None, mask=masksK[mask])
    cth.format(cth.labels[xkey], cth.labels[ykey], colorbar=None, title=
                fr'Counts Graphs for $E_b=$ {Ebeam_name} GeV')
    cth.savefig(f'figures/Q2={Q2}', f'histo2D_{xkey}_{ykey}_{mask}')
print(f"2D Histograms created: {time.time() - t2H:.2f} s")
    
# SIM PLOTTING

# coefficients
weights = {"1pi": resultsSIM['1pi']['Weight'] * cth.normfac['1pi'] / cth.ngen,
           "2pi": resultsSIM['2pi']['Weight'] * cth.normfac['2pi'] / cth.ngen,
           "norad": resultsSIM['norad']['Weight'] * cth.normfac['norad'] / cth.ngen}

weight_norm = 0.243

weights_rate = {"1pi": weights['1pi'] * cth.current,
                "2pi": weights['2pi'] * cth.current,
                "norad": weights['norad'] * cth.current * weight_norm}

weights_lum = {"1pi": weights_rate['1pi'] * cth.luminosity['1pi'],
               "2pi": weights_rate['2pi'] * cth.luminosity['2pi'],
               "norad": weights_rate['norad'] * cth.luminosity['norad']}

tSIM = time.time()
plt.figure()
for key in resultsSIM:
    cth.hist(resultsSIM[key]['mmnuc'], 100, weights_rate[key], mask=None, type='step')
cth.format(cth.labels['mmnuc'], cth.labels['Counts_s'], colorbar=None, title=
           fr'Missing Mass Rates Graphs for 1pi, 2pi, and norad Processes')

# above plot is... interesting. ask Holly
# answer: norad was just to compare radiative effects on the reaction, 
# plus, plot makes sense relatively. 

# (xkey, ykey, binsize, mask)
plotSIM = [("mmnuc", "Pm", 100, weights_rate["1pi"]),
        ("mmnuc", "Q2", 100, weights_rate["1pi"])]

for xkey, ykey, binsize, weight in plotSIM:

    x = np.asarray(resultsSIM['1pi'][xkey])
    y = np.asarray(resultsSIM['1pi'][ykey])
    w = np.asarray(weight)

    plt.figure()
    cth.hist2D(x, y, binsize, 
               weights=w, mask=None)
    cth.format(cth.labels[xkey], cth.labels[ykey], colorbar='Counts/s', title=
                fr'Counts Graphs for $E_b=$ {Ebeam_name} GeV')

plt.figure()
cth.hist(resultsSIM['1pi']['mmnuc'], 100, weights_rate['1pi'], mask=None, type='step')
cth.hist(resultsSIM['2pi']['mmnuc'], 100, weights_rate['2pi'], mask=None, type='stepfilled')
cth.format(cth.labels['mmnuc'], cth.labels['Counts_s'], colorbar=None, title=
           fr'Missing Mass Rates Graphs for 1pi, 2pi, and norad Processes')

print(f"SIM plots created: {time.time() - tSIM:.2f} s")

# plot:

# Mm spectrum for og (1pi) and 2pi, 
# so that 2pi tail corresponds to 1pi, 
# then no rad histo, which should be thinner # not really for ct, more for understanding
# then rho histo. which should be wider and lower
# All on one figure to see. 

# then, cut Mm below a certain threshold.

#target diff:

density_diff = abs(cth.density['1pi'] - cth.density['norad']) / cth.density['1pi']
length_diff = abs(cth.length['1pi'] - cth.length['norad']) / cth.length['1pi']

print(f"\n Target Difference: \n Density: {density_diff:.3f} \n Length {length_diff:.3f}\n")
print(f"1pi to norad Peak normalization: {weight_norm:.3f}\n")

# table:
# place cuts around left tail of 2pi, 
# and calculate the contamination of 2pi related to 1pi. (fraction of 2pi from 1pi)
# cut on MM | 2pi rate | signal rate (subtracted)

cut_slider = np.arange(11.5, 12.5, 0.1)

cut_results = []
tot1pi_results = []
tot2pi_results = []
contamination_results = []

for cut in cut_slider: # Missing mass

    mm_cut = {'1pi': (resultsSIM['1pi']['mmnuc'] < cut), 
              '2pi': (resultsSIM['2pi']['mmnuc'] < cut),
              'norad': (resultsSIM['norad']['mmnuc'] < cut)}
    
    h_tot = {'1pi': np.sum(weights['1pi'][mm_cut['1pi']]),
             '2pi': np.sum(weights['2pi'][mm_cut['2pi']])}
    
    frac = h_tot['2pi'] / h_tot['1pi']

    cut_results.append(cut)
    tot1pi_results.append(h_tot['1pi'])
    tot2pi_results.append(h_tot['2pi'])
    contamination_results.append(frac)

# talk to Hemma and Leo for plot format

# print(f'\n ------ CONTAMINATION RESULTS for Q2 = {Q2} ------')
# print(f'\n Cut | 1pi Total | 2pi Total | Contamination ')
# for i in range(len(contamination_results)):
#     print(f'{cut_results[i]:.1f} | {tot1pi_results[i]:.1f}      | {tot2pi_results[i]:.3f}     | {contamination_results[i]:.6f}')
# print('')

with open(f"figures/Q2={Q2}/contamination.txt", "w") as f:
        header = f"\n ----- CONTAMINATION RESULTS for Q2 = {Q2} -----"
        columns = f"\n{'Cut':>5} |{'1pi Total':>10} |{'2pi Total':>10} |{'Contamination':>14}"

        #print(header)
        #print(columns)
        f.write(header + "\n")
        f.write(columns + "\n")

        for cut, tot1, tot2, contam in zip(
            cut_results, tot1pi_results, tot2pi_results, contamination_results
        ):
            line = f"{cut:5.1f} |{tot1:10.1f} |{tot2:10.3f} |{contam:14.6f}"
            #print(line)
            f.write(line + "\n")

        #print()
        f.write("\n")

cut_results = np.array(cut_results)
tot1pi_results = np.array(tot1pi_results)
tot2pi_results = np.array(tot2pi_results)
contamination_results = np.array(contamination_results)

contamination_mask = (contamination_results >= 0.040) & (contamination_results <= 0.060)
plt.axvline(cut_results[contamination_mask], linestyle='--')
cth.label(r'$Q^2 = 5.0$ $\mathrm{(GeV/c^2)}$')
cth.savefig(f'figures/Q2={Q2}', f'MMr_contamination')

_running = False
print(f"Main runtime: {time.time()-t0:.2f} s\n")