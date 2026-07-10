"""
Used for comparison
"""

import sys
import time

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import CTBeagle as ctb
import CTSim as cts
import CTHelp as cth

def main(input):

    Q2, target = input
    A, _ = cth.getTarget(target)

    resultsSIM, weightsSIM, _, _, _, _ = cts.calc(input)

    """
    resultsSIM: {key: 1pi or 2pi, var: Q2, t, W, etc. -> histogram}
    weightsSIM: {key: 1pi or 2pi, var: Q2, t, W, etc. -> weights of histogram}
    """

    resultsBEAG = ctb.main(input)

    fitsBEAG = pd.read_csv('pionCT-beagle/fits_Q2tW.csv')

    plot1D = [
        ('Q2', 100, 'dQ2'),
        ('t', 100, 'dt'),
        ('W', 100, 'dW'),
        ('mmnuc', 100, 'dMM')
    ]

    for key in resultsSIM:
        resultsSIM[key]['t'] = -1*resultsSIM[key]['t']     # for coinciding BeAGLE and SIMC results to show a negative t.

    for var, binsize, label in plot1D:

        histrange = (-3, 0) if var == 't' else None     # restricting displayed range

        plt.figure()
        bin_widths = {}
        for key in resultsSIM:

            counts, bins, patches = plt.hist(resultsSIM[key][var], binsize, weights=weightsSIM['sigma'][key], 
                                             histtype='step', range=histrange)

            bin_width = (bins[-1] - bins[0]) / binsize

            bin_widths[key] = bin_width

        for key in resultsSIM:
            counts, bins, patches = plt.hist(resultsSIM[key][var], binsize, weights=weightsSIM['sigma1'][key], 
                                             histtype='step', label=cth.labels[key], range=histrange)

            bin_width = (bins[-1] - bins[0]) / binsize

            bin_widths[key] = bin_width

        cth.format(cth.labels[var], cth.labels['sigma'], colorbar=None, title=
                fr'Cross-section for 1pi + multipi background')
        cth.label(fr'$Q^2={Q2}$ GeV$^2$')
        plt.legend()
        cth.savefig(f'figures_{target}/Q2={Q2}/SIM', f'{A}{target}_{Q2}_SIM_sigma_{var}.png')

        plt.figure()

        weights_dsigma = {}

        for key in resultsSIM:

            weights_dsigma[key] = (weightsSIM['sigma'][key] / bin_widths[key])       # diff cross section in fb/GeV or fb/GeV^2

            plt.hist(resultsSIM[key][var], binsize, weights=weights_dsigma[key], 
                     histtype='step', range=histrange)

        for key in resultsSIM:

            weights_dsigma[key] = (weightsSIM['sigma1'][key] / bin_widths[key])       # diff cross section in fb/GeV or fb/GeV^2

            plt.hist(resultsSIM[key][var], binsize, weights=weights_dsigma[key], 
                     histtype='step', label=cth.labels[key], range=histrange)

        # if var == 'Q2' or var == 't' or var == 'W':

        #     fitBEAG = fitsBEAG[(fitsBEAG['Q2'] == Q2) & (fitsBEAG['targ'] == f'{A}{target}') & (fitsBEAG['var'] == var)].iloc[0]

        #     min, max, p0, p1, p2, p3, p4, p5, p6 = fitBEAG['min'], fitBEAG['max'], fitBEAG['p0'], fitBEAG['p1'], fitBEAG['p2'], fitBEAG['p3'], fitBEAG['p4'], fitBEAG['p5'], fitBEAG['p6']

        #     scale = 10e-1

        #     cth.poly(min, max, p0, p1, p2, p3, p4, p5, p6, scale)

        cth.format(cth.labels[var], cth.labels[label], colorbar=None, title=
                   'Differential cross-section for 1pi + multipi background')
        cth.label(fr'$Q^2={Q2}$ GeV$^2$')
        plt.legend()
        cth.savefig(f'figures_{target}/Q2={Q2}/SIM', f'{A}{target}_{Q2}_SIM_dsigma_{label}.png') 

if __name__ == "__main__":

    input = [5.0, 'C']
    empty = False

    if len(sys.argv) < 3:
        print("\nUsage:")
        print("python CTMain.py <Q2> <target>\n")
        empty = True

    if not empty:
        Q2 = round(float(sys.argv[1]), 1)
        target = sys.argv[2]
        input = [Q2, target]

    main(input)