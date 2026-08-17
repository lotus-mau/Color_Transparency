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

def plot(input):

    Q2, target = input
    A, _ = cth.getTarget(target)

    resultsSIM, weightsSIM, _, _, _, _ = cts.calc(input)

    """
    resultsSIM: {key: 1pi or 2pi, var: Q2, t, W, etc. -> histogram}
    weightsSIM: {key: 1pi or 2pi, var: Q2, t, W, etc. -> weights of histogram}
    """

    #resultsBEAG = ctb.main(input)

    fitsBEAG = {'1piB': pd.read_csv('pionCT-beagle/fits_Q2tW_single.csv'),
                '2piB': pd.read_csv('pionCT-beagle/fits_Q2tW_multi.csv')}

    plot1D = [
        ('Q2', 100, 'dQ2'),
        ('t', 100, 'dt'),
        ('W', 100, 'dW'),
        ('mmnuc', 100, 'dMM')
    ]

    for var, binsize, label in plot1D:
        
        plt.figure()
        bin_widths = {}
        for key in resultsSIM:

            x = np.asarray(resultsSIM[key][var])
            w = np.asarray(weightsSIM['sigma'][key])

            if var == 't':
        
                mask = x > -2.5
                x = x[mask]
                w = w[mask]

            plt.hist(x, binsize, weights=w, 
                     histtype=cth.histtypes[key], color=cth.colors[key], label=cth.labels[key])

            bin_widths[key] = cth.get_binwidth(x, binsize, w)

        # for key in resultsSIM:
        #     counts, bins, patches = plt.hist(resultsSIM[key][var], binsize, weights=weightsSIM['sigma1'][key], 
        #                                      histtype='step', label=cth.labels[key], range=histrange)

        #     bin_width = (bins[-1] - bins[0]) / binsize

        #     bin_widths[key] = bin_width

        cth.format(cth.labels[var], cth.labels['sigma_i'], colorbar=None, title=
                fr'Cross-section for 1pi + multipi background')
        cth.label(fr'$Q^2={Q2}$ GeV$^2$')
        plt.legend()
        plt.ylim(0, None)
        cth.savefig(f'figures_{target}/Q2={Q2}/SIM', f'{A}{target}_{Q2}_SIM_sigma_{var}.png')

        ## DIFF CROSS SECTION PLOTTING

        plt.figure()

        weights_dsigma = {}

        for key in resultsSIM:

            x = np.asarray(resultsSIM[key][var])
            w = np.asarray(weightsSIM['sigma'][key] / bin_widths[key])      # diff cross section in fb/GeV or fb/GeV^2

            if var == 't':
        
                mask = x > -2.5
                x = x[mask]
                w = w[mask]

            plt.hist(x, binsize, weights=w, 
                     histtype=cth.histtypes[key], label=cth.labels[key], color=cth.colors[key])

        # for key in resultsSIM:

        #     weights_dsigma[key] = (weightsSIM['sigma1'][key] / bin_widths[key])       # diff cross section in fb/GeV or fb/GeV^2

        #     plt.hist(resultsSIM[key][var], binsize, weights=weights_dsigma[key], 
        #              histtype='step', label=cth.labels[key], range=histrange)

        if var == 'Q2' or var == 't' or var == 'W':

            for key in fitsBEAG:
            
                fitBEAG = fitsBEAG[key][(fitsBEAG[key]['Q2'] == Q2) & (fitsBEAG[key]['targ'] == f'{A}{target}') & (fitsBEAG[key]['var'] == var)].iloc[0]

                min, max, p0, p1, p2, p3, p4, p5, p6 = fitBEAG['min'], fitBEAG['max'], fitBEAG['p0'], fitBEAG['p1'], fitBEAG['p2'], fitBEAG['p3'], fitBEAG['p4'], fitBEAG['p5'], fitBEAG['p6']

                # print(min, max, p0, p1, p2, p3, p4, p5, p6)

                x, y = cth.poly(min, max, p0, p1, p2, p3, p4, p5, p6)
                scale = 0.1 #if key == '2piB' else 1

                plt.plot(x, scale*y(x), label=cth.labels[key], color=cth.colors[key])

        cth.format(cth.labels[var], cth.labels[label], colorbar=None, title=
                   'Differential cross-section for 1pi + multipi background')
        cth.label(fr'$Q^2={Q2}$ GeV$^2$')
        plt.legend()
        plt.ylim(0, None)
        cth.savefig(f'figures_{target}/Q2={Q2}/SIM', f'{A}{target}_{Q2}_SIM_dsigma_{label}.png') 

def total(target):

    A, Z = cth.getTarget(target)

    binsize = 100

    fitsBEAG = {'1piB': pd.read_csv('pionCT-beagle/fits_Q2tW_single.csv'),
                '2piB': pd.read_csv('pionCT-beagle/fits_Q2tW_multi.csv')}

    plot1D = [
        ('Q2', 'dQ2'),
        ('t', 'dt'),
        ('W', 'dW'),
        ('mmnuc', 'dMM')
    ]

    for var, label in plot1D:

        fig, axs = plt.subplots(2, 2, figsize=(10, 8))

        for ax, q2val in zip(axs.flat, cth.q2vals):

            resultsSIM_i, weightsSIM_i, _, _, _, _ = cts.calc([q2val, target])

            bin_widths = {}

            for key in resultsSIM_i:

                labelhist = cth.labels[key] if q2val == 5.0 else None
            
                x = np.asarray(resultsSIM_i[key][var])
                w = np.asarray(weightsSIM_i['sigma'][key])
    
                bin_widths[key] = cth.get_binwidth(x, binsize, w)

                x = np.asarray(resultsSIM_i[key][var])
                w = np.asarray(weightsSIM_i['sigma'][key] / bin_widths[key])

                if var == 't' and target == 'C':
                            
                    mask = x > -3.0
                    x = x[mask]
                    w = w[mask]

                ax.hist(x, bins=binsize, weights=w, 
                        histtype=cth.histtypes[key], color=cth.colors[key], label=labelhist)

            if var == 'Q2' or var == 't' or var == 'W':
            
                for key in fitsBEAG:

                    labelfit = cth.labels[key] if q2val == 5.0 else None
                
                    fitBEAG = fitsBEAG[key][(fitsBEAG[key]['Q2'] == q2val) & (fitsBEAG[key]['targ'] == f'{A}{target}') & (fitsBEAG[key]['var'] == var)].iloc[0]
    
                    min, max, p0, p1, p2, p3, p4, p5, p6 = fitBEAG['min'], fitBEAG['max'], fitBEAG['p0'], fitBEAG['p1'], fitBEAG['p2'], fitBEAG['p3'], fitBEAG['p4'], fitBEAG['p5'], fitBEAG['p6']
    
                    # print(min, max, p0, p1, p2, p3, p4, p5, p6)
    
                    x, y = cth.poly(min, max, p0, p1, p2, p3, p4, p5, p6)
                    scale = 0.1 if target == 'H' else 1
    
                    ax.plot(x, scale*y(x), label=labelfit, color=cth.colors[key])

            ax.set_title(fr'$Q^2$ = {q2val}')
            ax.set_xlabel(cth.labels[var])
            ax.set_ylabel(cth.labels[label])
            ax.set_ylim(0, None)

        
        fig.legend(loc="upper right")
        fig.suptitle(f'Differential Cross-sections for {A}{target}')
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        cth.savefig(f'figures_{target}', f'{A}{target}_dsigmas_{label}')

def result(nucleus):

    totals = {}

    A, Z = cth.getTarget(nucleus)

    for target in cth.targets:

        totals[target] = cts.total(target, False)

    ct = {}
    ct_errs = {}

    for key in cth.keys:

        transparency = np.asarray(totals['C']['sigma'][key]) / (Z * np.asarray(totals['H']['sigma'][key]))

        trans_err = transparency * np.sqrt(
            (np.asarray(totals['C']['sigma_err'][key]) / np.asarray(totals['C']['sigma'][key]))**2 +
            (np.asarray(totals['H']['sigma_err'][key]) / np.asarray(totals['H']['sigma'][key]))**2
        )

        ct[key] = transparency; ct_errs[key] = trans_err

    plt.figure()
    for key in cth.keys:
        plt.errorbar(totals['C']['Q2'][key], ct[key], xerr=totals['C']['Q2_err'][key], yerr=ct_errs[key],
                    fmt='o', label=cth.labels[key], color=cth.colors[key])
    cth.format(cth.labels['Q2'], cth.labels['ct'],colorbar=None,
               title=f'Transparency per setting')
    plt.legend()
    plt.savefig(f'{A}{nucleus}_transparency.png')

    plt.close('all')

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

    Q2, target = input

    plot(input)

    total(target)

    result('C')