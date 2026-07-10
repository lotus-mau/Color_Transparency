# -*- coding: utf-8 -*-
"""
Simulation Code: for loading and processing SIMC files,
then analyze the histograms -> contamination from dpi background

@author: Lotus
"""

import sys
import time
import uproot

import numpy as np
import matplotlib.pyplot as plt

import CTHelp as cth

def load(file_root):
    # Open the ROOT file
    file = uproot.open(file_root)

    # Get the tree
    tree = file["h10"]

    """
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
    """

    return {v: tree[v].array() for v in cth.vars}

def calc(input):

    Q2, target = input
    
    A, Z = cth.getTarget(target)

    specs = cth.get_specs()

    results = {
            "1pi": load(f"pionCT-simc/pion_{cth.q2val_tags[Q2]}_{target}.root"), 
            "2pi": load(f"pionCT-simc/pion_{cth.q2val_tags[Q2]}_{target}_multipi.root")
        }

    # coefficients
    weights_og = {
        key: results[key]['Weight'] * specs[cth.q2val_tags[Q2]][key][target]['normfac'] / cth.ngen 
        for key in results
    }

    weights_lum = {
        key: results[key]['Weight'] * cth.luminosity(Q2,target)[key] / (specs[cth.q2val_tags[Q2]][key][target]['Ntried'] * cth.ngen)
        for key in results      
    }

    # print(cth.luminosity(Q2,target)['1pi'])
    # print(specs[cth.q2val_tags[Q2]]['1pi'][target]['luminosity'])

    weights_rate = {
        key: weights_og[key] * cth.current         
        for key in results      # [Counts/mC] * [mC/s] -> [Counts/s]
    }

    weights_sigma = {
        key: weights_rate[key] / cth.luminosity(Q2,target)[key]
        for key in results      # [Counts/s] / [s^-1*fb^-1] -> Counts [fb]
    }

    weights_sigma1 = {          # equivalent to above :D!
        key: weights_og[key] / specs[cth.q2val_tags[Q2]][key][target]['luminosity']
        for key in results      # [Counts/mC] / [mC^-1*fb^-1]
    }

    # bin_widths = {}
    # for key in results:
    #     counts, bins, patches = plt.hist(results[key]['Q2'], 100, weights=weights_sigma[key], 
    #                                             histtype='step', label=cth.labels[key])

    #     bin_width = (bins[-1] - bins[0]) / 100

    #     bin_widths[key] = bin_width

    # weights_dsigma = {
    #     key: weights_sigma[key] / bin_widths[key]
    #     for key in results
    # }

    if target == 'C':

        weights_rate['2pi'] = weights_rate['2pi'] * 2

    Ehad = {key: np.sqrt(results[key]['phad']**2 + cth.m_pi**2) for key in results}
    ppi = {5.0: 5.111, 6.5: 6.715, 7.5: 7.784, 8.5: 8.430}

    Ex = {key: results[key]['nu'] + specs[cth.q2val_tags[Q2]][key][target]['mass'] - Ehad[key] for key in results}
    reconmass1 = {key: np.sqrt(Ex[key]**2 - results[key]['Pm']**2) for key in results}              # equivalent to mmnuc
    reconmass2 = {key: np.sqrt(results[key]['Em']**2 - results[key]['Pm']**2) for key in results}   # equivalent to mmnuc

    weights = {'og': weights_og, 'lum': weights_lum, 'rate': weights_rate, 
               'sigma': weights_sigma, 'sigma1': weights_sigma1#, 'dsigma': weights_dsigma
               }
    
    return results, weights, Ehad, Ex, reconmass1, reconmass2

def main(input):
    #print("Alan is goated")
    Q2, target = input

    A, Z = cth.getTarget(target)

    results, weights, Ehad, Ex, reconmass1, reconmass2 = calc(input)

    binsize = 100

    tSIM = time.time()

    """ plt.figure() referenced below showed 1pi, 2pi, and norad -> not needed

    above plot is... interesting. ask Holly
    answer: norad was just to compare radiative effects on the reaction, 
    plus, plot makes sense relatively. """

    # PLOTTING MM spectrum of 1pi + multipi

    plt.figure() 
    for key in reconmass1:

        cth.hist(reconmass1[key], binsize, weights['rate'][key], mask=None, type='step')
        if target == 'H':
            plt.yscale('log')
        
    cth.format(cth.labels['mmnuc'], cth.labels['Counts_s'], colorbar=None, title=
            fr'Missing Mass Rates Graphs for 1pi + multipi background')

    """ plot:

    Mm spectrum for og (1pi) and 2pi, 
    so that 2pi tail corresponds to 1pi, 
    then no rad histo, which should be thinner # not really for ct, more for understanding
    then rho histo. which should be wider and lower
    All on one figure to see. 

    then, cut Mm below a certain threshold.

    target diff:

    density_diff = abs(cth.density['1pi'] - cth.density['norad']) / cth.density['1pi']
    length_diff = abs(cth.length['1pi'] - cth.length['norad']) / cth.length['1pi']

    print(f"\n Target Difference: \n Density: {density_diff:.3f} \n Length {length_diff:.3f}\n")
    print(f"1pi to norad Peak normalization: {weight_norm:.3f}\n")

    table:
    place cuts around left tail of 2pi, 
    and calculate the contamination of 2pi related to 1pi. (fraction of 2pi from 1pi)
    cut on MM | 1pi rate | 2pi rate | signal rate (subtracted) """

    # INTEGRATING HISTOGRAMS

    cut_slider = np.arange(np.min(results['1pi']['mmnuc']), np.max(results['1pi']['mmnuc']), 1/binsize) # integration range, plus steps

    cut_results = []
    cut1pi_results = []
    cut2pi_results = []
    contamination_results = []

    for cut in cut_slider: # Missing Mass cuts

        h_tot = {}

        for key in results:

            mask = results[key]['mmnuc'] < cut

            hist, edges = np.histogram(
                results[key]['mmnuc'][mask],
                bins=binsize,
                weights=weights['rate'][key][mask]
            )

            h_tot[key] = np.sum(hist)
            
        frac = h_tot['2pi'] / h_tot['1pi']

        cut_results.append(cut)
        cut1pi_results.append(h_tot['1pi'])
        cut2pi_results.append(h_tot['2pi'])
        contamination_results.append(frac)

    # talk to Hemma and Leo for plot format

    with open(f"figures_{target}/Q2={Q2}/contamination.txt", "w") as f:
            header = f"\n ----- CONTAMINATION RESULTS for Q2 = {Q2} -----"
            columns = f"\n{'Cut':>5} |{'1pi Total':>10} |{'2pi Total':>10} |{'Contamination':>14}"

            #print(header)
            #print(columns)
            f.write(header + "\n")
            f.write(columns + "\n")

            for cut, tot1, tot2, contam in zip(
                cut_results, cut1pi_results, cut2pi_results, contamination_results
            ):
                line = f"{cut:5.1f} |{tot1:10.3f} |{tot2:10.3f} |{contam:14.6f}"
                #print(line)
                f.write(line + "\n")

            #print()
            f.write("\n")

    cut_results = np.array(cut_results)
    cut1pi_results = np.array(cut1pi_results)
    cut2pi_results = np.array(cut2pi_results)
    contamination_results = np.array(contamination_results)

    #contamination_mask = (contamination_results >= 0.040) & (contamination_results <= 0.060)
    #plt.axvline(cut_results[contamination_mask], linestyle='--')
    cth.label(fr'$Q^2 = {Q2}$ $GeV/c^2$')
    cth.savefig(f'figures_{target}/Q2={Q2}/SIM', f'{A}{target}_{Q2}_SIM_contamination.png')

    plot1D = [
        ('epsilon', 100),
        ('Pm', 100),
        ('Em', 100),
        ('mmnuc', 100),
        ('Mhadron', 100),
        ('missmass', 100),
        ('phad', 100),
        ('nu', 100),
        ('t', 100),
        ('dpimm', 100),
        ('W', 100)
    ]

    # (xkey, ykey, binsize, weight)
    plot2D = [
        ("mmnuc", "Pm", 100, weights['rate']["1pi"]),
        ("mmnuc", "Q2", 100, weights['rate']["1pi"])
    ]

    for xkey, ykey, binsize, weight in plot2D:

        x = np.asarray(results['1pi'][xkey])
        y = np.asarray(results['1pi'][ykey])
        w = np.asarray(weight)

        plt.figure()
        cth.hist2D(x, y, binsize, weights=w, mask=None)
        cth.format(cth.labels[xkey], cth.labels[ykey], colorbar='Counts/s', 
                title=fr'Counts Graphs for $E_b=$ {Q2} GeV')
        cth.label(fr'$Q^2={Q2}$ GeV$^2$')
        cth.savefig(f'figures_{target}/Q2={Q2}/SIM', f'{A}{target}_{Q2}_SIM_{xkey}_{ykey}.png')

    for var, binsize in plot1D:
        plt.figure()
        for key in results:
            cth.hist(results[key][var], binsize, weights=weights['rate'][key], mask=None, type='step')
        cth.format(cth.labels[var], cth.labels['Counts_s'], colorbar=None, title=
                fr'Graphs for 1pi + multipi background')
        cth.label(fr'$Q^2={Q2}$ GeV$^2$')
        cth.savefig(f'figures_{target}/Q2={Q2}/SIM', f'{A}{target}_{Q2}_SIM_{var}.png')

    plt.figure()
    for key in Ex:
        cth.hist(Ex[key], 100, weights['rate'][key], mask=None, type='step')
    cth.format(cth.labels['Em'], cth.labels['Counts_s'], colorbar=None, title=
               fr'Graphs for 1pi + multipi background')
    cth.label(fr'$Q^2={Q2}$ GeV$^2$')
    cth.savefig(f'figures_{target}/Q2={Q2}/SIM', f'{A}{target}_{Q2}_SIM_Ex.png')

    plt.close('all')

    print(f'\n CTSim finished: {time.time()-tSIM:.2f} s\n')

    return results, {'cut': cut_results, '1pi': cut1pi_results, '2pi': cut2pi_results, 'contamination': contamination_results}

def test(target):

    A, _ = cth.getTarget(target)

    q2vals = [5.0, 6.5, 7.5, 8.5]

    keys = ['1pi', '2pi']

    plt.figure()

    sigmas = {key: {} for key in keys}
    sigma_errs = {key: {} for key in keys}
    Q2_aves = {key: {} for key in keys}
    Q2_errs = {key: {} for key in keys}

    for q2val in q2vals:

        input = [q2val, 'C']

        results, weights, _, _, _, _ = calc(input)

        for key in results:
            
            hist_sigma, _ = np.histogram(results[key]['Q2'], 100, weights=weights['sigma'][key])

            hist_err2, _ = np.histogram(results[key]['Q2'], 100, weights=weights['sigma'][key]**2)
            hist_err = np.sqrt(hist_err2)

            sigma = np.sum(hist_sigma) ; sigmas[key][q2val] = sigma
            sigma_err = np.sqrt(np.sum(hist_err**2)) ; sigma_errs[key][q2val] = sigma_err

            Q2_ave = np.average(results[key]['Q2'], weights=weights['sigma'][key])
            Q2_aves[key][q2val] = Q2_ave
            Q2_err = np.sqrt(np.average((results[key]['Q2'] - Q2_ave)**2, weights=weights['sigma'][key]))
            Q2_errs[key][q2val] = Q2_err

    sigmas = {key: [sigmas[key][q2val] for q2val in q2vals] for key in keys}
    sigma_errs = {key: [sigma_errs[key][q2val] for q2val in q2vals] for key in keys}
    Q2_aves = {key: [Q2_aves[key][q2val] for q2val in q2vals] for key in keys}
    Q2_errs = {key: [Q2_errs[key][q2val] for q2val in q2vals] for key in keys}

    for key in keys:
        plt.errorbar(Q2_aves[key], sigmas[key], xerr=Q2_errs[key], yerr=sigma_errs[key], 
                    fmt='o', label=cth.labels[key])
    
    cth.format(cth.labels['Q2'], cth.labels['sigma'],colorbar=None,
               title=f'Total Cross Section per setting, for 1pi + multipi')
    plt.legend()
    cth.savefig(f'figures_{target}', f'{A}{target}_totalsigmas.png')



if __name__ == "__main__":

    target = 'C'
    Q2 = 5.0

    input = [Q2, target]
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

    test(target)