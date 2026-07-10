# -*- coding: utf-8 -*-
"""
Beagle Code: used to process Beagle output csv files 
and analyze output, later compare to SIMC output in Main.

@author: Lotus
"""

import sys
import time

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import CTHelp as cth

def load(file_csv):

    df = pd.read_csv(file_csv)

    """
    Read csv columns and store in separate variables

    Variables:

    x                   # Bjorken X
    Q2                  # Momentum Transfer
    W                   # invariant mass
    miss_mass           # Missing Mass
    phih                # Azimuthal angle of Hadron
    pt2                 # transverse momentum squared
    t                   # hadronic momentum transfer
    secondary_pips      # Secondary Pi+
    secondary_pims      # Secondary Pi-
    secondary_pizs      # Secondary Pi0
    secondary_esum      # Secondary energy sum 
    accepted            # However many particles were accepted

    """

    return {col: df[col].to_numpy() for col in df.columns}

def calc(input):

    Q2, target = input

    A, _ = cth.getTarget(target)

    results = {}

    for i in range(1, 11):
        data = load(f'pionCT-beagle/{A}{target}_{Q2}_{i}.csv')

        for key, value in data.items():
            results.setdefault(key, []).append(value)

    return {key: np.concatenate(arrays) for key, arrays in results.items()}

def main(input):

    tBEAG = time.time()

    Q2, target = input

    results = calc(input)

    for key in results:
        plt.figure()
        cth.hist(results[key], bins=100, weights=None, mask=None, type='step')
        cth.format(f'{key}', f'Counts', colorbar=None,
               title=f'Beagle Plots for 1pi, Q2 = {Q2}')
        cth.savefig(f'figures_{target}/Q2={Q2}/BeAGLE', f'BEAG_{key}.png')

    plt.close('all')

    print(f'\nBeagle Analysis Done: {time.time() - tBEAG:.2f} s\n')
    
    return results

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

    # code graphing the fitted polynomials