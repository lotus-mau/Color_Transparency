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
    miss_mass           # Missing Mass
    phih                # Azimuthal angle of Hadron
    pt2                 # transverse momentum squared
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

    Q2, target = input

    results = calc(input)

    plt.figure()
    cth.hist(results['miss_mass'], bins=100, weights=None, mask=None, type='step')
    cth.format(cth.labels['missmass'], f'Counts', colorbar=None,
               title=f'Beagle Plots for 1pi, Q2 = {Q2}')
    cth.savefig(f'figures_{target}/Q2={Q2}', f'BEAG_missmass')

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