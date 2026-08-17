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
import CTBeagle as ctb

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

q2vals = [5.0]
beam = 10.7
t = -0.527
A = 12; Z = 6; target = 'C'

binsize = 100

results = {}        # : {Q2, what was run}

for Q2 in q2vals:

    tQ2 = time.time()
    print(f'\nAnalysis for Q2={Q2} commencing. \n')

    results[Q2] = {}

    inputKT = [
        Q2,       
        beam,           # Ebeam
        t               # t_target
    ]

    inputK = [
        np.array([beam]),                       # E_beam       
        np.arange(Q2-0.5, Q2+0.5, 0.01),        # Q2 range
        target                                  # target: C or H
    ]
    
    inputSIM = [
        Q2,
        target          # target: C or H
    ]

    inputBEAG = [
        Q2,
        target          # target: C or H
    ]
    
    results[Q2]['KT'] = ctkt.main(inputKT)
    results[Q2]['K'] = ctk.main(inputK)
    results[Q2]['SIM'], results[Q2]['Cont'] = cts.main(inputSIM)
    results[Q2]['BEAG'] = ctb.main(inputBEAG)

    _, weightsSIM, _, _, _, _ = cts.calc(inputSIM)

    plt.figure()
    cth.hist(results[Q2]['SIM']['1pi']['mmnuc'], binsize, 
                 weights=weightsSIM['rates']['1pi']*3500, mask=None, type='step')
    cth.hist(results[Q2]['SIM']['2pi']['mmnuc'], binsize, 
                 weights=weightsSIM['rates']['2pi']*30000, mask=None, type='step')
    cth.hist(results[Q2]['BEAG']['miss_mass'], binsize, weights=None, mask=None, type='step') 
    cth.format(cth.labels['mmnuc'], f'Counts', colorbar=None, title=
               f'Comparison of SIMC and BeAGLE, Q2={Q2}')
    cth.savefig(f'figures_{target}/Q2={Q2}', f'SIMCBEAG_mmnuc.png')

    print(f'\nAnalysis for Q2={Q2} finished: {time.time()-tQ2:.2f} s\n')



# NOTE: for integrating the count rates: 
# fit the histogram and then integrate that function

# THEN: Make a table and/or plot showing the count rate vs. Q2. -> Done

_running = False
print(f"\nMain runtime: {time.time()-t0:.2f} s\n")