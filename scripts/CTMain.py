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

q2vals = [5.0, 6.5, 7.5, 8.5]
beam = 10.7
t = -0.527
A = 12; Z = 6; target = 'C'

int_rates = []

for Q2 in q2vals:

    tQ2 = time.time()
    print(f'\nAnalysis for Q2={Q2} commencing. \n')

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
    
    resultsKT = ctkt.main(inputKT)
    resultsK = ctk.main(inputK)
    resultsSIM = cts.main(inputSIM)

    print(f'\nAnalysis for Q2={Q2} finished: {time.time()-tQ2:.2f} s\n')

# NOTE: for integrating the count rates: 
# fit the histogram and then integrate that function

# THEN: Make a table and/or plot showing the count rate vs. Q2. 



_running = False
print(f"\nMain runtime: {time.time()-t0:.2f} s\n")