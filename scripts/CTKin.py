# -*- coding: utf-8 -*-
"""
Kinematic Phase Space and Histogram script

@author: Lotus
"""

import sys
import time

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import CTHelp as cth

# CONSTANTS

m_p  = 0.9382720813     # proton mass
m_n = 0.93956563        # neutron mass
m_pi = 0.139570611      # charged pion mass

# HELPERS

def Theta_e(q2, E_beam, Eep):
    
    denom = 4.0 * E_beam * Eep
    
    # if denom <= 0.0: return None
    
    sin2half = q2 / denom
    
    # if sin2half < 1.0 or sin2half > 0.0: return None
    
    return 2.0 * np.arcsin(np.sqrt(sin2half))

def Q_squared(Ee, Eep, theta): # theta is in radians
    
    return 4*Ee*Eep*np.sin(theta/2)**2

def XBjorken(Ee, Eep, theta): # elastic (jlab)
    
    return Q_squared(Ee, Eep, theta)/(2*m_p*(Ee-Eep))
    
# MAIN

def main(input):

    tK = time.time()

    E_beam, q2, target = input

    Q2 = q2[0] + 0.5

    E_prime = np.arange(0.3, E_beam[-1], 0.01)

    # CONSTANTS

    if target == 'C':

        A = 12; Z = 6

    elif target == 'H':

        A = 1; Z = 1

    else:

        print('\nNot a valid target.\n')

    m_A = Z*m_p + (A-Z)*m_n # target mass in GeV/c^2
    m_A_1 = (A-1)*m_p       # hadronic system mass

    # RESULT ARRAYS

    results = []

    q2_results = []
    xb_results = []
    theta_e_results = []
    p_pi_results = []
    theta_pi_results = []
    k_pi_results = []
    t_results = []
    W_results = []
    Ex_results = []
    Px_results = []
    Mx_results = []
    MAx_results = []

    # INPUT ARRAYS

    Ebeam_results = []
    Eprime_results = []

    def Calc_kin(Ebeam, Eprime, q2, theta_e, forward=True):
        """
        Compute the hadronic momentum transfer t for e + p -> e' + π+ + n
        using two-body exclusive kinematics.
        forward=True gives t_min (pion emitted along q-vector).

        """
        # Energy and momentum transfer
        omega = Ebeam - Eprime
        q_abs = np.sqrt(q2 + omega**2)
        
        q_par = Ebeam - Eprime*np.cos(theta_e)
        q_perp = -1*Eprime*np.sin(theta_e)
        q_abs_comp = np.sqrt(q_par**2 + q_perp**2)

        q_vec = np.array([q_perp, 0, q_par])
        
        # Invariant mass of the hadronic system
        W2 = m_p**2 + 2*m_p*omega - q2
        if W2 <= (m_p + m_pi)**2:        
            raise ValueError("Unphysical kinematics: W too small")
        W = np.sqrt(W2)
        
        # CM pion momentum and energy
        ppi_star = np.sqrt((W2 - (m_p + m_pi)**2)*(W2 - (m_p - m_pi)**2)) / (2*W)
        Epi_star = np.sqrt(ppi_star**2 + m_pi**2)
        
        # Boost to lab frame
        beta = q_abs / (m_p + omega)
        gamma = (m_p + omega) / W
        
        # For forward pion emission (θ* = 0) → minimal |t|
        cos_theta_pq = 1.0 if forward else -1.0
        Epi_lab = gamma * (Epi_star + beta * ppi_star * cos_theta_pq)
        ppi_lab = np.sqrt(Epi_lab**2 - m_pi**2)
        
        #Pion angle and vector
        cos_q = q_par / q_abs_comp if q_abs_comp != 0 else 1.0
        theta_pi = (np.arccos(cos_q) + np.arccos(cos_theta_pq)) * 180 / np.pi
        p_pi_vec = np.array([ppi_lab * np.sin(np.deg2rad(theta_pi)), 0, ppi_lab * np.cos(np.deg2rad(theta_pi))])
        
        k_pi = np.sqrt(max(0.0, ppi_lab**2 + q_abs_comp**2 - 2.0 * ppi_lab * q_abs * cos_theta_pq))

        # Missing energy, momentum, and mass calculations of the nucleon
        E_x = omega - Epi_lab + m_p; E_A_x = omega - Epi_lab + m_A
        P_x_vec = q_vec - p_pi_vec; P_x2 = np.dot(P_x_vec, P_x_vec); P_x = np.sqrt(P_x2)
        M_x2 = E_x**2 - P_x2; M_x = np.sqrt(max(0.0, M_x2))
        M_A_x2 = E_A_x**2 - P_x2; M_A_x = np.sqrt(max(0.0, M_A_x2))
        
        # Compute t (minimal value)
        t = -q2 + m_pi**2 - 2 * (omega * Epi_lab - q_abs * ppi_lab * cos_theta_pq)
        
        return t, ppi_lab, theta_pi, k_pi, W, E_x, P_x, M_x, M_A_x

    for Ebeam in E_beam:
        #print(Ebeam)
        for Q2_val in q2:
            #print(Q2_val)
            for Eprime in E_prime:
                #print(Eprime)
                try:
                    #print(Eprime)
                    
                    theta_e = Theta_e(Q2_val, Ebeam, Eprime)
                    
                    q2_val = Q_squared(Ebeam, Eprime, theta_e)
                    
                    xb_val = XBjorken(Ebeam, Eprime, theta_e)
                    
                    theta_e *= 180/np.pi # radians conversion
                    
                    t, p_pi, theta_pi, k_pi, W, E_x, P_x, M_x, M_A_x = Calc_kin(Ebeam, Eprime, Q2_val, theta_e, forward=True)

                    t *= -1
                    
                    if np.isnan(t): continue
                    
                    if np.abs(t) > 1.0: continue
                    
                    #results
                    Ebeam_results.append(Ebeam); Eprime_results.append(Eprime)
                    q2_results.append(q2_val); xb_results.append(xb_val)
                    theta_e_results.append(theta_e); theta_pi_results.append(theta_pi)
                    t_results.append(t); W_results.append(W)
                    p_pi_results.append(p_pi); k_pi_results.append(k_pi)
                    Ex_results.append(E_x); Px_results.append(P_x)
                    Mx_results.append(M_x); MAx_results.append(M_A_x)

                    results.append([Q2_val, xb_val, t, Ebeam, theta_e, Eprime, theta_pi, p_pi, k_pi, W, E_x, P_x, M_x, M_A_x])
                    
                except Exception:
                    continue

    df = pd.DataFrame(results, columns=["Q2", "Bjorken X", "t", "Ebeam", "Theta_e", "Eprime", "Theta_p", "p_pi", "k_pi", "W", "Em", "Pm", "Mm", "MMa"])
    df.to_csv(f"figures_{target}/Q2={Q2}/t_scan_results.csv", index=False, float_format="%.6f")
    print("\n\n     Saved results to t_scan_results.csv with", len(df), "rows.\n")

    Ebeam_results = np.array(Ebeam_results)
    Eprime_results = np.array(Eprime_results)

    nu_results = Ebeam_results - Eprime_results

    q2_results = np.array(q2_results)
    xb_results = np.array(xb_results)
    theta_e_results = np.array(theta_e_results)
    p_pi_results = np.array(p_pi_results)
    theta_pi_results = np.array(theta_pi_results)
    k_pi_results = np.array(k_pi_results)
    t_results = np.array(t_results)
    W_results = np.array(W_results)
    Em_results = np.array(Ex_results)
    Pm_results = np.array(Px_results)
    Mm_results = np.array(Mx_results)
    MMa_results = np.array(MAx_results)

    # PLOTTING with constraints

    mask_xb = (xb_results > 0) & (xb_results < 0.8)
    mask_theta_e = (theta_e_results > 12) & (theta_e_results < 90)
    massk = (Mm_results >= 1e-6)

    # fixed constraints
    xb_fixed = (xb_results > 0.45) & (xb_results < 0.55)     # x ~ 0.5
    t_fixed = (t_results > 0.4) & (t_results < 0.45)         # t ~ 0.4 - 0.45

    masks = {"detection":   mask_xb & mask_theta_e,             # ranges for detection parameters
             "fix_xb":      xb_fixed & mask_xb & mask_theta_e,  # fixes xb to desired value
             "fix_t":       t_fixed & mask_xb & mask_theta_e,   # fixes t to desired value
             "mass":        massk,                              # eliminates missing mass artifacts
             "physical":    massk & mask_xb & mask_theta_e
             }
    
    # Variable dictionaries

    results = {"nu": nu_results,
            "Q2": q2_results,
            "W": W_results,
            "xb": xb_results,
            "Eprime": Eprime_results,
            "theta_e": theta_e_results,
            "t": t_results,
            "p_pi": p_pi_results,
            "theta_pi": theta_pi_results,
            "k_pi": k_pi_results,
            "Em": Em_results,
            "Pm": Pm_results,
            "Mm": Mm_results,
            "MMa": MMa_results
            }

    # KIN PLOTTING

    # (xkey, ykey, zkey, mask)
    plotPS = [
        ("xb", "Q2", "theta_e", "detection"),
        ("theta_pi", "p_pi", "Eprime", "detection"),
        ("theta_e", "Q2", "Eprime", "fix_xb"),
        ("t", "p_pi", "Eprime", "fix_xb"),
        ("t", "theta_pi", "Eprime", "fix_xb"),
        ("theta_pi", "k_pi", "Eprime", "fix_xb"),
        ("theta_e", "Q2", "t", "fix_xb"),
        ("theta_pi", "Q2", "t", "fix_xb"),
        ("W", "Q2", "t", "fix_xb"),
        ("Mm", "Pm", "Eprime", "fix_t"),
        ("Mm", "Em", "Eprime", "fix_t"),
        ("Em", "Pm", "Eprime", "fix_t"),
        ("Mm", "Q2", "Eprime", "fix_t")
    ]

    # (xkey, binsize, mask)
    plotH = [("Mm", 100, "physical"),
             ("MMa", 100, "physical")]

    # (xkey, ykey, binsize, mask)
    plot2H = [("Mm", "Pm", 100, "physical"),
            ("Mm", "Q2", 100, "physical")]
    
    def make_label(mask):

        if mask == 'fix_xb':

            label = (r'$x_b = 0.5$')

        elif mask == 'fix_t':

            label = (r'$x_b = 0.5$' '\n' r'$t = -0.4$ (GeV/c)$^2$')

        else:
            
            label = None

        if label is not None:

            cth.label(label)


    tPS = time.time()
    for xkey, ykey, zkey, mask in plotPS:

        fig, _ = plt.subplots()
        cth.scatter(results[xkey], results[ykey], results[zkey], 
                    masks[mask], fig)
        cth.format(cth.labels[xkey], cth.labels[ykey], cth.labels[zkey],
                    fr'Phase Space for $E_b=$ {Q2} GeV')
        cth.savefig(f'figures_{target}/Q2={Q2}', f'PS_{xkey}_{ykey}_{mask}')

        make_label(mask)
    print(f"\n PS plots created: {time.time() - tPS:.2f} s")

    tH = time.time()
    for key, binsize, mask in plotH:

        plt.figure()
        cth.hist(x=results[key], bins=binsize, 
                weights=None, mask=masks[mask], type='bar')
        cth.format(cth.labels[key], ylabel='Counts', colorbar=None, title=
                fr'Counts Graphs for $E_b=$ {Q2} GeV')
        cth.savefig(f'figures_{target}/Q2={Q2}', f'histo_{key}_{mask}')

        make_label(mask)
    print(f" Histograms created: {time.time() - tH:.2f} s")

    t2H = time.time()
    for xkey, ykey, binsize, mask in plot2H:

        plt.figure()
        cth.hist2D(results[xkey], results[ykey], binsize, 
                weights=None, mask=masks[mask])
        cth.format(cth.labels[xkey], cth.labels[ykey], colorbar=None, title=
                    fr'Counts Graphs for $E_b=$ {Q2} GeV')
        cth.savefig(f'figures_{target}/Q2={Q2}', f'histo2D_{xkey}_{ykey}_{mask}')

        make_label(mask)
    print(f" 2D Histograms created: {time.time() - t2H:.2f} s\n")

    print(f"\n CTKin finished: {time.time() - tK:.2f} s\n")

    plt.close('all')

    return results

if __name__ == "__main__":

    input = [np.array([10.7]),              # E_beam       
             np.arange(4.5, 5.5, 0.01),     # Q2
             'C'                            # target: C or H
    ]
    empty = False

    if len(sys.argv) < 3:
        print("\nUsage:")
        print("python CTMain.py <E_beam> <Q2> <target>\n")
        empty = True

    if not empty:
        E_beam = round(float(sys.argv[1]), 1)
        Q2 = (float(sys.argv[2]), 1)
        Q2_range = np.arange(Q2-0.5, Q2+0.5, 0.01)
        target = sys.argv[3]
        input = [E_beam, Q2_range, target]

    main(input)