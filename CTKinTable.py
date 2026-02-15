# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 17:17:59 2025

@author: Lotus
"""

import numpy as np
import math

# Physical constants (GeV units)
M_p  = 0.9382720813   # proton mass
M_pi = 0.139570611    # charged pion mass
PI   = math.pi

# HELPERS

def pstar_from_W(W):
    """Two-body CM momentum for p + pi final state given W (GeV)."""
    s = W*W
    term = (s - (M_p + M_pi)**2) * (s - (M_p - M_pi)**2)
    if term <= 0.0:
        return 0.0
    return 0.5/W * math.sqrt(term)

# 4-vector helpers: represent as numpy arrays [px, py, pz, E]
def vec3_mag(v):
    return math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)

def lorentz_boost_z(fourvec, beta):
    """Boost a 4-vector along +z by velocity beta. fourvec = [px,py,pz,E]."""
    if abs(beta) >= 1.0: 
        raise ValueError("beta >= 1 in boost")
    gamma = 1.0 / math.sqrt(1.0 - beta*beta)
    px, py, pz, E = fourvec
    pz_new = gamma*(pz + beta*E)
    E_new  = gamma*(E + beta*pz)
    return np.array([px, py, pz_new, E_new])

def rotate_y(fourvec, alpha):
    """Rotate spatial components of 4-vector around y-axis by angle alpha (radians)."""
    px, py, pz, E = fourvec
    ca = math.cos(alpha); sa = math.sin(alpha)
    px_r =  ca*px + sa*pz   # note rotation that maps +z into new vector
    pz_r = -sa*px + ca*pz
    return np.array([px_r, py, pz_r, E])

def compute_t_and_pionLab(Ebeam, Q2, Eprime, thetaStar):
    """
    Compute t for given Ebeam, Q2, Eprime, thetaStar (radians).
    Returns: t (GeV^2), pionLab (4-vector np.array), thetaPiLab (radians), qlab (4-vector), ppi_pt_wrt_q
    If unphysical, returns (np.nan, None, None, None, None).
    """
    nu = Ebeam - Eprime
    if nu <= 0.0:
        return np.nan, None, None, None, None

    W2 = M_p*M_p + 2.0*M_p*nu - Q2
    if W2 <= (M_p + M_pi)**2:
        return np.nan, None, None, None, None
    W = math.sqrt(W2)

    # CM pion momentum
    pstar = pstar_from_W(W)
    if pstar <= 0.0:
        return np.nan, None, None, None, None
    Epi_star = math.sqrt(pstar*pstar + M_pi*M_pi)

    # virtual photon magnitude
    qmag = math.sqrt(nu*nu + Q2)

    # hadronic CM velocity along q direction (lab frame)
    denom_beta = (M_p + nu)
    if denom_beta == 0.0:
        return np.nan, None, None, None, None
    betaCM = qmag / denom_beta
    if betaCM >= 1.0:
        return np.nan, None, None, None, None

    # pion 4-vector in CM (choose phi=0)
    px_star = pstar * math.sin(thetaStar)
    pz_star = pstar * math.cos(thetaStar)
    pionStar = np.array([px_star, 0.0, pz_star, Epi_star])

    # Boost from CM to "q-aligned lab" along +z by betaCM
    pionLab = lorentz_boost_z(pionStar, betaCM)

    # virtual photon 4-vector in q-aligned lab (aligned with +z)
    qlab = np.array([0.0, 0.0, qmag, nu])

    # Now rotate from q-aligned lab into true lab where beam is +z
    # compute electron scattering angle thetaE from Q2 and energies
    denom = 4.0 * Ebeam * Eprime
    if denom <= 0.0:
        return np.nan, None, None, None, None
    sin2half = Q2 / denom
    if sin2half < 0.0 or sin2half > 1.0:
        return np.nan, None, None, None, None
    thetaE = 2.0 * math.asin(math.sqrt(sin2half))

    # electron k' components in true lab
    kx_prime = Eprime * math.sin(thetaE)
    kz_prime = Eprime * math.cos(thetaE)
    qx_lab = - kx_prime
    qz_lab = Ebeam - kz_prime

    # rotation angle about y that maps +z onto q-hat: alpha = atan2(qx, qz)
    alpha = math.atan2(qx_lab, qz_lab)

    # rotate q and pionLab into true lab frame
    qlab = rotate_y(qlab, alpha)
    pionLab = rotate_y(pionLab, alpha)

    # compute t = (q - p_pi)^2 with metric (+, -, -, -) => (E_q - E_pi)^2 - |p_q - p_pi|^2
    E_diff = qlab[3] - pionLab[3]
    p_diff = qlab[:3] - pionLab[:3]
    t = E_diff*E_diff - np.dot(p_diff, p_diff)

    # lab pion kinematics
    p3 = pionLab[:3]
    ppi_lab = vec3_mag(p3)
    thetaPiLab = 0.0
    if ppi_lab > 0.0:
        thetaPiLab = math.acos(p3[2] / ppi_lab)  # polar angle relative to beam +z

    # compute p_perp wrt q direction
    q3 = qlab[:3]
    qmag3 = vec3_mag(q3)
    if qmag3 > 0.0:
        qhat = q3 / qmag3
        p_par_along_q = np.dot(p3, qhat)
        p_par_vec = p_par_along_q * qhat
        p_perp_vec = p3 - p_par_vec
        ppi_pt_wrt_q = vec3_mag(p_perp_vec)
    else:
        ppi_pt_wrt_q = 0.0

    return t, pionLab, thetaPiLab, qlab, ppi_pt_wrt_q

def find_thetaStar_for_t_given_Eprime(Ebeam, Q2, Eprime, t_target, Nscan=360):
    """Scan thetaStar in [0,pi] to find thetaStar that gives t_target. 
    Returns (found, thetaStar, pionLab, thetaPiLab)."""
    thmin = 1e-6
    thmax = PI - 1e-6
    prev_th = thmin
    prev_t, prev_pion, prev_thpi, _, _ = compute_t_and_pionLab(Ebeam, Q2, Eprime, prev_th)
    if not np.isfinite(prev_t):
        prev_t = 1e9
    for i in range(1, Nscan+1):
        th = thmin + (thmax - thmin) * i / float(Nscan)
        tval, pion, thpi, _, _ = compute_t_and_pionLab(Ebeam, Q2, Eprime, th)
        if not np.isfinite(tval):
            prev_th = th
            prev_t = tval
            continue
        fprev = prev_t - t_target
        fcur  = tval - t_target
        if fprev * fcur <= 0.0:   # bracket found
            a = prev_th; b = th; fa = fprev; fb = fcur
            for _ in range(60):
                c = 0.5 * (a + b)
                fc, pionc, thpic, _, _ = compute_t_and_pionLab(Ebeam, Q2, Eprime, c)
                if not np.isfinite(fc):
                    break
                fc = fc - t_target
                if abs(fc) < 1e-6:
                    return True, c, pionc, thpic
                if fa * fc <= 0.0:
                    b = c; fb = fc
                else:
                    a = c; fa = fc
            # fallback midpoint
            cmid = 0.5 * (a + b)
            tmid, pionc, thpic, _, _ = compute_t_and_pionLab(Ebeam, Q2, Eprime, cmid)
            return True, cmid, pionc, thpic
        prev_th = th
        prev_t = tval
    return False, None, None, None

# HELPERS

def main(input):
    """
    Main routine: scan E' and find forward-most thetaStar solution matching t_target.
    Prints results.
    """
    Q2, Ebeam, t_target = input

    # scanning parameters (tune as needed)
    Emin = 0.05
    Emax = Ebeam - 0.05
    dE = 0.0025   # GeV step (2.5 MeV)
    solutions = []

    Eprime = Emin
    while Eprime <= Emax + 1e-12:
        denom = 4.0 * Ebeam * Eprime
        if denom <= 0:
            Eprime += dE; continue
        sin2half = Q2 / denom
        if sin2half <= 0.0 or sin2half > 1.0:
            Eprime += dE; continue

        nu = Ebeam - Eprime
        W2 = M_p*M_p + 2.0*M_p*nu - Q2
        if W2 <= (M_p + M_pi)**2:
            Eprime += dE; continue

        ok, thetaStar_sol, pionLab_sol, thetaPiLab_sol = find_thetaStar_for_t_given_Eprime(Ebeam, Q2, Eprime, t_target)
        if ok:
            W = math.sqrt(W2)
            pstar = pstar_from_W(W)
            t_found, _, _, qlab, ppi_pt_wrt_q = compute_t_and_pionLab(Ebeam, Q2, Eprime, thetaStar_sol)
            if np.isfinite(t_found):
                solutions.append({
                    'Eprime': Eprime,
                    'thetaStar': thetaStar_sol,
                    't': t_found,
                    'pionLab': pionLab_sol,
                    'thetaPiLab': thetaPiLab_sol,
                    'W': W,
                    'pstar': pstar,
                    'qlab': qlab,
                    'ppi_pt_wrt_q': ppi_pt_wrt_q
                })
        Eprime += dE

    if len(solutions) == 0:
        print(f"No solution found for Q2={Q2}, Ebeam={Ebeam}, t={t_target}")
        print("Try increasing E' scan range, decreasing dE, or relax t tolerance.")
        return

    # choose forward-most: smallest thetaStar
    best = min(solutions, key=lambda s: s['thetaStar'])

    Eprime = best['Eprime']
    nu = Ebeam - Eprime
    W = best['W']
    pstar = best['pstar']
    pionLab = best['pionLab']
    p3 = pionLab[:3]
    ppi_lab = vec3_mag(p3)
    theta_pi_lab_deg = best['thetaPiLab'] * 180.0 / PI
    ppi_lab_pt_beam = math.sqrt(p3[0]**2 + p3[1]**2)

    # recompute q components for chosen Eprime
    denom_e = 4.0 * Ebeam * Eprime
    thetaE = 2.0 * math.asin(math.sqrt(Q2/denom_e))
    kx_p = Eprime * math.sin(thetaE)
    kz_p = Eprime * math.cos(thetaE)
    qx = - kx_p
    qz = Ebeam - kz_p
    qmag = math.sqrt(qx*qx + qz*qz)
    q3 = np.array([qx, 0.0, qz])
    qhat_mag = vec3_mag(q3)
    qhat = q3 / qhat_mag if qhat_mag > 0 else q3

    cos_pq = 0.0
    if ppi_lab > 0.0 and qhat_mag > 0.0:
        cos_pq = np.dot(p3, qhat) / ppi_lab

    kpi_formula = math.sqrt(max(0.0, ppi_lab*ppi_lab + qmag*qmag - 2.0 * ppi_lab * qmag * cos_pq))

    # t check and sqrt(-t)
    qlab = best['qlab']
    t_check = None
    if qlab is not None:
        E_diff = qlab[3] - pionLab[3]
        p_diff = qlab[:3] - pionLab[:3]
        t_check = E_diff*E_diff - np.dot(p_diff, p_diff)
        sqrt_minus_t = math.sqrt(-t_check) if t_check < 0.0 else 0.0
    else:
        t_check = float('nan')
        sqrt_minus_t = 0.0

    thetaE_deg = thetaE * 180.0 / PI

    t = best['t']
    thetaStar = best['thetaStar']*180.0/PI

    # Print
    print("=== Solution (chosen forward-most) ===")
    print(f"Input:  Q2 = {Q2:.5f} GeV^2, Ebeam = {Ebeam:.5f} GeV, t_target = {t_target:.5f} GeV^2")
    print(f"--->Computed W       = {W:.6f} GeV")
    print(f"--->Electron angle   = {thetaE_deg:.6f} deg")
    print(f"--->Scattered E'     = {Eprime:.6f} GeV")
    print(f"--->theta_pi (lab)   = {theta_pi_lab_deg:.6f} deg")
    print(f"--->p_pi (lab)       = {ppi_lab:.6f} GeV/c")
    print(f"--->k_pi (your formula) = {kpi_formula:.6f} GeV")
    print(f"nu (q0)          = {nu:.6f} GeV, |q| = {qmag:.6f} GeV")
    print(f"p*_pi (CM)       = {pstar:.6f} GeV/c")
    print(f"t (found)        = {t:.6f} GeV^2")
    print(f"theta*_CM (deg)  = {thetaStar:.6f} deg")
    print("-------------------------------------")
    print("Note: If you need more precision, reduce dE and/or increase Nscan in the theta* scanning function.\n")

    results = {"Q2": Q2, 
               "Ebeam": Ebeam, 
               "t_target": t_target, 
               "W": W, 
               "theta_e": thetaE_deg, 
               "Eprime": Eprime, 
               "theta_pi": theta_pi_lab_deg, 
               "p_pi": ppi_lab, 
               "k_pi": kpi_formula, 
               "nu": nu, 
               "q": qmag, 
               "p_star": pstar, 
               "t": t, 
               "theta_star": thetaStar
               }
    
    print("Kinematic Calculations Process Finished\n")

    return results

# Example usage:
if __name__ == "__main__":

    input = [5.0,       # Q2
             11.0,      # Ebeam
             -0.527     # t_target
             ]

    main(input)