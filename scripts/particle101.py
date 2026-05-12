# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 17:11:44 2025

@author: Lotus
"""

import numpy as np
import matplotlib.pyplot as plt

c = 3e8 # Speed of light, scientific notation

m_pi = 0.140 # GeV/c^2
m_k = 0.493 
m_p = 0.938

l = 6 # m
t_sig = 0.2e-9 # s

#%% Plotting beta vs. momentum

beta = lambda p, m: p / np.sqrt(p**2 + m**2)
beta_sig = lambda beta: (c*beta**2 / l) * t_sig

p = np.arange(0, 5, 0.01)

plt.figure()
plt.errorbar(p, beta(p, m_pi), beta_sig(beta(p, m_pi)), label=r'$m_\pi$', color='purple')
plt.errorbar(p, beta(p, m_k), beta_sig(beta(p, m_k)), label=r'$m_k$', color='blue')
plt.errorbar(p, beta(p, m_p), beta_sig(beta(p, m_p)), label=r'$m_p$', color='red')
plt.xlabel(r'$p$'); plt.ylabel(r'$\beta$'); plt.title(r'Velocity-Momentum Relation for $\pi$, $k$, $p$')
plt.legend()
plt.show()