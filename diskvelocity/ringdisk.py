"""
use a disk with rings
for all the inputs, the radius should be in cm
"""
import pdb
import numpy as np
import matplotlib.pyplot as plt
import natconst
from natconst import au

# shared parameters
par = {'ms':0.5 * natconst.ms, 
       'r':[40*au, 70*au], 
       'w':[5*au, 5*au],
       'sig':[10, 10], 
       'q':0.5, 
        }

# disk radial profiles
def temperature(r):
    return 200 * (r / au)**(-par['q'])

def cs(r):
    cs = np.sqrt(natconst.kk * temperature(r) / 2.3 / natconst.mp)
    return cs

def vth(r):
    vth = np.sqrt(8 * natconst.kk * temperature(r) / np.pi / 2.3 / natconst.mp)
    return vth

def surface_density(r):
    """
    takes the power-law index of the surface density
    sigma = sig0 * (R / R0)**(-p)
    """
    cen = par['r']
    width = par['w']
    amp = par['sig']

    sig = np.zeros_like(r)
    for i in range(len(cen)):
        sig += amp[i] * np.exp(-0.5 * ((r - cen[i]) / width[i])**2)

    return sig

def omega(r):
    return np.sqrt(natconst.gg * par['ms'] / r**3)

def pressure_H(r):
    """
    pressure scale height
    """
    return cs(r) / omega(r)

def midplane_density(r):
    """
    midplane gas density
    """
    sigma = surface_density(r)
    H = pressure_H(r)
    rho = sigma / np.sqrt(2*np.pi) / H

    return rho

def beta(r):
    """
    the pressure gradient
    d lnP / d lnR
    """
    sigma = surface_density(r)
    P = midplane_density(r) * cs(r)**2
    lnP = np.log(P)
    lnr = np.log(r)
    beta_cell = (lnP[1:] - lnP[:-1]) / (lnr[1:] - lnr[:-1])
    beta = np.zeros_like(r)
    beta[1:-1] = 0.5 * (beta_cell[1:] + beta_cell[:-1])
    beta[0] = beta_cell[0]
    beta[-1] = beta_cell[-1]
    return beta

def H_over_R(r):
    return pressure_H(r) / r

def eta(r):
    """
    (H/R)**2 * beta
    """
    return - H_over_R(r)**2 * beta(r)


