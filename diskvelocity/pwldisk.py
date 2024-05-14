"""
use this to create a standard disk

for all the inputs, the radius should be in cm
"""
import pdb
import numpy as np
import matplotlib.pyplot as plt
import natconst
from natconst import au

# shared parameters
par = {'ms':0.5 * natconst.ms, 
       'mdisk':0.005 * natconst.ms, 
       'rin':0.1*au, 
       'rout':100*au, 
       'q':0.5, 
       'p':1.0, 
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
    R0 = 1 * au
    mdisk= par['mdisk']
    p = par['p']
    rin = par['rin']
    rout = par['rout']

    if p == 2:
        sig0 = mdisk / (2*np.pi*R0**2) / np.log(rout/rin)
    else:
        sig0 = mdisk * (2-p) / (2*np.pi*R0**p) / (rout**(2-p) - rin**(2-p))

    return sig0 * (r / R0)**(-p)

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
    cs**2 should be proportional to:    -q
    H should be proportional to:        -q/2+1.5
    rho should be proportional to:      -p + q/2 - 1.5

    pressure should be proportional to: -p - q/2 - 1.5
    """
    return -par['p'] - par['q']/2 - 1.5 + np.zeros_like(r)

def H_over_R(r):
    return pressure_H(r) / r

def eta(r):
    """
    (H/R)**2 * beta
    """
    return - H_over_R(r)**2 * beta(r)


