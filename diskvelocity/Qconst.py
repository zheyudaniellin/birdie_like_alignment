"""
use this to create a standard disk

I will assume a Toomre Q = 1 disk

for all the inputs, the radius should be in cm
"""
import numpy as np
import matplotlib.pyplot as plt
import natconst
from natconst import au

# shared parameters
par = {'ms':0.5 * natconst.ms, 
        'ToomreQ':1, 
        'q':0.5
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
    assume Toomre Q = 1
    """
    wk = omega(r)
    return cs(r) * wk / np.pi / natconst.gg / par['ToomreQ']

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
#    sigma = surface_density(r)
#    H = pressure_H(r)
#    rho = sigma / np.sqrt(2*np.pi) / H

    # the above should equal
    rho = 1./np.pi/np.sqrt(2*np.pi) * par['ms'] / par['ToomreQ'] / r**3
    return rho

def beta(r):
    """
    the pressure gradient
    """
    return -3 - par['q']

def H_over_R(r):
    return pressure_H(r) / r

def eta(r):
    """
    (H/R)**2 * beta
    """
    return - H_over_R(r)**2 * beta(r)


