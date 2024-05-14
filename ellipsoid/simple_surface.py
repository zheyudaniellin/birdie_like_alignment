"""
just learn to integrate for the surface area

The spheroidal grain 
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import pdb

def fn(w, a, c):
    return a * c * np.sqrt(1 + ((a/c)**2 - 1) * w**2)

# ==== settings ====
a = 0.5
c = 1.0

# normalized z
nw = 11
w = np.linspace(-1, 1, nw)
ws = 0.5 * (w[1:] + w[:-1])
dw = np.diff(w)

# azimuthal angle
nphi = 15
phi = np.linspace(0, 2 * np.pi, nphi)
dphi = np.diff(phi)

#integrand = a * c * np.sqrt(1 + ((a/c)**2 - 1) * ws**2)
#S = np.sum(dphi) * np.sum(integrand * dw)
out = quad(fn, -1, 1, args=(a,c))
S = 2 * np.pi * out[0]

# expected 
e = np.sqrt(c**2 - a**2) / c
ana = 2 * np.pi * (a**2 + a * c / e * np.arcsin(e))

diff = (S - ana) / ana

print('diff = %.2e'%(diff*1e2))
pdb.set_trace()

