"""
model_tensor2.py

Calculate the results based on the tensor formalism

The center of mass offset can be within the b1-b3 plane of the ellipsoid (perpendicular to the axis of oscillation)
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quadrature as quad
import pdb
rad = np.pi / 180

# =========================================================
# integration constants 
# =========================================================
def fn_W(x, abar=1):
    """
    """
    return x**2 / np.sqrt(1 + (abar**2 - 1) * x**2)

def get_W(abar):
    """
    """
    res = quad(fn_W, -1, 1, args=abar)

    return res[0]

def fn_H(x, abar=1):
    """
    the function for integration
    """
    return (1 - x**2) / np.sqrt(1 + (abar**2 - 1) * x**2)

def get_H(abar):
    """
    use integratation
    """
    res = quad(fn_H, -1, 1, args=abar)

    return res[0]

def fn_L(x, abar=1):
    """
    the function for integration for L
    """
    return (1 - x**2) * x**2 / np.sqrt(1 + (abar**2 - 1) * x**2)

def get_L(abar):
    """
    """
    res = quad(fn_L, -1, 1, args=abar)

    return res[0]

# =========================================================
def get_I(c, abar, rhos, g):
    """
    calculate the moment of inertia
    """
    a = abar * c
    I_perp = 4 * np.pi / 15. * a**2 * c * (a**2 + c**2) * rhos

    # total mass
    M = 4 * np.pi / 3 * a**2 * c * rhos

    return I_perp + g**2 * M

def coefficients(rhog, vth, v, c, abar, rhos, g, psi=0):
    """
    Wrapper function to calculate all the necessary coefficients

    Most of the parameters can take a float or an array. If you want to give two parameters with an array, the array has to be of the same size. 

    Parameters
    ----------
    rhog : float or ndarray
    vth : float or ndarray
    v : float or ndarray
        The drift velocity of gas relative to the dust
    c : ndarray
        The half length of the axis of symmetry. 
    abar : float
        This can only be a number because of we need to calculate some integration constants which cannot take arrays. 
    g : ndarray
    psi : ndarray
        The location from the center of mass forms an angle with the b_3 basis
    """
    # H constant
    H = get_H(abar)

    # L constant
    L = get_L(abar)

    # W constant
    W = get_W(abar)

    # ellipsoid
    a = abar * c

    # calculate I
    I = get_I(c, abar, rhos, g)

    # calculate D
    D = 2 * rhog * vth * np.pi * a * c * (
            g**2 * H * np.cos(psi)**2 
            + g**2 * (2 * abar**2 * W) * np.sin(psi)**2 
            + c**2 * (1 - abar**2)**2 * L
            )

    # calculate Ks
    Ks = 2 * rhog * vth * v * np.pi * a * c * H * g * np.cos(psi)

    # calculate Kc
    Kc = 2 * rhog * vth * v * np.pi * a * c * (2*abar**2*W) * g * np.sin(psi)

    return {'I':I, 'D':D, 'Ks':Ks, 'Kc':Kc}

def timescales(I, D, K):
    """
    some timescales and frequencies
    """
    # damping time scale
    t_d = 2 * I / D

    # oscillation angular frequency
    w_o = np.sqrt(K / I)

    # oscillation time scale
    t_o = 2 * np.pi / w_o

    return {'t_d':t_d, 't_o':t_o, 'w_o':w_o}

def main():
    """
    take a look at the values
    """
    c = 0.1 # [cm]
    abar = 0.9

if __name__ == '__main__':
    main()

