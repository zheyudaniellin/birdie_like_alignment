"""
model_tensor1.py

Calculate the results based on the tensor formalism
Assume that the center of mass is offset along the axis of symmetry
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quadrature as quad
import pdb
rad = np.pi / 180

# =========================================================
# integration constants 
# =========================================================
def fn_H(x, abar=1):
    """
    the function for integration
    """
    return (1 - x**2) / np.sqrt(1 + (abar**2 - 1) * x**2)

def get_H(abar):
    """
    use integratation
    let this function take multidimensional inputs
    """
    # check if the input is a numpy array
    if type(abar) == np.ndarray:
        out = np.zeros_like(abar.flatten())
        for i, _abar in enumerate(abar):
            res = quad(fn_H, -1, 1, args=_abar)
            out[i] = res[0]
        out = np.reshape(out, abar.shape)
    else:
        res = quad(fn_H, -1, 1, args=abar)
        out = res[0]

    return out

def fn_L(x, abar=1):
    """
    the function for integration for L
    """
    return (1 - x**2) * x**2 / np.sqrt(1 + (abar**2 - 1) * x**2)

def get_L(abar):
    """
    """
    # check if the input is a numpy array
    if type(abar) == np.ndarray:
        out = np.zeros_like(abar.flatten())
        for i, _abar in enumerate(abar):
            res = quad(fn_L, -1, 1, args=_abar)
            out[i] = res[0]
        out = np.reshape(out, abar.shape)
    else:
        res = quad(fn_L, -1, 1, args=abar)
        out = res[0]

    return out

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

def coefficients(rhog, vth, v, c, abar, rhos, gbar):
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
    gbar : ndarray
        The offset of the geometric center from the center of mass g, but normalized by c
    """
    # H constant
    H = get_H(abar)

    # L constant
    L = get_L(abar)

    # ellipsoid
    a = abar * c

    g = gbar * c

    # calculate I
    I = get_I(c, abar, rhos, g)

    # calculate D
    D = 2 * rhog * vth * np.pi * a * c * (g**2 * H + c**2 * (1 - abar**2)**2 * L)

    # calculate K
    K = 2 * rhog * vth * v * np.pi * a * c * H * g

    return {'I':I, 'D':D, 'K':K}

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
    abar = np.array([0.5, 0.9, 1.5])

    H = get_H(abar)

    pdb.set_trace()


if __name__ == '__main__':
    main()

