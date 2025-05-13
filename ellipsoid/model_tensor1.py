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
def afn_E_0(abar):
    """
    analytical function for E[1]

    abar : float
    """
    if 0 < abar < 1:
        e = np.sqrt(1 - abar**2)
        res = 2 * np.arcsin(e) / e
    elif abar == 1:
        res = 2
    elif abar > 1:
        f = np.sqrt(abar**2 - 1)
        res = 2 * np.arcsinh(f) / f
    else:
        raise ValueError('abar unknown')
    return res

def afn_E_2(abar):
    """
    analytical function for E[x^2]
    """
    if 0 < abar < 1:
        e = np.sqrt(1 - abar**2)
        res = - np.sqrt(1 - e**2) / e**2 + np.arcsin(e) / e**3
    elif abar == 1:
        res = 2/3.
    elif abar > 1:
        f = np.sqrt(abar**2 - 1)
        res = np.sqrt(1 + f**2) / f**2 - np.arcsinh(f) / f**3
    else:
        raise ValueError('abar unknown')
    return res


def afn_E_4(abar):
    """
    analytical solution to E[x^4]
    """
    if 0 < abar < 1:
        e = np.sqrt(1 - abar**2)
        res = - np.sqrt(1 - e**2) * (2 * e**2 + 3) / 4 / e**4 + 3 * np.arcsin(e) / 4 / e**5
    elif abar == 1:
        res = 2/5.
    elif abar > 1:
        f = np.sqrt(abar**2 - 1)
        res = np.sqrt(1 + f**2) * (2 * f**2 - 3) / 4 / f**4 + 3 * np.arcsinh(f) / 4 / f**5
    else:
        raise ValueError('abar unknown')
    return res


def fn_H(x, abar=1):
    """
    the function for integration
    E[1 - x^2]
    """
    return (1 - x**2) / np.sqrt(1 + (abar**2 - 1) * x**2)

def get_H(abar, mode='analytical'):
    """
    use integratation
    let this function take multidimensional inputs

    E[1 - x^2]
    """
    # check if the input is a numpy array
    if type(abar) == np.ndarray:
        out = np.zeros_like(abar.flatten())

        for i, _abar in enumerate(abar):
            if mode == 'numerical':
                res = quad(fn_H, -1, 1, args=_abar)
                out[i] = res[0]
            elif mode == 'analytical':
                out[i] = afn_E_0(_abar) - afn_E_2(_abar)
            else:
                raise ValueError('mode unknown')

        out = np.reshape(out, abar.shape)
    else:
        if mode == 'numerical':
            res = quad(fn_H, -1, 1, args=abar)
            out = res[0]
        elif mode == 'analytical':
            out = afn_E_0(abar) - afn_E_2(abar)
        else:
            raise ValueError('mode unknown')

    return out

def fn_L(x, abar=1):
    """
    the function for integration for L
    """
    return (1 - x**2) * x**2 / np.sqrt(1 + (abar**2 - 1) * x**2)

def get_L(abar, mode='analytical'):
    """
    E[x^2 - x^4]
    """
    # check if the input is a numpy array
    if type(abar) == np.ndarray:
        out = np.zeros_like(abar.flatten())
        for i, _abar in enumerate(abar):
            if mode == 'numerical':
                res = quad(fn_L, -1, 1, args=_abar)
                out[i] = res[0]
            elif mode == 'analytical':
                out[i] = afn_E_2(_abar) - afn_E_4(_abar)
        out = np.reshape(out, abar.shape)
    else:
        if mode == 'numerical':
            res = quad(fn_L, -1, 1, args=abar)
            out = res[0]
        elif mode == 'analytical':
            out = afn_E_2(abar) - afn_E_4(abar)
        else:
            raise ValueError('mode unknown')

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
    D = rhog * vth * np.pi * a * c * (g**2 * H + c**2 * (1 - abar**2)**2 * L)

    # calculate K
    K = rhog * vth * v * np.pi * a * c * H * g

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

def other_timescales(rhog, vth, v, c, rhos):
    t_stop = rhos / rhog * c / vth
    t_osc_c = np.sqrt(rhos * c**2 / rhog / v / vth)

    return {'t_stop':t_stop, 't_osc_c':t_osc_c}

def main():
    """
    take a look at the values
    """
    abar = np.array([0.5, 0.9, 1.5])

    H = get_H(abar)

    pdb.set_trace()


if __name__ == '__main__':
    main()

