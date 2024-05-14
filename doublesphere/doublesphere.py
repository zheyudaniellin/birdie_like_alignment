"""
also a double sphere model like twosphere.py 
but different in implementation. 
I utilize ratios of mass and ratios of radius

"""

import numpy as np
import matplotlib.pyplot as plt
import pdb
import natconst

def coefficients(l, rhos, epsilon, kappa, lam, rhog, vth, v):
    """
    calculate the coefficients for the equation of motion

    Parameters
    ----------
    l : ndarray
        the length between the centers of the two sphere in cm
    rhos : ndarray
        the specific weight of the grain
    epsilon : ndarray
        ratio of the sizes a1 / a2
    kappa : ndarray
        ratio of the density rho1 / rho2
    lam : ndarray
        fraction of l that the radii occupies: a1 + a2 = lam * l
    """
    # radius
    a_1 = epsilon * lam / (1 + epsilon) * l
    a_2 = lam / (1 + epsilon) * l

    # volume
    vol_1 = 4 * np.pi / 3 * a_1**3
    vol_2 = 4 * np.pi / 3 * a_2**3

    # total mass
    m = (vol_1 + vol_2) * rhos

    # individual mass
    m_1 = epsilon**3 * kappa / (1 + epsilon**3 * kappa) * m
    m_2 = 1 / (1 + epsilon**3 * kappa) * m

    # specific weight
#    rho_1 = (1 + epsilon**3) * kappa / (1 + epsilon**3 * kappa) * rhos
#    rho_2 = (1 + epsilon**3) / (1 + epsilon**3 * kappa) * rhos
    rho_1 = m_1 / vol_1
    rho_2 = m_2 / vol_2

    # lever arms
    r_1 = 1 / (1 + epsilon**3 * kappa) * l
    r_2 = epsilon**3 * kappa / (1 + epsilon**3 * kappa) * l

    # moment of inertia
    I = m_1 * (r_1**2 + 0.4 * a_1**2) + m_2 * (r_2**2 + 0.4 * a_2**2)

    D = 4*np.pi/3 * rhog * vth * (a_1**2 * r_1**2 + a_2**2 * r_2**2)
    P = 4*np.pi/3 * rhog * vth * v * (-a_1**2 * r_1 + a_2**2 * r_2)

    return {'I':I, 'D':D, 'P':P}

def timescales(I, D, P):
    """
    some timescales and frequencies
    """
    # damping time scale
    t_d = 2 * I / D

    # oscillation angular frequency
    w_o = np.sqrt(P / I)

    # oscillation time scale
    t_o = 2 * np.pi / w_o

    # damping ratio
    q = 0.5 * D / np.sqrt(I * P)

    return {'t_d':t_d, 't_o':t_o, 'w_o':w_o, 'q':q}

def drift_threshold(T, rhog, l, epsilon, kappa, lam):
    """
    the minimum level of drift to trap a grain into oscillation
    """
    fac = drift_threshold_factor(epsilon, kappa, lam)
    vthres = np.sqrt(2.3*natconst.mp * natconst.kk * T) / rhog / l**3
    return vthres

def drift_threshold_factor(epsilon, kappa, lam):
    """
    the dimensionless factor
    """
    fac = 3/16/np.sqrt(2*np.pi) * (1+epsilon)**2 * (1+epsilon**3*kappa) / lam**2 / epsilon**2 / (epsilon * kappa - 1)
    return fac

