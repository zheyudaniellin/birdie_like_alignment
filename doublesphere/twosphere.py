"""
twosphere.py

function of the model
"""
import numpy as np
import matplotlib.pyplot as plt
import pdb

def coefficients(l, m, f, g, s, rhog, vth, v):
    """
    calculate the quantities

    Parameters
    ----------
    l : ndarray
        the length between the centers of the two spheres in cm
    m : ndarray
        the total mass in grams
    f : ndarray
        the different in the center of mass
    g : ndarray
        the different in size
    s : ndarray
        the fraction of empty space
    rhog : ndarray
        gas density in g/cm^3
    vth : ndarray
        gas median(?) sound speed in cm/s
    v : ndarray
        gas velocity in the frame of the grain
    """
    # mass
    m_1 = (1 + f) / 2 * m
    m_2 = (1 - f) / 2 * m

    # radius
    h = 1. - s
    a_1 = (h - g) / 2 * l
    a_2 = (h + g) / 2 * l

    # distance to the center of mass
    l_1 = (1 - f) / 2 * l
    l_2 = (1 + f) / 2 * l

    # moment of inertia
    I = m_1 * (l_1**2 + 0.4 * a_1**2) + m_2 * (l_2**2 + 0.4 * a_2**2)

    # coefficients
    D = 4*np.pi/3 * rhog * vth * (a_1**2 * l_1**2 + a_2**2 * l_2**2)
    K = 4*np.pi/3 * rhog * vth * v * (-a_1**2 *l_1 + a_2**2 * l_2)

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


