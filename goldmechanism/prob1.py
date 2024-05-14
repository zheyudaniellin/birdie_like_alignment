"""
I'm redoing Goldreich's calculation, but utilize a formalism to calculate the polarization angle
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import dblquad
import pdb
rad = np.pi / 180

def p_function(ang, p0):
    """
    the intrinsic level of polarization from a grain 

    Parameteres
    -----------
    ang : float, ndarray
        the angle between the alignment direction and viewer
    p0 : float, ndarray
    """
    # approximation when p0 << 1
#    p = p0 * np.sin(ang)**2

    # complete
    p = p0 * np.sin(ang)**2 / (1 - p0 * np.cos(ang)**2)

    return p

def sdir(thet, phi):
    return np.array([
        np.sin(thet) * np.cos(phi), 
        np.sin(thet) * np.sin(phi), 
        np.cos(thet)
        ])

def func_Q(thet, phi, thet_s, p0):
    """
    integrate for Q
    """
    # alignment direction
    a = sdir(thet, phi)

    # sky direction
    phi_s = np.pi / 2
    n = sdir(thet_s, phi_s)

    # viewing angle
    cos_ang = np.dot(a, n)

    ang = np.arccos(cos_ang)

    # the stokes rotation
    cos_eta = (np.cos(thet) - cos_ang * np.cos(thet_s)) / (np.sin(ang) * np.sin(thet_s))

    sin_eta = np.sin(phi_s - phi) * np.sin(thet) / np.sin(ang)

    cos_2eta = 1 - 2 * sin_eta**2
    sin_2eta = 2 * sin_eta * cos_eta

    p = p_function(ang, p0)
    Q = p * cos_2eta

    return Q

def func_U(thet, phi, thet_s, p0):
    """
    integrate for U
    """
    # alignment direction
    a = sdir(thet, phi)

    # sky direction
    phi_s = np.pi / 2
    n = sdir(thet_s, phi_s)

    # viewing angle
    cos_ang = np.dot(a, n)

    ang = np.arccos(cos_ang)

    # the stokes rotation
    cos_eta = (np.cos(thet) - cos_ang * np.cos(thet_s)) / (np.sin(ang) * np.sin(thet_s))

    sin_eta = np.sin(phi_s - phi) * np.sin(thet) / np.sin(ang)

    cos_2eta = 1 - 2 * sin_eta**2
    sin_2eta = 2 * sin_eta * cos_eta

    p = p_function(ang, p0)
    U = p * sin_2eta

    return U

def solve(thet_s, p0):
    """
    solve the end polarization fraction
    """

    Q, Qres = dblquad(func_Q, 0, np.pi, 0, 2*np.pi, args=(thet_s, p0))
    U, Ures = dblquad(func_U, 0, np.pi, 0, 2*np.pi, args=(thet_s, p0))

    return Q/2./np.pi**2, U/2/np.pi**2

def main():
    # ==== settings ====
    p0 = 0.1

    # determine the direction of the sky
    # note that phi_s is assumed to be at 90 deg
    thet_s = 45 * rad

    Q, U = solve(thet_s, p0)

    print('Q = %f'%Q)
    print('U = %f'%U)

    pdb.set_trace()
if __name__ == '__main__':
    main()

