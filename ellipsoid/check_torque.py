""" 
check_torque.py

Take a look at the torque
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import dblquad
import model_integ
import model_tensor1
import pdb
from natconst import au, year, ms, gg, mp, kk
import sys
sys.path.append('../analytical')
import twosphere
rad = np.pi / 180

def main():
    # ==== settings ====
    # ellipsoid shape
    c = 0.23 # [cm]
    a = 1.2 * c
    rho_s = 3.0

    # location of the center of mass in body frame
    gb = - 0.1 * c

    # velocity of gas
    vg = 1.0 # [cm/s]

    # number density
    N = 1e10

    # average mass
    m = 2.3 * mp

    # temperature
    T = 20.0
    vth = np.sqrt(8 * kk * T / np.pi / m)

    # spin
    omega = 1e-4

    # ==== calculate the torque ====
    theta = np.linspace(0, 2*np.pi, 15)
    tq_integ = np.zeros_like(theta)

    # first the one that requires integration
    for i in range(len(theta)):
        # first the one that requires integration
        args = (a, c, gb, theta[i], omega, vg, N, m, T)
        tq_integ[i] = model_integ.do_grid(a, c, gb, theta[i], omega, vg, N, m, T, nphi=20, nw=10)

    # second, the analytical expression
    theta_tensor = np.linspace(0, 2*np.pi, 50)
    coeff = model_tensor1.coefficients(N*m, vth, vg, c, a/c, rho_s, abs(gb))
    tq_tensor = - coeff['D'] * omega - coeff['K'] * np.sin(theta_tensor)

    # ==== use the twosphere model for comparison ====
    coeff2s = twosphere.coefficients(c, 4*np.pi/3*a*c*rho_s, 0.01, 0, 0, N*m, vth, vg)

    # ==== plotting ====
    ax = plt.gca()
    ax.plot(theta/rad, tq_integ, 'ko', label='grid-based')
    ax.plot(theta_tensor/rad, tq_tensor, 'C0', label='analytical')
    ax.axhline(y=0, color='k', linestyle='--')
    ax.set_xticks(np.arange(0, 361, 45))
    ax.set_xlabel(r'$\theta$ [$^{\circ}$]')
    ax.legend()
    plt.show()

    pdb.set_trace()

if __name__ == '__main__':
    main()
