""" 
prob_omega.py

By taking vg=0, we can see that the torque always has an opposite sign to omega. It always resists the motion.

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import dblquad
import prob1
import pdb
from natconst import au, year, ms, gg, mp, kk
rad = np.pi / 180

def main():
    # ==== settings ====
    # prolate shape
    c = 0.1 # [cm]
    a = 0.9 * c
    rho_s = 3.0

    vol = 4.*np.pi / 3 * a**2 * c
    Iy = vol * rho_s / 5. * (a**2 + c**2)

    # location of the center of mass in body frame
    gb = - 0.1 * c

    # velocity of gas
    vg = 0.0 # [cm/s]

    # number density
    N = 1e8

    # average mass
    m = 2.3 * mp

    # temperature
    T = 20.0

    # characteristic omega from temperature
    omc = np.sqrt(2 * kk * T / Iy)

    # ==== 
    theta = 90 * rad

    omega = np.linspace(-1e-5, 1e-5, 50)

    torque = np.zeros_like(omega)

    for i in range(len(omega)):
        args = (a, c, gb, theta, omega[i], vg, N, m, T)
        res = dblquad(prob1.integrand, 0, 2*np.pi, -1, 1, args=args)
        torque[i] = res[0] * a * c


    ax = plt.gca()
    ax.plot(omega/rad, torque)
    ax.axhline(y=0, color='k')
    plt.show()

    pdb.set_trace()

if __name__ == '__main__':
    main()
