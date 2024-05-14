"""
Plot the timescales as a function of radius in a disk
"""
import numpy as np
import matplotlib.pyplot as plt
import pdb
import natconst
au = natconst.au
year = natconst.year
rad = np.pi / 180

def plot1(r, t_orbit, t_stable, t_osc, t_stop):
    ax = plt.gca()

    ax.plot(r / au, t_orbit/year, label='orbit')
    ax.plot(r / au, t_stable/year, label='stable')
    ax.plot(r / au, t_osc/year, label='oscillation')
    ax.plot(r / au, t_stop/year, label='stop')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('r [au]')
    ax.set_ylabel('t [year]')

    ax.legend()

    plt.show()

def main():
    # ==== settings ====
    ms = 1 * natconst.ms # mass of star
    r = np.geomspace(0.1, 200, 20) * au

    # temperature
    q = -0.5
    T = 20 * (r / 100 / au)**(q)

    # gas surface density
    p = -1
    sigma_g = 0.1 * (r / 100 / au)**(p)

    # details of grain
    Mc = 3.0 # cork mass [g]
    Mb = 2.0 # [g]
    S = 28. # cm^2
    a = np.sqrt(S / np.pi)
    rho_s = (Mb + Mc) / a**3

    Cd = 0.44 # drag coefficient

    l_GC = 2.0 # cm

    # ==== calculate profiles ====
    cs = np.sqrt(natconst.kk * T / natconst.mp / 2.3)
    wk = np.sqrt(natconst.gg * ms / r**3)
    Hg = cs / wk
    rhog = sigma_g / np.sqrt(2*np.pi) / Hg

    # thermal velocity of the molecules
    vth = np.sqrt(8 * natconst.kk * T / np.pi / natconst.mp / 2.3)

    # gas rotation frequency
    w_g = wk * (1 + 0.5 * (Hg / r)**2 * (p + q))

    eta = - (Hg / r)**2 * (p + q)

    Q = wk * cs / np.pi / natconst.gg / sigma_g # toomre Q

    # stopping time
    t_stop = rho_s / rhog * a / vth

    # Stokes number
#    St = np.pi / 2 * a * rho_s / sigma_g
    St = t_stop * wk

    # velocity difference
    dvphi = - St**2 / (1 + St**2) * eta * wk * r

    # aerodynamic length scale
    L = 2 * (Mb + Mc) / rhog / S / Cd

    # stabilizing time
    t_stable = 2 * Mb * (1 + Mb / Mc) / rhog / S / Cd / abs(dvphi)

    # oscillation time
    t_osc = 2 * np.pi * np.sqrt(L * l_GC) / abs(dvphi)

    # orbital time
    t_orbit = 2 * np.pi / wk

    # ==== plotting ====
    plot1(r, t_orbit, t_stable, t_osc, t_stop)

    pdb.set_trace()

if __name__ == '__main__':
    main()

