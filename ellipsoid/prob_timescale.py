"""
prob_timescale.py

Take a look at the alignment timescale in a disk as a function of radius
"""
import numpy as np
import matplotlib.pyplot as plt
from natconst import au, year
import natconst
import model_tensor1
import diskvel
import pdb
import sys
sys.path.append('../diskvelocity')
import standard_disk

def main():
    # ==== settings ====
    ms = 0.5 * natconst.ms
    r = np.geomspace(0.1, 200, 20) * au
    alpha = 1e-3

    # temperature
    q = 0.5
    T = 20 * (r / 100 / au)**(-q)
    # sound speed
    cs = np.sqrt(natconst.kk * T / natconst.mp / 2.3)

    # keplerian frequency
    wk = np.sqrt(natconst.gg * ms / r**3)

    # gas density
    Q = 1
    sigma_g = cs * wk / np.pi / G / Q
#    p = 1.
#    sigma_g = 10. * (r / 100 / au)**(-p)

    # details of grain
    c = 0.1 # [mm]
    abar = 0.9
    rho_s = 3.0
    gbar = 0.01 # g / c

    # ==== some profiles ====
    vk = np.sqrt(natconst.gg * ms / r)
    t_kep = 2 * np.pi / wk

    # pressure scale height
    h = cs / wk

    # toomre Q parameter
    Q = cs * wk / np.pi / natconst.gg / sigma_g

    # mean thermal velocity 
    vth = np.sqrt(8 * natconst.kk * T / np.pi / natconst.mp / 2.3)

    # midplane density
    rho_g = sigma_g / h / np.sqrt(2*np.pi)
    rho_pwl = 1.5 - q/2 - p

    # stopping time
    t_stop = rho_s / rho_g * c / vth

    # Stokes number
    St = t_stop * wk 

    # relative velocity
    dvphi = diskvel.rel_v_phi(St, h/r, vk, rho_pwl, q)
    dvr = diskvel.rel_v_r(St, alpha, h/r, vk, rho_pwl, q)
    dv = np.sqrt(dvr**2 + dvphi**2)

    # calculate oscillation of the grain
    coef = model_tensor1.coefficients(rho_g, vth, dv, c, abar, rho_s, gbar)
    time = model_tensor1.timescales(coef['I'], coef['D'], coef['K'])
    t_damp = time['t_d']
    t_osc = time['t_o']

    # ==== plotting ====
    fig, axgrid = plt.subplots(2, 2, sharex=False, sharey=False, 
        squeeze=False,
        figsize=(10,7))
    axes = axgrid.flatten()

    ax = axes[0]
    ax.plot(r/au, t_damp/year, label='damping')
    ax.plot(r/au, t_osc/year, label='oscillation')
    ax.plot(r/au, t_kep/year, 'k--', label='orbital')
    ax.plot(r/au, t_stop/year, 'k:', label='stopping')
    ax.set_yscale('log')
    ax.legend()
    ax.set_title('time [year]')

    ax = axes[1]
    ax.plot(r/au, St)
    ax.set_yscale('log')
    ax.set_title('Stokes number')

    ax = axes[2]
    ax.plot(r/au, Q)
    ax.set_yscale('log')
    ax.set_title('Toomre Q')

    ax = axes[3]
    ax.plot(r/au, - dvphi, label=r'$- v_{\phi}$')
    ax.plot(r/au, dvr, label=r'$v_{r}$')
    ax.legend()
    ax.set_yscale('log')
    ax.set_title('relative velocity [cm/s]')

    for ax in axes:
        ax.set_xlabel('r [au]')
        ax.set_xscale('log')

    fig.tight_layout()
    plt.show()

    pdb.set_trace()

if __name__ == '__main__':
    main()

