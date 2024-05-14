"""
calculate the timescale for an ellipsoid
"""
import pdb
import numpy as np
import matplotlib.pyplot as plt
import natconst
from natconst import au, year
import diskvel_alpha
#import Qconst as disk
import pwldisk as disk
#import ringdisk as disk
import sys
sys.path.append('../ellipsoid')
import model_tensor1
rad = np.pi / 180

def plot1(r, t_kep, t_d, t_o):
    """
    """
    fig, ax = plt.subplots()
    pltx = r / au

    ax.plot(pltx, t_kep/year, 'k:', label='Keplerian')
    ax.plot(pltx, t_d/year, 'k', label='damping')
    ax.plot(pltx, t_o/year, 'k--', label='oscillation')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('r [au]')
    ax.set_ylabel('time [year]')
    ax.set_title('Timescales')
    ax.legend()

    fig.tight_layout()
    plt.show()

def plot2(r, t_kep, t_d, t_o, ang, St, v, d_vr, d_vphi):
    """
    row 1: time
    row 2: angle
    row 3: St
    """
    fig, axgrid = plt.subplots(2,2, sharex=True,sharey=False,
            squeeze=False, figsize=(8,5))
    axes = axgrid.flatten()

    pltx = r / au

    ax = axes[0]
#    ax.plot(pltx, t_kep/year, 'k:', label='Keplerian')
#    ax.plot(pltx, t_d/year, 'k', label='damping')
#    ax.plot(pltx, t_o/year, 'k--', label='oscillation')
#    ax.set_ylabel('time [year]')

    ax.plot(pltx, t_d/t_kep, 'k', label='damping')
#    ax.plot(pltx, t_o/t_kep, 'k--', label='oscillation')
    ax.set_ylabel('t_d / t_kep')

    ax.set_yscale('log')
    ax.legend()

    # angle
    ax = axes[1]
    ax.plot(pltx, ang)
    ax.set_ylabel('polarization angle')

    # St
    ax = axes[2]
    ax.plot(pltx, St)
    ax.set_yscale('log')
    ax.set_ylabel('St')

    # drift velocity
    ax = axes[3]
    ax.plot(pltx, v, 'k', lw=2, label='v')
    ax.plot(pltx, d_vr, label=r'$\delta v_{R}$')
    ax.plot(pltx, -d_vphi, label=r'-$\delta v_{\phi}$')
    ax.set_yscale('log')
    ax.set_ylabel('wind [cm/s]')
    ax.legend()

    # xlabel
    ax.set_xlabel('r [au]')
    ax.set_xscale('log')

    fig.tight_layout()
    plt.show()

def main():
    # ==== settings ====
    r = np.geomspace(1, 100, 100) * au
    alpha = 1e-2

    # details of grain
    c = 0.1
    abar = 0.9
    rho_s = 3.0
    gbar = 0.01

    # ==== disk profiles ====
    rho_g = disk.midplane_density(r)
    T = disk.temperature(r)
    vth = disk.vth(r)
    wk = disk.omega(r)
    t_kep = 2 * np.pi / wk
    beta = disk.beta(r)

    # drift
    hr = disk.pressure_H(r) / r
    vk = wk * r

    tstop = rho_s * c / rho_g / vth
    St = tstop * wk
    d_vr = diskvel_alpha.rel_v_r(St, alpha, beta, hr, vk)
    d_vphi = diskvel_alpha.rel_v_phi(St, alpha, beta, hr, vk)
    v = np.sqrt(d_vr**2 + d_vphi**2)

    # wind angle
    ang = diskvel_alpha.get_angle(d_vr, d_vphi) / rad
    ang = np.mod(ang, 180)

    # apply the grain model
    coef = model_tensor1.coefficients(rho_g, vth, v, c, abar, rho_s, gbar)
    time = model_tensor1.timescales(coef['I'], coef['D'], coef['K'])

    # ==== plotting ====
    plot1(r, t_kep, time['t_d'], time['t_o'])
    plot2(r, t_kep, time['t_d'], time['t_o'], ang, St, v, d_vr, d_vphi)
    pdb.set_trace()

if __name__ == '__main__':
    main()

