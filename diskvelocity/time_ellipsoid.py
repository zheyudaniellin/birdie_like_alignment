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
alphabet = 'abcdefghijklmnopqrstuv'
lincol = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']

def get_mfp(n):
    """
    calculate the mean free path

    n : ndarray
        number density in particles per cm^3
    """
    # collision cross section in cm^2
    sig = 2e-15 # H2

    return 1. / n / sig

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

    xticklabels = [1, 10, 100]
    ax.set_xticklabels(xticklabels)

    fig.tight_layout()
    plt.show()

def plot_column(r, t_kep, t_d, t_o, ang, St, v, d_vr, d_vphi, mfp, figname=None):
    """
    plot things in a colum

    row 1: St
    row 2: gas wind
    row 3: time
    """
    nrow, ncol, figsize = 4, 1, (5, 10)
#    nrow, ncol, figsize = 2, 2, (8, 5)
    fig, axgrid = plt.subplots(nrow,ncol, sharex=True,sharey=False,
            squeeze=False, figsize=figsize)
    axes = axgrid.flatten()

    pltx = r / au

    # mfp
    ax = axes[0]
    ax.plot(pltx, mfp, color='k', lw=2)
    ax.set_yscale('log')
    ax.set_ylabel(r'$\lambda_{mfp}$ [cm]')

    # St
    ax = axes[1]
    ax.plot(pltx, St, color='k', lw=2, label='St')
    ax.set_yscale('log')
    ax.set_ylabel('St')

    # angle
#    ax = axes[1]
#    ax.plot(pltx, ang)
#    ax.set_ylabel('polarization angle')

    # drift velocity
    ax = axes[2]
    ax.plot(pltx, v, color='k', lw=2, label=r'$A$')
    ax.plot(pltx, abs(d_vr), color=lincol[0], linestyle='--', label=r'|$A_{R}$|')
    ax.plot(pltx, abs(d_vphi), color=lincol[1], linestyle=':', label=r'|$A_{\Phi}$|')
    ax.set_yscale('log')
    ax.set_ylabel(r'$A$ [cm/s]')
    ax.legend(loc='lower right')

    # time
    ax = axes[-1]
    ax.plot(pltx, t_kep/year, 'k:', label='Keplerian', lw=2)
    ax.plot(pltx, t_d/year, 'k', label='damping', lw=2)
    ax.plot(pltx, t_o/year, 'k--', label='oscillation', lw=2)
    ax.set_ylabel('time [year]')

    ax.set_yscale('log')
    ax.legend()

    # xlabel
    ax = axes[-1]
    ax.set_xlabel('R [au]')
    ax.set_xscale('log')
    ax.set_xlim(pltx[0], pltx[-1])
    xticks = [1, 10, 100]
    ax.set_xticks(xticks, labels=['{:d}'.format(itick) for itick in xticks])

    # plot labels
    for i, ax in enumerate(axes):
        txt = '({:s})'.format(alphabet[i])
        ax.text(0.02, 0.98, txt, color='k', va='top', ha='left',
                transform=ax.transAxes)

    fig.tight_layout()

    if figname is not None:
        fig.savefig(figname)

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

    # mean free path
    n_g = rho_g / 2.3 / natconst.mp
    mfp = get_mfp(n_g)

    # ==== plotting ====
#    plot1(r, t_kep, time['t_d'], time['t_o'])

    figname = 'results/time_profile_spheroid.pdf'
    plot_column(r, t_kep, time['t_d'], time['t_o'], ang, St, v, d_vr, d_vphi, mfp, figname=figname)

    # ==== print info ====
    # surface density
    sig0 = disk.surface_density(1 * au)
    print('sig0 = %f'%(sig0))

    # beta
    print('beta = %f'%beta[0])
    pdb.set_trace()

if __name__ == '__main__':
    main()

