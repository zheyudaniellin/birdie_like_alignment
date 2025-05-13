"""
prob_phase.py

calculate the phase portrait 
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import phaseportrait
import doublesphere
from natconst import kk, mp
import natconst
import pdb
lincol = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']

rad = np.pi / 180

# ==== dynamical functions ====
def func(theta, omega, *, D=1, K=1, I=1):
    """
    the system of equation:
    I theta'' + D theta' + K sin(theta) = 0

    don't do any normalization
    """

    return [omega, - (D * omega + K * np.sin(theta)) / I]

def func2(hat_theta, hat_omega, *, D=1, K=1, I=1):
    """
    The system of equations specifically of the form
    I theta'' + D theta' + K sin(theta) = 0
    We normalized the equation. 
    omega_0 = sqrt(K / I)
    The two first order ODE are:
    d hat_theta / d hat_t = hat_omega
    d hat_omega / d hat_t = - 2 pi (D / sqrt(K I) * hat_omega + sin(2 pi hat_theta)
    """
    omega_0 = np.sqrt(K / I)
    return [omega_0 / 2 / np.pi * hat_omega, 
            -D/I * hat_omega - np.sqrt(K/I) * np.sin(2*np.pi*hat_theta)
            ]
#    return [hat_omega, - 2 * np.pi * (D/np.sqrt(K*I) * hat_omega + np.sin(2*np.pi*hat_theta))]

# ==== plotting====
def create_constant_cmap():
    """
    """
    N = 256
    vals = np.ones((N, 4))
    vals[:,:3] *= 0
    newcmp = mpl.colors.ListedColormap(vals, name='my_constant')
#    plt.register_cmap('my_constant', newcmp)
    mpl.colormaps.register(cmap=newcmp)

def plot1(I, D, K, hat_theta_lim=[-1,1], hat_omega_lim=[-3,3], figname=None):
    """
    get the phase diagram given the coefficients for func

    hat_theta_lim : list of 2 floats
        this is theta normalized by 2 pi
    hat_omega_lim : list of 2 floats
        this is omega normalized by omega_0
    """
    omega_0 = np.sqrt(K / I)

    theta_lim = 2*np.pi * np.array(hat_theta_lim)
    omega_lim = omega_0 * np.array(hat_omega_lim)

    # we want the customized cmap
    cmap = plt.colormaps['my_constant']

    # get the phase diagram
    pp = phaseportrait.PhasePortrait2D(func,
        [theta_lim, omega_lim],
        dF_args={'D':D, 'K':K, 'I':I}, 
        density=2,
        MeshDim=30,
        color=cmap
        )

    fig, ax = pp.plot()

    # xlabel. originally in radians
    deg_theta = np.arange(-360, 361, 90)
    rad_theta = deg_theta * rad
    ax.set_xticks(rad_theta)
    ax.set_xticklabels(['%d'%ival for ival in deg_theta])
    ax.set_xlabel(r'$\theta$ $[^{\circ}]$')

    # ylabel. originally in radians/sec
    # but we want to normalize it
    yticklabels = np.arange(-3, 3.1, 1)
    ax.set_yticks(yticklabels * omega_0)
    ax.set_yticklabels(['%d'%ival for ival in yticklabels])
    ax.set_ylabel(r'$\omega / \omega_{0}$')
#    ax.set_title('Double-sphere phase portrait')

    # special points
    # attractors
    for ival in [-360, 0, 360]:
        ax.plot([ival*rad], [0], 'ko')
    # unstable points
    for ival in [-180, 180]:
        ax.plot([ival*rad], [0], 'kx')

    # save the figure
    fig.set_size_inches(12, 7)
    if figname is not None:
        fig.savefig(figname)

    plt.show()

def plot2(I, D, K, hat_theta_lim=[-1,1], hat_omega_lim=[-3,3], figname=None):
    """
    get the phase diagram given the coefficients for func2
    """
    # we want the customized cmap
    cmap = plt.colormaps['my_constant']

    # get the phase diagram
    pp = phaseportrait.PhasePortrait2D(func2,
        [hat_theta_lim, hat_omega_lim],
        dF_args={'D':D, 'K':K, 'I':I}, 
        density=2,
        MeshDim=50,  # 40
        color=cmap)

#    fig, ax = pp.plot(color=plt.cm.viridis)
    fig, ax = pp.plot()

    # plot the limiting separatrix
    x_sep_hat = np.linspace(hat_theta_lim[0], hat_theta_lim[1], 300)
    x_sep_rad = x_sep_hat * 2*np.pi
    y_sep_hat = np.sqrt(2 * (1 + np.cos(x_sep_rad)))
    for _sign in [-1, 1]:
        ax.plot(x_sep_hat, _sign * y_sep_hat, color=lincol[2], zorder=0.5, linewidth=3)

    # xlabel
#    deg_theta = np.arange(-hat_theta_lim[0]*360, hat_theta_lim[1]*360+1, 90)
    deg_theta = np.array([-360, -180, 0, 180, 360])
    xticks = deg_theta / 360.0
    ax.set_xticks(xticks)
    ax.set_xticklabels(['%d'%itick for itick in deg_theta])
    ax.set_xlabel(r'$\theta$ [$^{\circ}$]')

    # ylabel
    ax.set_ylabel(r'$\omega / \omega_{o}$')
#    ax.set_title('Double-sphere phase portrait')

    # also show omega with respect to the thermal rate
    """
    ax2 = ax.twinx()
    ax2.set_ylim(hat_omega_lim*timescale['w_o']/w_t)
    ax2.set_ylabel(r'$\omega / \omega_{T}$')
    """

    # special points
    # attractors
    for ival in [-1, 0, 1]:
        ax.plot([ival], [0], 'o', color='k', markersize=10)

    # unstable points
    for ival in [-0.5, 0.5]:
        ax.plot([ival], [0], '+', markersize=15, mew=2, markeredgecolor='k')

    # save the figure
    fig.set_size_inches(9, 6)
    fig.tight_layout()

    if figname is not None:
        fig.savefig(figname)

    plt.show()

# ==== panel ====
def main():
    # ==== settings ====
    hat_theta_lim = np.array([-1, 1])
    hat_omega_lim = np.array([-1, 1]) * 3.5

    # grain parameters
    l = 0.1 # [cm]
    rhos = 3.0
    epsilon = 1.0
    kappa = 1.01
    lam = 1.0
    rhog = 1e9 * (2.3 * mp)

    T = 100 # [kelvin]
    W = 2000.0 # [cm/s]
    
    # thermal velocity
#    vth = np.sqrt(8 * kk * T / 2.3 / mp)
    vth = 20000.0 # [cm/s]

    # ==== calculation ====
    # calculate the grain model
    coeff = doublesphere.coefficients(l, rhos, epsilon, kappa, lam, rhog, vth, W)
    timescale = doublesphere.timescales(coeff['I'], coeff['D'], coeff['P'])

    # print the damping ratio
    print('damping ratio = %.3e'%timescale['q'])

    # thermal rotation rate
    w_t = np.sqrt(2 * kk * T / coeff['I'])

    # now calculate the portrait
    figname = 'results/doublesphere_phase_portrait.pdf'
#    plot1(coeff['I'], coeff['D'], coeff['P'], figname=figname)
    plot2(coeff['I'], coeff['D'], coeff['P'], hat_theta_lim=hat_theta_lim, hat_omega_lim=hat_omega_lim, figname=figname)

def main2():
    """
    create another phase portrait where the damping is stronger
    """
    # ==== settings ====
    hat_theta_lim = np.array([-1, 1])
    hat_omega_lim = np.array([-1, 1]) * 3.5

    # grain parameters
    l = 0.1 # [cm]
    rhos = 3.0
    epsilon = 1.0
    kappa = 1.01
    lam = 1.0
    rhog = 1e12 * (2.3 * mp)

    T = 100 # [kelvin]
    W = 5e-5 # [cm/s]

    # thermal velocity
#    vth = np.sqrt(8 * kk * T / 2.3 / mp)
    vth = 20000.0 # [cm/s]

    # ==== calculation ====
    # calculate the grain model
    coeff = doublesphere.coefficients(l, rhos, epsilon, kappa, lam, rhog, vth, W)
    timescale = doublesphere.timescales(coeff['I'], coeff['D'], coeff['P'])

    # print the damping ratio
    print('damping ratio = %.3e'%timescale['q'])

    # thermal rotation rate
    w_t = np.sqrt(2 * kk * T / coeff['I'])

    # now calculate the portrait
    figname = 'results/doublesphere_phase_portrait_2.pdf'
    plot2(coeff['I'], coeff['D'], coeff['P'], hat_theta_lim=hat_theta_lim, hat_omega_lim=hat_omega_lim, figname=figname)


if __name__ == '__main__':
    # create the cmap before calling any pipeline
    create_constant_cmap()

    # ==== pipelines ====
    main()
    main2()

