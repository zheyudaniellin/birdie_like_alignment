"""
prob_phase.py

calculate the phase portrait 
"""

import numpy as np
import matplotlib.pyplot as plt
import phaseportrait
import model_tensor1
from natconst import kk, mp
import natconst
import pdb

rad = np.pi / 180

def func(hat_theta, hat_omega, *, D=1, K=1, I=1):
    """
    The system of equations specifically of the form
    I theta'' + D theta' + K sin(theta) = 0
    We normalized the equation. 
    omega_0 = sqrt(K / I)
    The two first order ODE are:
    d hat_theta / d hat_t = hat_omega
    d hat_omega / d hat_t = - 2 pi (D / sqrt(K I) * hat_omega + sin(2 pi hat_theta)
    """

    return [hat_omega, - 2 * np.pi * (D/np.sqrt(K*I) * hat_omega + np.sin(2*np.pi*hat_theta))]

def main():
    # ==== settings ====
    hat_theta_lim = np.array([-1, 1])
    hat_omega_lim = np.array([-1, 1]) * 3

    # grain parameters
    c = 0.1 # [cm]
    abar = 0.5
    gbar = 0.1
    rhog = 1e-12
    T = 100 # [kelvin]
    v = 1e-3 # [cm/s]
    rho_s = 3.0
    
    # thermal velocity
    vth = np.sqrt(8 * kk * T / 2.3 / mp)

    # ==== calculation ====
    # calculate the grain model
    coeff = model_tensor1.coefficients(rhog, vth, v, c, abar, rho_s, gbar)
    timescale = model_tensor1.timescales(coeff['I'], coeff['D'], coeff['K'])

    # thermal rotation rate
    w_t = np.sqrt(2 * kk * T / coeff['I'])

    # now calculate the portrait

    # get the phase diagram
    pp = phaseportrait.PhasePortrait2D(func, 
        [hat_theta_lim, hat_omega_lim], 
        dF_args={'D':coeff['D'], 'K':coeff['K'], 'I':coeff['I']}, 
        density=2, 
        MeshDim=30)
    fig, ax = pp.plot(color=plt.cm.viridis)

    # xlabel
    hat_theta = ax.get_xticks()
    deg_theta = hat_theta * 360
    ax.set_xticklabels(deg_theta)
    ax.set_xlabel(r'$\theta$ $[^{\circ}]$')

    # ylabel
    ax.set_ylabel(r'$\omega / \omega_{0}$')
    ax.set_title('phase diagram')

    # also show omega with respect to the thermal rate
    ax2 = ax.twinx()
    ax2.set_ylim(hat_omega_lim*timescale['w_o']/w_t)
    ax2.set_ylabel(r'$\omega / \omega_{T}$')

    fig.set_size_inches(12, 7)

    figname = 'results/spheroid_phase_portrait.pdf'
    plt.show()

    pdb.set_trace()
if __name__ == '__main__':
    main()

