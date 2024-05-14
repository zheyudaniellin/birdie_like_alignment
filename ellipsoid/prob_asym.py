"""
prob_asym.py

get a feel of how much asymmetry of the grain we need for the damping time and oscillation time

The asymmetry here means abar or gbar
"""
import numpy as np
import matplotlib.pyplot as plt
import natconst
import model_tensor1
import pdb

def plota(abar, gbar, t, tnorm):
    """
    plot the time as a function of abar for different gbar
    Normalize the time by tnorm. 
    Don't plot. just return the ax
    """
    ax = plt.gca()
    for i in range(len(gbar)):
        label = r'$\bar{g}$ = %.3f'%(gbar[i])
        ax.plot(abar, t[:,i]/tnorm, label=label)

    ax.set_xlabel(r'$\bar{a}$')
    ax.legend()

    return ax

def plotg(abar, gbar, t, tnorm):
    """
    plot the time as a function of gbar for different abar
    Normalize the time by tnorm.
    Don't plot. just return the ax
    """
    ax = plt.gca()
    for i in range(len(abar)):
        label = r'$\bar{a}$ = %.3f'%(abar[i])
        ax.plot(gbar, t[i,:]/tnorm, label=label)

    ax.set_xlabel(r'$\bar{g}$')
    ax.legend()

    return ax

def main():
    # ==== settings ====
    c = 0.1 # [cm]
    abar = np.linspace(0.5, 2.0, 50)
#    gbar = np.geomspace(1e-3, 0.5, 5)
    gbar = np.array([0, 0.01, 0.05, 0.1, 0.5])
    rho_s = 3.0

    # environmental inputs
    # these actually don't matter, since we'll normalize it later
    rho_g = 1e-15
    T = 50 # [Kelvin]
    vth = np.sqrt(8 * natconst.kk * T / np.pi / natconst.mp)
    dv = 0.1 # [cm/s]

    # ==== calculations ====
    # stopping time
    t_stop = rho_s / rho_g * c / vth

    # characteristic oscillation time
    t_osc_c = np.sqrt(rho_s * c**2 / rho_g / dv / vth)

    # calculate the grain
    t_damp = np.zeros([len(abar), len(gbar)])
    t_osc = np.zeros_like(t_damp)

    for i in range(len(abar)):
        for j in range(len(gbar)):
            coef = model_tensor1.coefficients(rho_g, vth, dv, c, abar[i], rho_s, gbar[j])
            time = model_tensor1.timescales(coef['I'], coef['D'], coef['K'])

            t_damp[i,j] = time['t_d']
            t_osc[i,j] = time['t_o']

    # ==== plotting ====
    # plot the damping time
    ax = plota(abar, gbar, t_damp, t_stop)
    ax.set_title(r't_damp / t_stop')
    ax.set_yscale('log')
    plt.show()

    # plot the damping time
    ax = plotg(abar, gbar, t_damp, t_stop)
    ax.set_title(r't_damp / t_stop')
    ax.set_yscale('log')
    plt.show()

    # plot the oscillation time
    ax = plota(abar, gbar, t_osc, t_osc_c)
    ax.set_title(r't_osc / t_c')
    ax.set_yscale('log')
    plt.show()

    # plot the oscillation time
    ax = plotg(abar, gbar, t_osc, t_osc_c)
    ax.set_title(r't_osc / t_c')
    ax.set_yscale('log')
    plt.show()


    pdb.set_trace()

if __name__ == '__main__':
    main()

