"""
plot the dimensionless time factors
"""
import numpy as np
import matplotlib.pyplot as plt
import natconst
import model_tensor1
import pdb
alphabet = 'abcdefghijklmn'
lincol = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']

def plot1(abar, gbar, tdamp, tosc, figname=None):
    """
    plot the two timescales as a function of gbar and iterate for different abar

    Assume that number of abar is odd and the central one is the spherical case
    """
    def get_color(abar):
        """
        get the color depending on abar
        """
        if abar == 1:
            color = 'k'
        else:
            color = None
        return color

    def get_linestyle(abar):
        """
        """
        if abar < 1:
            linestyle = '--'
        elif abar == 1:
            linestyle = '-'
        else:
            linestyle = '-.'
        return linestyle

    # create figure
    fig, axes = plt.subplots(2, 1, sharex=False, sharey=False, squeeze=True,
            figsize=(5, 8))

    # plot damping time
    ax = axes[0]
    for i in range(len(abar)):
        col = get_color(abar[i])
        if col is None:
            col = lincol[i]
        linestyle = get_linestyle(abar[i])
        label = r'$\bar{a}$=%.1f'%abar[i]
        ax.plot(gbar, tdamp[i,:], label=label, 
                color=col,
                linestyle=linestyle, 
                linewidth=2)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'$\hat{t}_{d,s}$')
    ax.legend()

    # plot oscillation time
    ax = axes[1]
    for i in range(len(abar)):
        col = get_color(abar[i])
        if col is None:
            col = lincol[i]
        linestyle = get_linestyle(abar[i])
        label = r'$\bar{a}$=%.1f'%abar[i]
        ax.plot(gbar, tosc[i,:], label=label, 
                color=col,
                linestyle=linestyle, 
                linewidth=2)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'$\hat{t}_{o,s}$')
    ax.legend()

    for i, ax in enumerate(axes):
        ax.set_xlabel(r'$\bar{g}$')
        ax.set_xlim(gbar[0], gbar[-1])
        ax.text(0.03, 0.98, '({:s})'.format(alphabet[i]), va='top', ha='left', transform=ax.transAxes)

    fig.tight_layout()
    if figname is not None:
        fig.savefig(figname)

    plt.show()

def main():
    # ==== settings ====
    c = 0.1 # [cm]
#    abar = np.linspace(0.5, 0.9, 5)
    abar = np.array([0.5, 0.9, 1.0, 1.1, 2])
    #gbar = np.geomspace(1e-3, 0.5, 5)
    gbar = np.concatenate((
            np.geomspace(5e-3, 0.1, 50, endpoint=False),
            np.linspace(0.1, 0.5, 50)
            ))

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

    # plot the two profiles
    figname = 'results/spheroid_time_factor.pdf'
    plot1(abar, gbar, t_damp/t_stop, t_osc / t_osc_c, figname=figname)
    pdb.set_trace()

if __name__ == '__main__':
    main()

