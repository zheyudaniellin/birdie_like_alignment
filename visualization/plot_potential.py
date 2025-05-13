"""
plot_potential.py

show a schematic of the potential from the restoring torque
"""
import numpy as np
import matplotlib.pyplot as plt
rad = np.pi / 180

def main():
    # ==== settings ====
    K = 1.0
    thet = np.linspace(-180, 180, 101)
    U = - K * np.cos(thet*rad)

    fig = plt.figure(figsize=(9, 4))
    ax = fig.gca()
    ax.plot(thet, U, color='k', linestyle='-', linewidth=2)

    # add the potential range
    ax.errorbar([thet[-1]+5],[0], yerr=[K], color='C1', capsize=5, linewidth=2)
    ax.text(thet[-1], 0, r'$2P$', ha='right', va='center')

    # x-axis
    ax.set_xticks(np.arange(thet[0], thet[-1]+1, 90))
    ax.axvline(0, color='k', linestyle=':')
#    ax.plot([0,0], [-K,K], color='k', linestyle=':')
    for ithet in [-90, 90]:
        ax.axvline(ithet, color='k', linestyle='--')

    ax.set_xlim(thet[0], thet[-1]+10)
    ax.set_xlabel(r'$\theta$ [$^{\circ}$]')

    # y-axis
#    ax.get_yaxis().set_visible(False)
    ax.set_yticks([])
#    ax.set_yticklabels(['-K', '0', 'K'])
    ax.set_ylabel(r'Potential Energy $U$')
    ax.set_ylim(-K*1.05, K*1.05)

    ax.spines[['right', 'top', 'left']].set_visible(False)

    fig.tight_layout()

    figname = 'results/potential.pdf'
    fig.savefig(figname)
    plt.show()

if __name__ == '__main__':
    main()

