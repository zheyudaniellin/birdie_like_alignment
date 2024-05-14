"""
plot the alignment angle of the general spheroid

"""
import numpy as np
import matplotlib.pyplot as plt
import model_tensor2
rad = np.pi / 180
alphabet = 'abcdefghijklmn'
lincol = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']

def plot1(psi, abar, aln, figname=None):
    fig, ax = plt.subplots(figsize=(5,5))
    for i in range(len(abar)):
        label = r'$\bar{a}$=%.1f'%abar[i]
        ax.plot(psi, aln[i,:], label=label, color=lincol[i], linewidth=2)

    ax.set_xlabel(r'$\psi$ [$^{\circ}$]')
    ax.set_xlim(0, 90)
    ax.set_xticks(np.arange(0, 91, 15))

    ax.set_ylabel(r'$\theta_{\text{align}}$ [$^{\circ}$]')
    ax.set_ylim(-90, 0)
    ax.set_yticks(np.arange(-90, 1, 15))

    # plot a diagonal
    ax.plot([0,90], [0, -90], color='k', linestyle=':')

    ax.set_aspect('equal')
    ax.legend()

    fig.tight_layout()

    if figname is not None:
        fig.savefig(figname)

    plt.show()

def plot2(psi, abar, aln, aln_gbar, figname=None):
    """
    plot the alignment angle and also the angle for g
    """
    fig, axes = plt.subplots(2, 1, figsize=(5,8), squeeze=True)

    ax = axes[0]
    for i in range(len(abar)):
        label = r'$\bar{a}$=%.1f'%abar[i]
        ax.plot(psi, aln[i,:], label=label, color=lincol[i], linewidth=2)

    ax.set_xlabel(r'$\psi$ [$^{\circ}$]')
    ax.set_xlim(0, 90)
    ax.set_xticks(np.arange(0, 91, 15))

    ax.set_ylabel(r'$\theta_{\text{align}}$ [$^{\circ}$]')
    ax.set_ylim(-90, 0)
    ax.set_yticks(np.arange(-90, 1, 15))

    # plot a diagonal
    ax.plot([0,90], [0, -90], color='k', linestyle=':')

    ax.legend(loc='lower left')

    # gbar
    ax = axes[1]
    for i in range(len(abar)):
        label = r'$\bar{a}$=%.1f'%abar[i]
        ax.plot(psi, aln_gbar[i,:], label=label, color=lincol[i], linewidth=2)

    ax.set_xlabel(r'$\psi$ [$^{\circ}$]')
    ax.set_xlim(0, 90)
    ax.set_xticks(np.arange(0, 91, 15))

    ax.set_ylabel(r'$\theta_{\text{g}}$ [$^{\circ}$]')
#    ax.set_ylim(-90, 0)
#    ax.set_yticks(np.arange(-90, 1, 15))

    # labels
    for i, ax in enumerate(axes):
        txt = '({:s})'.format(alphabet[i])
        ax.text(0.98, 0.98, txt, va='top', ha='right', 
                transform=ax.transAxes)

    fig.tight_layout()

    if figname is not None:
        fig.savefig(figname)

    plt.show()

def main():
    # ==== settings ====
    abar = np.array([0.5, 0.9, 1, 1.1, 2])
    psi = np.linspace(0, 90, 50)

    rhog = 1
    vth = 1
    v = 1
    c = 1
    rhos = 1
    g = 0.01

    # calculate the alignment angle
    Ks = np.zeros([len(abar), len(psi)])
    Kc = np.zeros([len(abar), len(psi)])

    for i in range(len(abar)):
        for j in range(len(psi)):
            coef = model_tensor2.coefficients(rhog, vth, v, c, abar[i], rhos, g, psi=psi[j]*rad)
            Ks[i,j] = coef['Ks']
            Kc[i,j] = coef['Kc']

    # alignment angle
    aln = np.arctan2(-Kc, Ks) / rad

    # angle of the offset
    aln_gbar = aln + psi[None,:]

    # ==== plotting ====
    figname = 'results/theta_align.pdf'
#    plot1(psi, abar, aln, figname=figname)

    plot2(psi, abar, aln, aln_gbar, figname=figname)


if __name__ == '__main__':
    main()

