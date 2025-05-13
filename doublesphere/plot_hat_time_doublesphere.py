"""
plot how the dimensionless asymmetry time factors depend on asymmetry
"""
import numpy as np
import matplotlib.pyplot as plt
import doublesphere
alphabet = 'abcdefghijk'
lincol = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']

def plot1(epsilon, kappa, hat_t_d, hat_t_o, hat_v_t, figname=None):
    """
    plot as a function of epsilon and kappa
    """
    fig, axes = plt.subplots(2, 1, sharex=False, sharey=False, squeeze=True,
            figsize=(5, 8))

    # plot damping time
    ax = axes[0]
    for i in range(len(kappa)):
        label = r'$\kappa$={:.2f}'.format(kappa[i])
        ax.plot(epsilon, hat_t_d[:,i], label=label, color=lincol[i], 
                linewidth=2)

#    ax.set_xscale('log')
#    ax.set_yscale('log')
    ax.set_ylabel(r'$\breve{t}_{d,d}$')
    ax.legend(loc='center right')

    # plot oscillation time
    ax = axes[1]
    for i in range(len(kappa)):
        label = r'$\kappa$={:.2f}'.format(kappa[i])
        ax.plot(epsilon, hat_t_o[:,i], label=label, color=lincol[i], 
                linewidth=2)

#    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'$\breve{t}_{o,d}$')
    ax.legend()

    for i, ax in enumerate(axes):
        ax.set_xlabel(r'$\epsilon$')
        ax.set_xlim(epsilon[0], epsilon[-1])
        ax.text(0.02, 0.97, '({:s})'.format(alphabet[i]), va='top', ha='left', transform=ax.transAxes)

    fig.tight_layout()
    if figname is not None:
        fig.savefig(figname)

    plt.show()

def main():
    # ==== settings ====
    rhog = 1
    vth = 1
    v = 1
    l = 1
    rhos = 1

    lam = 1.0
    epsilon = np.concatenate((
            np.geomspace(1, 1.1, 50, endpoint=False), 
            np.linspace(1.1, 1.5, 50)
            ))
    kappa = np.array([1, 1.01, 1.1, 1.5])

    # ==== calculations ====
    coef = doublesphere.coefficients(l, rhos, epsilon[:,None], kappa[None,:], lam, rhog, vth, v)

    time = doublesphere.timescales(coef['I'], coef['D'], coef['P'])

    hat_t_d = time['t_d']
    hat_t_o = time['t_o']
    hat_q = time['q']

    # threshold drift
    hat_v_t = doublesphere.drift_threshold_factor(epsilon[:,None], kappa[None,:], lam)

    # ==== plotting ====
    figname = 'results/doublesphere_time_factor.pdf'
    plot1(epsilon, kappa, hat_t_d, hat_t_o, hat_v_t, figname=figname)

def calculate_value():
    # ==== settings ====
    rhog = 1
    vth = 1
    v = 1
    l = 1
    rhos = 1
    lam = 1

    epsilon = 1.01
    kappa = 1.0

    # ==== calculations ====
    coef = doublesphere.coefficients(l, rhos, epsilon, kappa, lam, rhog, vth, v)

    time = doublesphere.timescales(coef['I'], coef['D'], coef['P'])

    hat_t_d = time['t_d']
    hat_t_o = time['t_o']
    hat_q = time['q']

    hat_v_t = doublesphere.drift_threshold_factor(epsilon, kappa, lam)

    print(hat_t_d)
    print(hat_t_o)
    print(hat_q)
    print(hat_v_t)

if __name__ == '__main__':
    main()
#    calculate_value()

