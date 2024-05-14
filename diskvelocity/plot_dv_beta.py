"""
plot the relative velocities using diskvel_alpha prescription

I want to show the gas drift velocity as a function of St for different values of alpha and beta

"""
import pdb
import numpy as np
import matplotlib.pyplot as plt
import natconst
from natconst import au
import diskvel_alpha
rad = np.pi / 180
alphabet = 'abcdefghijklmnopqrstuv'

def plot_grid(St, alpha, v_g_r, v_g_phi, v_d_r, v_d_phi, d_vr, d_vphi,
        beta, hr, vk, figname=None):
    """
    row 1: radial direction
    row 2: azimuthal direction
    row 3: gas wind
    row 4: wind angle
    """
    # angle of the drift velocity
    ang = diskvel_alpha.get_angle(d_vr, d_vphi)
    ang = np.mod(ang, np.pi)

    nrow, ncol = 4, len(beta)
    fig, axgrid = plt.subplots(nrow,ncol,sharex=True,sharey='row',
            squeeze=False, figsize=(12,8))
    axes = axgrid.flatten()

    lw = 2
    gas_color = 'C2'
    dust_color = 'k'
    # iterate through alpha
    for i, _beta in enumerate(beta):
        # radial direction
        ax = axgrid[0,i]
        norm = hr**2 * vk
        ax.axhline(v_g_r/norm, color=gas_color, label='gas', 
                linewidth=lw)
        ax.plot(St, v_d_r[:,i]/norm, color=dust_color, label='dust', 
                linewidth=lw)

        # azimuthal direction
        ax = axgrid[1,i]
        norm = hr**2 * vk
        ax.axhline((v_g_phi[i]-vk)/norm, color=gas_color, label='gas', 
                linewidth=lw)
        ax.plot(St, (v_d_phi[:,i]-vk)/norm, color=dust_color, label='dust', 
                linewidth=lw)
        ax.axhline(0, color='k', linestyle=':')

        # gas wind
        ax = axgrid[2,i]
        norm = hr**2 * vk
        ax.plot(St, d_vr[:,i]/norm, label=r'$\delta v_{r}$', color='k', 
                linestyle='--', linewidth=lw)
        ax.plot(St, d_vphi[:,i]/norm, label=r'$\delta v_{\phi}$', color='k', 
                linewidth=lw)
        ax.axhline(0, color='k', linestyle=':')

        # alignment angle
        ax = axgrid[3,i]
        norm = rad
        ax.plot(St, ang[:,i]/norm, color='k', linewidth=lw)

        ax.axhline(-90, color='k', linestyle=':')

    # xscale
    for ax in axes:
        ax.set_xscale('log')
        ax.set_xlim(St[0], St[-1])

    # legends
    axgrid[0,0].legend(loc='lower right')
    axgrid[1,0].legend(loc='center right')
    axgrid[2,0].legend(loc='upper right')

    # x labels
    for ax in axgrid[-1,:]:
        ax.set_xlabel('St')

    # y labels
    axgrid[0,0].set_ylabel(r'$v_{r}$'+'\n'+r'[$(H/R)^{2}$ $v_{K}$]')
    axgrid[1,0].set_ylabel(r'$v_{\phi}-v_{K}$'+'\n'+r'[$(H/R)^{2}$ $v_{K}$]')
    axgrid[2,0].set_ylabel(r'Gas Wind'+'\n'+r'[$(H/R)^{2}$ $v_{K}$]')
    axgrid[3,0].set_ylabel(r'$\chi$'+'\n'+r'[$^{\circ}$]')

    for i, ax in enumerate(axgrid[0,:]):
        ax.set_title(r'$\beta$=%.2f'%beta[i])

    # additional modification for wind angle
#    for ax in axgrid[-1,:]:
#        ax.set_yticks(np.arange(-105, 1, 15))

    # move the ylim a bit for the labels
#    axgrid[0,0].set_ylim(None, 0.2)
#    axgrid[1,0].set_ylim(None, 0.1)
#    axgrid[-1,0].set_ylim(None, 15)

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
    hr = 0.1
    ms = 1 * natconst.ms
    St = np.geomspace(1e-3, 1e3, 50)

    alpha = 0.01
    beta = np.array([0.0])

    # ==== calculate the velocity ====
    vk = np.sqrt(natconst.gg * 1*ms / (100 * au))

    # gas wind
    d_vr = diskvel_alpha.rel_v_r(St[:,None], alpha, beta[None,:], hr, vk)
    d_vphi = diskvel_alpha.rel_v_phi(St[:,None], alpha, beta[None,:], hr, vk)

    # gas
    v_g_r = diskvel_alpha.get_v_g_r(alpha, hr, vk)
    v_g_phi = diskvel_alpha.get_v_g_phi(beta, hr, vk)

    # dust
    v_d_r = diskvel_alpha.get_v_d_r(St[:,None], alpha, beta[None,:], hr, vk)
    v_d_phi = diskvel_alpha.get_v_d_phi(St[:,None], alpha, beta[None,:], hr, vk)
    
    # ==== plotting ====
#    figname = 'results/St_alpha_grid.pdf'
    figname = None
    plot_grid(St, alpha, v_g_r, v_g_phi, v_d_r, v_d_phi, d_vr, d_vphi, beta, hr, vk, figname=figname)

if __name__ == '__main__':
    main()

