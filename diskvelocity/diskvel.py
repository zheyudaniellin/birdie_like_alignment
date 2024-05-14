"""
diskvel.py

Solutions based on Nakagawa which includes dust back reaction on the gas

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from natconst import gg, ms, au
import pdb
rad = np.pi / 180

def get_f(p, q):
    """
    calculate the fraction for eta
    Or it's simply 11/4 for some typical disk
    """
    f = p + q
    return f

def get_eta(hr, p, q):
    """
    calculate the eta value
    """
    f = get_f(p, q)
    eta = f * hr**2
    return eta


def v_d_r(St, d2g):
    """
    dust radial velocity
    This is relative to eta * vk

    Parameters
    ----------
    St : stokes number
    d2g : dust to gas mass ratio
    """
    v = -2 / (St + 1./St * (1 + d2g)**2)
    return v

def v_d_phi(St, d2g):
    """
    dust azimuthal velocity
    """
    v = - (1 + d2g) / ((1 + d2g)**2 + St**2)
    return v

def v_g_r(St, d2g):
    """
    gas radial velocity
    """
    v = 2 * d2g / (St + 1./St * (1 + d2g)**2)
    return v

def v_g_phi(St, d2g):
    """
    gas azimuthal velocity
    """
    v = d2g / (1 + d2g) / (1 + St**2 / (1+d2g)**2) - 1.0
    return v

def rel_v_r(St, d2g):
    """
    relative velocity of the gas from the frame of the dust
    """
    dv = v_g_r(St, d2g) - v_d_r(St, d2g)
    return dv

def rel_v_phi(St, d2g):
    """
    relative velocity of the gas from the frame of the dust
    """
    dv = v_g_phi(St, d2g) - v_d_phi(St, d2g)
    return dv

def plot1(St, d2g):
    # === plotting ====
    # top row: vr and vphi for each dust and gas
    # bottom row: difference of vr and vphi
    fig, axgrid = plt.subplots(2,2,sharex=False, sharey=False,squeeze=False,
            figsize=(10,7))
    axes = axgrid.flatten()

    linestyle = ['-', '--', ':']
    cols = ['C0', 'C1', 'C2']

    for i in range(len(d2g)):
        # radial direction
        ax = axgrid[0,0]
        ax.plot(St, v_d_r(St, d2g[i]), color=cols[i], 
                linestyle=linestyle[0], label=r'$\epsilon$=%.2f'%d2g[i])
        ax.plot(St, v_g_r(St, d2g[i]), color=cols[i], 
                linestyle=linestyle[1], label=r'$\epsilon$=%.2f'%d2g[i])

        # azimuthal direction
        ax = axgrid[0,1]
        ax.plot(St, v_d_phi(St, d2g[i]), color=cols[i],
                linestyle=linestyle[0], label=r'$\epsilon$=%.2f'%d2g[i])
        ax.plot(St, v_g_phi(St, d2g[i]), color=cols[i],
                linestyle=linestyle[1], label=r'$\epsilon$=%.2f'%d2g[i])

        # difference in radial direction
        ax = axgrid[1,0]
        ax.plot(St, rel_v_r(St, d2g[i]), color=cols[i],  
                linestyle=linestyle[0], label=r'$\epsilon$=%.2f'%d2g[i])

        # difference in azimuthal direction
        ax = axgrid[1,1]
        ax.plot(St, rel_v_phi(St, d2g[i]), color=cols[i], 
                linestyle=linestyle[0], label=r'$\epsilon$=%.2f'%d2g[i])

    # common legend element
    legend_element = [
            Line2D([0],[0],color='k', linestyle=linestyle[0], label='dust'),
            Line2D([0],[0],color='k', linestyle=linestyle[1], label='gas'),
            ]
    for i in range(len(d2g)):
        legend_element.append(
                Line2D([0],[0], color=cols[i], label=r'$\epsilon$=%.2f'%d2g[i])
                    )

    ax = axgrid[0,0]
    ax.set_ylabel('radial velocity [$\eta v_{K}$]')
    ax.legend(handles=legend_element, )

    ax = axgrid[0,1]
    ax.set_ylabel('azimuthal velocity from $v_{K}$ [$\eta v_{K}$]')
    ax.legend(handles=legend_element, )

    ax = axgrid[1,0]
    ax.set_ylabel('radial velocity difference [$\eta v_{K}$]')
    ax.legend()

    ax = axgrid[1,1]
    ax.set_ylabel('azimuthal velocity difference [$\eta v_{K}$]')
    ax.legend()

    for ax in axes:
        ax.set_xlabel('St')
        ax.axvline(x=1, color='k', linestyle=':')
        ax.axhline(y=0, color='k', linestyle=':')
        ax.set_xscale('log')

    fig.tight_layout()
    plt.show()

def plot2(St, d2g):
    """
    calculate the angle of the polarization 
    """
    cols = ['C0', 'C1', 'C2']

    fig = plt.figure()
    ax = fig.gca()
    for i in range(len(d2g)):
        dv_r = rel_v_r(St, d2g[i])
        dv_phi = rel_v_phi(St, d2g[i])
        ang = np.arctan2(dv_phi, dv_r) / rad
        ax.plot(St, ang, color=cols[i], label=r'$\epsilon$=%.2f'%d2g[i])

    ax.legend()
    ax.set_xscale('log')
    ax.set_xlabel('St')

    ax.set_ylabel('angle [degrees]')
    ylim = ax.get_ylim()
    ax.set_yticks(np.arange(0, -91, -15))

    ax.set_title('angle from the radial direction')

    fig.tight_layout()
    plt.show()

def main():
    """
    some demonstration
    """
    St = np.geomspace(1e-3, 1e2, 100)
    d2g = np.array([1e-2, 0.1, 1])

    # ==== plotting ====
    plot1(St, d2g)
    plot2(St, d2g)

if __name__ == '__main__':
    main()

