"""
diskvel_alpha.py
Place to keep the equations for disk velocity structure
"""
import numpy as np
import matplotlib.pyplot as plt
from natconst import gg, ms, au
rad = np.pi / 180

# ==== gas velocity ====
def get_v_g_r(alpha, hr, vk):
    """
    gas velocity in the radial direction
    """
    return - alpha * hr**2 * vk

def get_v_g_phi(beta, hr, vk):
    """
    gas velocity in the azimuthal direction
    """
    eta = - beta * hr**2
    return vk * np.sqrt(1 - eta)

# ==== relative velocity
def rel_v_r(St, alpha, beta, hr, vk):
    """
    the velocity of the gas relative to the grain in the radial direction

    Parameters
    ----------
    beta : float, ndarray
        the pressure gradient d ln P / d ln r
    hr : float, ndarray
        the height/radius ratio
    vk : float, ndarray
        the keplerian velocity
    """
    diff = - (alpha * St + beta) / (St + 1./St) * hr**2 * vk
    return diff

def rel_v_phi(St, alpha, beta, hr, vk):
    """
    the velocity of the gas relative to the grain in the azimuthal direction

    Parameters
    ----------
    beta : float, ndarray
        the pressure gradient d ln P / d ln r
    """
    diff = 0.5 * (-alpha + beta * St) / (St + 1./St) * hr**2 * vk
    return diff

# ==== dust velocity ====
def get_v_d_r(St, alpha, beta, hr, vk):
    # gas velocity
    v_g_r = get_v_g_r(alpha, hr, vk)

    # drift
    d_vr = rel_v_r(St, alpha, beta, hr, vk)

    # dust 
    v_d_r = v_g_r - d_vr

    return v_d_r

def get_v_d_phi(St, alpha, beta, hr, vk):
    # gas
    v_g_phi = get_v_g_phi(beta, hr, vk)

    # drift
    d_vphi = rel_v_phi(St, alpha, beta, hr, vk)

    # dust velocity
    v_d_phi = v_g_phi - d_vphi

    return v_d_phi

def get_angle(vec_r, vec_phi):
    """
    calculate the angle from the radial direction and increases towards the azimuthal direction for any vector

    results are in radians
    """
    ang = np.arctan2(vec_phi, vec_r)

    return ang

# ==== plotting ====
def plot1(St, alpha, d_vr, d_vphi, hr, vk):
    """
    plot the relative velocities only
    as a function of St
    """
    nrow = len(alpha)
    ncol = 3
    fig, axgrid = plt.subplots(nrow,ncol,sharex=False, sharey=False, 
            squeeze=False, figsize=(12, 7))
    axes = axgrid.flatten()

    # iterate
    for i in range(len(alpha)):
        # plot radial 
        ax = axgrid[i,0]
        ax.plot(St, d_vr[:,i] / vk / hr**2, color='k')
        ax.axhline(y=0, color='k', linestyle=':')
#        ax.set_yscale("log")

        # plot azimuthal
        ax = axgrid[i,1]
        ax.plot(St, d_vphi[:,i] / vk / hr**2, color='k')
#        ax.set_yscale("log")

        # plot the angle
        ax = axgrid[i,2]
        ang = get_angle(d_vr[:,i], d_vphi[:,i]) / rad
        ax.plot(St, ang, color='k')

    for ax in axes:
        ax.axvline(x=1, linestyle=':', color='k')
        ax.set_xscale("log")
#        ax.set_xlabel("St")
        ax.set_xlim(St[0], St[-1])

    axgrid[0,0].set_title(r'$\delta v_{r}$')
    axgrid[0,1].set_title(r'$\delta v_{\phi}$')
    axgrid[0,2].set_title(r'$\psi$ [deg]')

    for ax in axgrid[-1,:]:
        ax.set_xlabel("St")
    for i, ax in enumerate(axgrid[:,0]):
        ax.set_ylabel(r'$\alpha$ = %.1e'%alpha[i])

    fig.tight_layout()
    plt.show()

def plot2(St, alpha, d_vr, d_vphi):
    """
    plot the relative velocities only
    as a function of St
    """
    ax = plt.gca()
    for i in range(len(alpha)):
        col = 'C%d'%i
        tag = r'$\Delta v_{\phi}$, $\alpha$=%.1e'%alpha[i]
        ax.plot(St, d_vphi[:,i] / vk / hr**2, label=tag, color=col)

        tag = r'$\Delta v_{r}$, $\alpha$=%.1e'%alpha[i]
        ax.plot(St, d_vr[:,i] / vk / hr**2, label=tag, color=col, linestyle='--')

    ax.axvline(x=1, color='k', linestyle=':')
    ax.legend()
    ax.set_xlabel('St')
    ax.set_ylabel('relative velocity [vk (h/r)^2]')
    ax.set_xscale('log')
#    ax.set_yscale('log')

    plt.show()

def plot_vels(St, v_g_r, v_g_phi, d_vr, d_vphi, hr, vk, title=None):
    """
    plot the drift velocity and the actual gas and dust velocities

    plot as a function of St

    v_g_r : float
        the radial velocity of the gas. doesn't depend on St
    v_g_phi : float
    """
    # dust velocity
    v_d_r = v_g_r - d_vr
    v_d_phi = v_g_phi - d_vphi

    # angle of the drift velocity
    ang = get_angle(d_vr, d_vphi)

    # ==== plotting ====
    nrow, ncol = 4, 1
    fig, axgrid = plt.subplots(nrow,ncol,sharex=True,sharey=False, 
            squeeze=False, figsize=(10, 7))
    axes = axgrid.flatten()

    gcol = 'C0'
    dcol = 'C1'
    # radial direction
    ax = axes[0]
    ax.axhline(v_g_r / hr**2 / vk, color=gcol, label='gas')
    ax.plot(St, v_d_r / hr**2 / vk, color=dcol, label='dust')
    ax.axhline(0, color='k', linestyle=':')
    ax.set_ylabel(r'radial v [$(h/r)^2$ vk]')
    ax.legend()

    # azimuthal direction
    ax = axes[1]
    ax.axhline(v_g_phi / vk, color=gcol, label='gas')
    ax.plot(St, v_d_phi / vk, color=dcol, label='dust')
    ax.set_ylabel('azimuthal v [vk]')
    ax.legend()

    # drift velocity
    ax = axes[2]
    ax.plot(St, d_vr / hr**2 / vk, color='C2', label='d_vr')
    ax.plot(St, d_vphi / hr**2 / vk, color='C3', label='d_vphi')
    ax.set_ylabel('drift velocity [$(h/r)^2$ vk]')
    ax.legend()

    # angle
    ax = axes[3]
    ax.plot(St, ang/rad, color='k')
    ax.set_ylabel('angle [deg]')

    ax.set_xscale('log')
    ax.set_xlim(St[0], St[-1])
    axes[-1].set_xlabel('St')

    fig.suptitle(title)
    fig.tight_layout()
    plt.show()

# ==== settings for different pressure gradient ====
def get_pressure_gradient(mode = 'powerlaw'):
    """
    fetch the pressure gradient
    d ln P / d ln r
    also denoted as beta
    """
    if mode == 'powerlaw':
        beta = beta_powerlaw()
    elif mode == 'gaussian':
        beta = bet_gaussian()
    else:
        raise ValueError('mode unknown')
    return beta

def beta_powerlaw():
    """
    pressure gradient for a power-law disk
    """
    # density power-law index
    p = 3.0

    # temperature power-law index
    q = 0.5

    beta = -p - q
    beta = -1.5
    return beta

# ==== main pipeline ====
def main():
    """
    some demonstration
    """
    # ==== settings ====
    hr = 0.1
    vk = np.sqrt(gg * 1*ms / (100 * au))
    St = np.geomspace(1e-3, 1e3, 50)
    alpha = np.array([0, 1e-2, 1e-1])

    # ==== calculate pressure gradient ====
    beta = get_pressure_gradient(mode='powerlaw')

    # ==== calculate the velocity ====
    d_vr = rel_v_r(St[:,None], alpha[None,:], beta, hr, vk)
    d_vphi = rel_v_phi(St[:,None], alpha[None,:], beta, hr, vk)

    v_g_r = get_v_g_r(alpha, hr, vk)
    v_g_phi = get_v_g_phi(beta, hr, vk)

    # ==== plotting ====
    plot1(St, alpha, d_vr, d_vphi, hr, vk)
    for i in range(len(alpha)):
        title = r'$\alpha$=%.1e, $\beta$=%.2f'%(alpha[i], beta)
        plot_vels(St, v_g_r[i], v_g_phi, d_vr[:,i], d_vphi[:,i], hr, vk, title=title) 

if __name__ == '__main__':
    main()

