"""
check_ode.py

compare the evolution using the grid-based model and the tensor model
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, dblquad
import pdb
import model_integ
import model_tensor1
from natconst import au, year, ms, gg, mp, kk
rad = np.pi / 180

# ==== plotting ====
def plot1(sols, tags, style):
    """
    plot the evolution of theta and omega

    sols : list of dict
        each dict should have 
        t = time
        theta, omega
    """
    nrow, ncol = 2, 1
    fig, axgrid = plt.subplots(nrow,ncol,sharex=False,sharey=False,
        squeeze=False)
    axes = axgrid.flatten()

    for i in range(len(sols)):
        t = sols[i]['t']
        theta = sols[i]['theta']
        omega = sols[i]['omega']

        # ==== theta ====
        ax = axes[0]
        ax.plot(t/year, theta/rad, style[i], label=tags[i])
#        ax.set_yticks(np.arange(-180, 181, 45))

        # ==== omega ====
        ax = axes[1]
        ax.plot(t/year, omega/rad, style[i], label=tags[i])

    # ==== labels ====
    ax = axes[0]
    ax.set_ylabel(r'$\theta$ [$^{\circ}$]')
    ax.axhline(y=0, color='k', linestyle=':')
    ax.legend()

    ax = axes[1]
    ax.set_ylabel(r'$\omega$')
    ax.axhline(y=0, color='k', linestyle=':')

    for ax in axes:
        ax.set_xlabel('t [years]')
        ax.set_xlim(0, None)

    fig.tight_layout()
    plt.show()

def plot2(t, theta, omega):
    """
    plot the phase
    """
    fig = plt.figure()

    ax = fig.gca()

    ax.plot(theta/rad, omega/rad)

    ax.set_xlabel(r'$\theta$ [$^{\circ}$]')

    ax.set_ylabel(r'$\omega$')

    fig.tight_layout()
    plt.show()


# ==== integration function ====
def fn(y, t, a, c, Iy, gb, vg, N, m, T):
    """
    the system of equation
    """
    theta, omega = y

    # calculate the torque
    args = (a, c, gb, theta, omega, vg, N, m, T)
#    res = dblquad(model_integ.integrand, 0, 2*np.pi, -1,1, args=args)
#    torque = res[0] * a * c
    torque = model_integ.do_grid(a, c, gb, theta, omega, vg, N, m, T, nphi=20, nw=10)

    dydt = [omega, torque / Iy]

    return dydt

def fn2(y, t, I, D, K):
    """
    system of equation if we know the coefficients
    """
    theta, omega = y
    dydt = [omega, - D/I * omega - K/I * np.sin(theta)]
    return dydt

# =================================
def main():
    # ==== settings ====
    # ellipsoid shape
    c = 0.1 # [cm]
    a = 0.9 * c
    rho_s = 3.0

    # location of the center of mass in body frame
    gb = -0.01 * c

    # velocity of gas
#    vg = 10.0 # [cm/s]
    vg = 1e-5

    # number density
#    N = 1e8
    N = 1e12

    # average mass
    m = 2.3 * mp

    # temperature
    T = 100.0
    vth = np.sqrt(8 * kk * T / np.pi / m)

    # calculate the coefficients based on the analytical model
    coef = model_tensor1.coefficients(N*m, vth, vg, c, a/c, rho_s, abs(gb))

    # some timescale
    R = 10 * au
    wk = np.sqrt(gg * ms / R**3)
    Tk = 2. * np.pi / wk

    T_damp = rho_s * c / N / m / vth
    T_osc = np.sqrt(rho_s * c**2 / N / m / vth / vg)

    # ==== settings ====
    y0 = [5 * rad, 0]
    t = np.linspace(0, 2e2, 50) * T_osc
#    t = np.linspace(0, 1e-4, 100 ) * T_damp
#    t = np.linspace(0, 1, 100) * 1e6 # of order 1e5

    # ==== solve grid-based ====
    args = (a, c, coef['I'], gb, vg, N, m, T)
    sol = odeint(fn, y0, t, args=args)
    sol_g = {'t':t, 'theta':sol[:,0], 'omega':sol[:,1]}

    # ==== tensor analytical model ====
    t = np.linspace(t[0], t[-1], 300)
    args = (coef['I'], coef['D'], coef['K'])
    sol = odeint(fn2, y0, t, args=args)
    sol_t = {'t':t, 'theta':sol[:,0], 'omega':sol[:,1]}

    # ==== plotting ====
    plot1([sol_g, sol_t], ['Grid-based', 'Analytical'], ['ko', 'C0-'])

#    plot2(t/T_damp, sol_g[:,0], sol_g[:,1])

    pdb.set_trace()

if __name__ == '__main__':
    main()
    
