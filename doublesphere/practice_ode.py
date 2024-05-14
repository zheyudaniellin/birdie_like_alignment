"""
 practice_ode.py
solve ode using the scipy odeint function
https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html

The example gives the right form that we want
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
rad = np.pi / 180

def pend(y, t, b, c):
    """
    The system of equations:
    phi'' + b * phi' + c sin(phi) = 0
    The second order differential can be written as two first order ODE:
    phi' = omega
    omega' = -b * omega - c * sin(phi)
    """
    phi, omega = y

    dydt = [omega, -b * omega -c * np.sin(phi)]

    return dydt

def plot1(t, phi, omega):
    """
    plot
    """
    nrow, ncol = 2, 1
    fig, axgrid = plt.subplots(nrow,ncol,sharex=False,sharey=False, 
        squeeze=False)
    axes = axgrid.flatten()

    # ==== phi ====
    ax = axes[0]
    ax.plot(t, phi/rad)
    ax.set_ylabel(r'$\phi$')
    ax.axhline(y=0, color='k', linestyle=':')
    ax.set_yticks(np.arange(-180, 181, 45))

    ax = axes[1]
    ax.plot(t, omega/rad)
    ax.set_ylabel(r'$\phi$ dot')

    for ax in axes:
        ax.set_xlabel('t')
        ax.set_xlim(0, None)

    fig.tight_layout()
    plt.show()

def plot2(t, phi, omega):
    fig = plt.figure()

    ax = fig.gca()

    ax.plot(phi/rad, omega/rad)

    fig.tight_layout()
    plt.show()

def main():
    # ==== settings ====
    b = 0.1
    c = 2.
    y0 = [np.pi-0.01, -1]
    t = np.linspace(0, 10, 101)

    # ==== solve equations ====
    sol = odeint(pend, y0, t, args=(b, c))
    phi = sol[:,0]
    omega = sol[:,1]

    # ==== plotting ====
    plot1(t, phi, omega)
    plot2(t, phi, omega)

if __name__ == '__main__':
    main()

