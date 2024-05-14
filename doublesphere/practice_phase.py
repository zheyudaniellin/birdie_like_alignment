"""
practice using the phaseportrait code
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import phaseportrait

rad = np.pi / 180

def pend(phi, omega, *, b=1, c=1):
    """
    The system of equations:
    phi'' + b * phi' + c sin(phi) = 0
    The second order differential can be written as two first order ODE:
    phi' = omega
    omega' = -b * omega - c * sin(phi)
    """
    return [omega, -b * omega - c * np.sin(phi)]


def main():
    # ==== settings ====
    phi_lim = [-2*np.pi, 2*np.pi]
    omega_lim = [-np.pi, np.pi]
    pp = phaseportrait.PhasePortrait2D(pend, [phi_lim, omega_lim], 
        density=2, 
        Title='Damped Pendulum', xlabel=r"$\phi$", ylabel=r"$\dot{\phi}$")
    pp.add_slider('b', valinit=1)
    pp.add_slider('c', valinit=1)
    pp.plot()

    plt.show()

if __name__ == '__main__':
    main()

