"""
plot the mean free path of gas
"""
import pdb
import numpy as np
import matplotlib.pyplot as plt
import natconst
from natconst import au, year
import diskvel_alpha
#import Qconst as disk
import pwldisk as disk

def get_mfp(n):
    """
    calculate the mean free path

    n : ndarray
        number density in particles per cm^3
    """
    # collision cross section in cm^2
    sig = 2e-15 # H2

    return 1. / n / sig

def plot1(r, mfp):
    """
    """
    fig, ax = plt.subplots()
    pltx = r / au

    ax.plot(pltx, mfp, 'k')

    ax.axhline(y=0.1, color='k', linestyle='--')

    ax.set_xscale('log')
    ax.set_yscale('log')

    xticks = [1, 10, 100]
    ax.set_xticks(xticks, labels=['{:d}'.format(itick) for itick in xticks])

    ax.set_xlim(pltx[0], pltx[-1])

    ax.set_xlabel('R [au]')
    ax.set_ylabel(r'$\lambda_{\text{mfp}}$ [cm]')

    plt.show()

def main():
    # ==== settings ====
    r = np.geomspace(1, 100, 100) * au

    # ==== disk profiles ====
    rho_g = disk.midplane_density(r)
    n_g = rho_g / 2.3 / natconst.mp

    mfp = get_mfp(n_g)

    # ==== plotting ====
    plot1(r, mfp)

if __name__ == "__main__":
    main()

