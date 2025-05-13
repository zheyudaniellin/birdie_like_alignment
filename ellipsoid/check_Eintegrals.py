"""
check the E-integrals against the analytical solutions
"""
import numpy as np
import matplotlib.pyplot as plt
import model_tensor1

def main():
    """
    """
    # abar we want to calculate
    abar = np.concatenate((np.geomspace(0.01, 1, 20, endpoint=False), np.geomspace(1, 1e2, 20)))

    # calculate the constants
    ana = model_tensor1.get_H(abar, mode='analytical')
    num = model_tensor1.get_H(abar, mode='numerical')

    ax = plt.gca()
    ax.plot(abar, ana, 'C0', label='analytical')
    ax.plot(abar, num, 'kx', label='numerical')
    ax.legend()
    ax.set_title(r'$E[1-x^{2}]$')
    ax.set_xlabel(r'$\bar{a}$')
    ax.set_xscale('log')
    ax.axvline(1, color='k')
    plt.show()

def main2():
    """
    """
    # abar we want to calculate
    abar = np.concatenate((np.geomspace(0.01, 1, 20, endpoint=False), np.geomspace(1, 1e2, 20)))

    # calculate the constants
    ana = model_tensor1.get_L(abar, mode='analytical')
    num = model_tensor1.get_L(abar, mode='numerical')

    ax = plt.gca()
    ax.plot(abar, ana, 'C0', label='analytical')
    ax.plot(abar, num, 'kx', label='numerical')
    ax.legend()
    ax.set_title(r'$E[x^{2}-x^{4}]$')
    ax.set_xlabel(r'$\bar{a}$')
    ax.set_xscale('log')
    ax.axvline(1, color='k')
    plt.show()
if __name__ == '__main__':
    main()
    main2()
