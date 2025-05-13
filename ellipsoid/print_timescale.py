"""
adopt a few parameter and see the timescale

Also take a look at the threshold flow speed
"""
import numpy as np
import matplotlib.pyplot as plt
import natconst
import model_tensor1
import pdb
alphabet = 'abcdefghijklmn'
lincol = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']

def get_hat_W(abar, gbar):
    """
    calculate the geometry factor
    """
    E = model_tensor1.get_H(abar)

    hat_W = 1. / (16 * abar * gbar * E)
    return hat_W

def get_results(c, rho_s, gbar, abar, rho_g, W, vth):
    # calculate coefficients ====
    coef = model_tensor1.coefficients(rho_g, vth, W, c, abar, rho_s, gbar)
    time = model_tensor1.timescales(coef['I'], coef['D'], coef['K'])
    other = model_tensor1.other_timescales(rho_g, vth, W, c, rho_s)

    print('stopping time = %e [year]'%(other['t_stop'] / natconst.year))
    print('damping time = %e [years]'%(time['t_d'] / natconst.year))
    print('characteristic oscillation time = %e [hours]'%(other['t_osc_c'] / natconst.hour))
    print('oscillation time = %e [hours]'%(time['t_o'] / natconst.hour))

    # calculate the threshold flow speed
    n = rho_g / 2.3 / natconst.mp
    hat_W = get_hat_W(abar, gbar)
    W_t = vth / n / c**3 * hat_W
    print('threshold flow = %e [cm/s]'%(W_t))

def main():
    """
    ISM conditions
    """
    # ==== settings ====
    c = 1e-5 # [cm]
    rho_s = 3.0
    gbar = 0.01
    abar = 0.9

    # environmental inputs
    n = 20.0
    rho_g = n * 2.3 * natconst.mp

    W = 1 # cm/s]

    # vth part
    T = 10 # [kelvin]
#    vth = np.sqrt(8 * natconst.kk * T / np.pi / 2.3 / natconst.mp)
#    print('vth = %e [cm/s]'%vth)
    vth = 0.3e5

    print('==== ISM conditions ====')
    get_results(c, rho_s, gbar, abar, rho_g, W, vth)
    print(' ')

def main_cores():
    """
    dense molecular cores
    """
    # ==== settings ====
    c = 1e-4 # [cm]
    rho_s = 3.0
    gbar = 0.01
    abar = 0.9

    # environmental inputs
    n = 1e5
    rho_g = n * 2.3 * natconst.mp

    W = 1 # cm/s]

    # vth part
    T = 10 # [kelvin]
#    vth = np.sqrt(8 * natconst.kk * T / np.pi / 2.3 / natconst.mp)
#    print('vth = %e [cm/s]'%vth)
    vth = 0.3e5

    print('==== core conditions ====')
    get_results(c, rho_s, gbar, abar, rho_g, W, vth)
    print(' ')

def main_disk():
    """
    disk conditions
    """
    # ==== settings ====
    c = 1e-1 # [cm]
    rho_s = 3.0
    gbar = 0.01
    abar = 0.9

    # environmental inputs
    n = 1e9
    rho_g = n * 2.3 * natconst.mp

    W = 1 # cm/s]

    # vth part
    T = 100 # [kelvin]
#    vth = np.sqrt(8 * natconst.kk * T / np.pi / 2.3 / natconst.mp)
#    print('vth = %e [cm/s]'%vth)
    vth = 1e5

    print('==== disk conditions ====')
    get_results(c, rho_s, gbar, abar, rho_g, W, vth)


if __name__ == '__main__':
    main()
    main_cores()
    main_disk()

