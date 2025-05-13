"""
calculate the threshold for W

"""
import numpy as np
import matplotlib.pyplot as plt
import natconst
import model_tensor1
import pdb

def get_hat_W(abar, gbar):
    """
    calculate the geometry factor
    """
    E = model_tensor1.get_H(abar)

    hat_W = 1. / (16 * abar * gbar * E)
    return hat_W

def main():
    # ==== settings ====
    vth = 1e5 # [cm /s]
    n = 20.0 # [cm^-3]
    c = 1e-5 # [cm]

    abar = 0.9
    gbar = 0.01

    hat_W = get_hat_W(abar, gbar)

    W_t = vth / n / c**3 * hat_W

    print('hat_W = %e'%hat_W)
    print('W_t = %e [cm/s]'%W_t)

if __name__ == '__main__':
    main()

