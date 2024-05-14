"""
prob_integ_const.py

Take a look at the integration constants

Compare three different methods of calculation.
1. integrate with a finite difference scheme using legendre quadrature
2. integrate using the scipy.integrate.quadrature method
3. analytical solution
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quadrature as quad
import pdb
import time

def H_fn(x, abar=1):
    """
    the function 
    """
    return (1 - x**2) / np.sqrt(1 + (abar**2 - 1) * x**2)

def H_quad(abar):
    """
    use integratation
    """
    res = quad(H_fn, -1, 1, args=abar)

    return res[0]

def H_integ(abar, nx=21):
    """
    integrate yourself
    """
    x, dx = np.polynomial.legendre.leggauss(nx)

    y = H_fn(x, abar=abar)

    H = np.sum(y * dx)

    return H

def H_ana(abar):
    """
    analytical solution to the integral
    For abar = 1, the solution is 4/3
    """
    if abar == 1:
        H = 4/3.
    else:
        H = (2*abar**2 - 1) * np.arcsinh(np.sqrt(abar**2-1)) / (abar**2-1)**1.5 - abar / (abar**2 - 1)

    return H

# =============================================
def L_fn(x, abar=1):
    """
    the function
    """
    return (1 - x**2) * x**2 / np.sqrt(1 + (abar**2 - 1) * x**2)

def L_quad(abar):
    """
    use integratation
    """
    res = quad(L_fn, -1, 1, args=abar)

    return res[0]

def L_integ(abar, nx=21):
    """
    integrate yourself
    """
    x, dx = np.polynomial.legendre.leggauss(nx)

    y = L_fn(x, abar=abar)

    L = np.sum(y * dx)

    return L

def L_ana(abar):
    """
    analytical solution to the integral
    For abar = 1, the solution is 4/3
    """
    if abar == 1:
        L = 4. / 15
    else:
        L = - (np.sqrt(abar**2 - 1) * (4*abar**2-1) * np.arcsinh(np.sqrt(abar**2-1)) - 2 * abar**5 + abar**3 + abar )/ 4 / (abar**6 - 3*abar**4 + 3*abar**2 - 1)

    return L

# =============================================
def W_fn(x, abar=1):
    """
    the function
    """
    return x**2 / np.sqrt(1 + (abar**2 - 1) * x**2)

def W_quad(abar):
    """
    use integratation
    """
    res = quad(W_fn, -1, 1, args=abar)

    return res[0]

def W_integ(abar, nx=21):
    """
    integrate yourself
    """
    x, dx = np.polynomial.legendre.leggauss(nx)

    y = W_fn(x, abar=abar)

    res = np.sum(y * dx)

    return res

def W_ana(abar):
    """
    analytical solution to the integral
    For abar = 1, the solution is 4/3
    """
    if 0 < abar < 1:
        e = np.sqrt(1 - abar**2)
#        res = 1./e**3 * np.log((1+e)/(1-e)) - 2/e**2
        res = (np.arcsin(e) - np.arcsin(-e)) / 2. / e**3
    elif abar == 1:
        res = 2. / 3
    elif abar > 1:
        f = np.sqrt(abar**2 - 1)
        res = np.sqrt(1+f**2)/f**2 - np.arcsinh(f) / f**3
    else:
        raise ValueError('abar must be a float')

    return res

# =============================================
def main():
    """
    """
    abar = np.linspace(0.5, 1.5, 5)

    const = 'W'

    # ==== integration ====
    fn = globals()['%s_integ'%const]
    tick = time.time()
    res_1 = np.zeros_like(abar)
    for i in range(len(abar)):
        res_1[i] = fn(abar[i])
    tock = time.time()
    time_integ = tock - tick

    # ==== scipy quadrature ====
    fn = globals()['%s_quad'%const]
    tick = time.time()
    res_2 = np.zeros_like(abar)
    for i in range(len(abar)):
        res_2[i] = fn(abar[i])
    tock = time.time()
    time_quad = tock - tick

    # ==== analytical ====
    fn = globals()['%s_ana'%const]
    tick = time.time()
    res_3 = np.zeros_like(abar)
    for i in range(len(abar)):
        res_3[i] = fn(abar[i])
    tock = time.time()
    time_ana = tock - tick

    print("Time for integration = %f"%time_integ)
    print("Time for scipy quadratic = %f"%time_quad)
    print("Time for analytical = %f"%time_ana)

    # ==== plotting ====
    plt.plot(abar, res_1, 'o', label='%s integration'%const)
    plt.plot(abar, res_2, '+', label='%s quad'%const)
    plt.plot(abar, res_3, label='%s analytical'%const)
    plt.xlabel(r'$\bar{a}$')
    plt.ylabel(r'$%s$'%const)
    plt.legend()
    plt.show()

if  __name__ == '__main__':
#    main_H()
#    main_L()
#    main_W()
    main()

