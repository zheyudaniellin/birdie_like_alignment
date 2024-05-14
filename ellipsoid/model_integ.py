"""
prob1.py
calculate the pressure at a given location
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import dblquad
import pdb
rad = np.pi / 180
kk = 1.3807e-16  # Bolzmann's constant [erg/K]
mp = 1.6726e-24  # Mass of proton [g]

def body_to_lab(xb,yb,zb, theta, gb):
    x = np.cos(theta) * xb + np.sin(theta) * (zb - gb)
    y = yb * 1
    z = - np.sin(theta) * xb + np.cos(theta) * (zb - gb)
    return x, y, z

def integrand(w, phi, a, c, gb, theta, omega, vg, N, m, T):
    """
    calculate the integral
    """

    # ==== calculations ====
    # at a particular location in the body frame
    zb = w * c
    Rw = a * np.sqrt(1. - w**2)
    xb = Rw * np.cos(phi)
    yb = Rw * np.sin(phi)

    # calculate the lab frame
    x, y, z = body_to_lab(xb,yb,zb,theta,gb)

    # calculate velocity vector
    V = np.array([omega * z, 0, -(omega*x + vg)])
    Vmag = np.sqrt(np.sum(V**2))

    # calculate body unit vectors
    xb_hat = np.array([np.cos(theta), 0, -np.sin(theta)])
    yb_hat = np.array([0, 1, 0])
    zb_hat = np.array([np.sin(theta), 0, np.cos(theta)])
    phi_hat = - np.sin(phi) * xb_hat + np.cos(phi) * yb_hat

    # surface frame unit vectors
#    norm = 2./a**2 * (xb*xb_hat + yb*yb_hat) + 2./c**2 * zb*zb_hat
#    xs_hat = norm / np.sqrt(np.sum(norm**2))
    a_bar = a / c
    xs_hat = (xb * xb_hat + yb*yb_hat + a_bar**2 * zb * zb_hat) / a / np.sqrt(1 + (a_bar**2 - 1) * w**2)

    ys_hat = phi_hat
    zs_hat = np.cross(xs_hat, ys_hat)

    # components of the velocity vector in surface frame
    alpha = np.dot(V, xs_hat) / Vmag
    beta = np.dot(V, ys_hat) / Vmag
    gamma = np.dot(V, zs_hat) / Vmag

    # ==== calculate the forces ====
    cs = np.sqrt(kk * T / m)

    # in the x-direction
    # calculate the fs without the -2 * N * m * cs
    alpha_p = xs_hat[0]
    fs_x_part = alpha_p * (1 + np.sqrt(8/np.pi) * Vmag/cs * alpha)

    # in the z-direction
    # calculate the fs without the -2 * N * m * cs
    alpha_p = xs_hat[2]
    fs_z_part = alpha_p * (1 + np.sqrt(8/np.pi) * Vmag/cs * alpha)

    # ==== calculate the integrand ====
    spart = np.sqrt(1 + ((a/c)**2 - 1)*w**2)

    out = spart * (z * fs_x_part - x * fs_z_part)

    return out

def do_dblquad(a, c, gb, theta, omega, vg, N, m, T):
    """
    use dblquad
    """
    # some physical quantities
    cs = np.sqrt(kk * T / m)

    args = (a, c, gb, theta, omega, vg, N, m, T)
    res = dblquad(integrand, 0, 2*np.pi, -1,1, args=args)
    torque = res[0] * a * c * (-2 * N * m * cs**2) 

    return torque

def do_grid(a, c, gb, theta, omega, vg, N, m, T, nphi=50, nw=40):
    """
    instead of calling dblquad, use a grid for integration
    """
    # some physical quantities
    cs = np.sqrt(kk * T / m)

    # grid walls for phi
    phi = np.linspace(0, 2*np.pi, nphi, endpoint=False)
    dphi = phi[1] - phi[0]

    # grid walls for w
    w = np.linspace(-1, 1, nw)
    w_s = 0.5 * (w[1:] + w[:-1])
    dw = np.diff(w)

    sdwdphi = np.zeros([nw-1, nphi])
    for i in range(nw-1):
        for j in range(nphi):
            sdwdphi[i,j] = integrand(w_s[i], phi[j], a, c, gb, theta, omega, vg, N, m, T)

    # integrate over phi
    sdw = np.sum(sdwdphi * dphi, axis=1)

    # integrate over w
    s = np.sum(sdw * dw, axis=0)

    torque = s * a * c * (-2 * N * m * cs**2)

    return torque

def main():
    # ==== settings ====
    # prolate shape
    c = 0.1 # [cm]
    a = 0.9 * c

    # location of the center of mass in body frame
    gb = -0.1 * c

    # angle of prolate
    theta = 90 * rad

    # angular velocity
    omega = 0

    # velocity of gas
    vg = 10.0 # [cm/s]

    # number density
    N = 1e8

    # average mass
    m = 2.3 * mp

    # temperature
    T = 20.0

    # calculate torques
    torque_quad = do_dblquad(a, c, gb, theta, omega, vg, N, m, T)
    torque_grid = do_grid(a, c, gb, theta, omega, vg, N, m, T)
#    print(torque)

    cs = np.sqrt(kk * T / m)
    ana = N * m * cs**2 * 4*np.pi*a*c * a
    print(torque_quad / ana)
    print(torque_grid / ana)

    pdb.set_trace()

if __name__ == '__main__':
    main()



