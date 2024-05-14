# prob3.py
# take a look at the timscale
# use the normalized equation
import numpy as np
import matplotlib.pyplot as plt
import pdb
import natconst
au = natconst.au
year = natconst.year
rad = np.pi / 180

def plot1(r, t_orbit, t_stable, t_osc, t_stop):
    ax = plt.gca()

    ax.plot(r / au, t_orbit/year, label='orbit')
    ax.plot(r / au, t_stable/year, label='stable')
    ax.plot(r / au, t_osc/year, label='oscillation')
#    ax.plot(r / au, t_stop/year, label='stop')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('r [au]')
    ax.set_ylabel('t [year]')

    ax.legend()

    plt.show()

def main():
    # ==== settings ====
    ms = 1 * natconst.ms # mass of star
    r = np.geomspace(1, 150, 20) * au

    # temperature
    q = -0.5
    T = 20 * (r / 100 / au)**(q)

    # gas surface density
    p = -1
    sigma_g = 1. * (r / 100 / au)**(p)

    # ==== details of grain ====
    f = 0.01 # difference in center of mass
    g = 0 # difference in size

    # length between the two spheres of the grain
    l = 0.1 # [cm]

    # amount of empty space
    s = 0.1

    # average density 
    rho_avg = 3.0

    # center of mass distance
    la = (1 - f) / 2 * l
    lb = (1 + f) / 2 * l

    # radii
    ra = (1 - s - g) / 2 * l
    rb = (1 - s + g) / 2 * l

    # total volume
    vol_a = 4 * np.pi / 3 * ra**3
    vol_b = 4 * np.pi / 3 * rb**3
    vol = vol_a + vol_b

    # determine total mass
    m = vol * rho_avg

    # mass of each sphere
    ma = (1 + f) / 2 * m
    mb = (1 - f) / 2 * m

    rho_a = ma / vol_a
    rho_b = mb / vol_b

    # ==== calculate profiles ====
    cs = np.sqrt(natconst.kk * T / natconst.mp / 2.3)
    wk = np.sqrt(natconst.gg * ms / r**3)
    Hg = cs / wk
    rhog = sigma_g / np.sqrt(2*np.pi) / Hg

    # thermal velocity of the molecules
    vth = np.sqrt(8 * natconst.kk * T / np.pi / natconst.mp / 2.3)

    # gas rotation frequency
    w_g = wk * (1 + 0.5 * (Hg / r)**2 * (p + q))

    eta = - (Hg / r)**2 * (p + q)

    Q = wk * cs / np.pi / natconst.gg / sigma_g # toomre Q

    # stopping time
    t_stop = rho_avg / rhog * l / vth

    # Stokes number
#    St = np.pi / 2 * a * rho_s / sigma_g
    St = t_stop * wk

    # velocity difference
    dvphi = - St**2 / (1 + St**2) * eta * wk * r 

    # calculate coefficients
    D = 4 * np.pi / 3 * rhog * vth * (ra**2 * la**2 + rb**2 * lb**2) / (ma*la**2 + mb*lb**2)

    K = 4 * np.pi / 3 * rhog * vth * abs(dvphi) * (-ra**2 * la + rb**2 * lb) / (ma * la**2 + mb*lb**2)

    # stabilizing time
    t_stable = 1. / D

    # oscillation time
    t_osc = 2. * np.pi / np.sqrt(K)

    # orbital time
    t_orbit = 2 * np.pi / wk

    # ==== plotting ====
    plot1(r, t_orbit, t_stable, t_osc, t_stop)

    pdb.set_trace()

if __name__ == '__main__':
    main()

