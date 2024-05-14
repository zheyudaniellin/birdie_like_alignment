"""
directions of the velocity for some given pressure profile

I want to calculate the gas drift field
"""
import numpy as np
import matplotlib.pyplot as plt
from natconst import gg, ms, au
import pdb
from diskvel_alpha import rel_v_r, rel_v_phi, get_angle
rad = np.pi / 180

# ==== plotting ====
def plot1(grid, press, dvr, dvphi, ang):
    """
    plot the profiles
    """
    nrow, ncol = 4, 1
    fig, axgrid = plt.subplots(nrow,ncol,sharex=True,sharey=False,
            squeeze=False)
    axes = axgrid.flatten()

    r = grid['r']

    # pressure
    ax = axes[0]
    ax.plot(r, press)
    ax.set_yscale('log')
    ax.set_ylabel('Pressure []')

    # dvr
    ax = axes[1]
    ax.plot(r, dvr)
    ax.set_ylabel('dvr []')

    # dvphi
    ax = axes[2]
    ax.plot(r, dvphi)
    ax.set_ylabel('dvphi')

    # angle
    ax = axes[3]
    ax.plot(r, ang/rad)
    ax.set_ylabel('angle [deg]')

    ax.set_xscale('log')

    fig.tight_layout()
    plt.show()
    
# ==== models of ln P ====
"""
keep different models of the pressure gradient 
basically we need to know ln P as a function of ln r
Note that P = rho * cs^2

since we only want to know d ln P / d ln r, we don't need to know the constants. 
ln P = ln rho + ln T + constants

Thus, 
d ln P / d ln r = d ln rho / d ln r + d ln T / d ln r
"""
def model1(r, p=1, q=0.5):
    """
    just calculate the pressure profile
    """
    # density
    rho= 1 * r**(-p)
    # temperature
    T = 1 * r**(-q)

    return rho * T

def model2(r, ringpar=None,):
    """
    disk with rings with constant temperature

    ringpar : 2d ndarray
        each ringpar[i,:] is a list of [r, w]
        where r is the radius of the ring center and w is the ring width. 
        the units should be the same as the input grid
        
    """
    rloc = ringpar[:,0]
    width = ringpar[:,1]

    # density on the walls
    rho = np.sum(np.exp(-0.5*((r[:,None] - rloc[None,:]) / width[None,:])**2), axis=1)

    # temperature
    T = 1

    return rho * T

def get_pressure_gradient(grid, fn, fnkwargs=None):
    """
    calculate the pressure
    """
    # ln_r at the walls
    ln_r_w = np.log(grid['r_w'])
    d_ln_r = ln_r_w[1:] - ln_r_w[:-1]

    # pressure at the walls
    ln_P_w = np.log(fn(grid['r_w'], **fnkwargs))

    d_ln_P = ln_P_w[1:] - ln_P_w[:-1]

    return d_ln_P / d_ln_r

# ==== main pipeline ====
def main():
    """
    some demonstration
    """
    # ==== settings ====
    alpha = 0
    St = 10.

    # determine the walls
    nr = 61
    r_w = np.geomspace(1, 100, nr+1) # in au
    r = (r_w[1:] * r_w[:-1])**0.5
#    r_w = np.linspace(20, 80, nr+1)
#    r = 0.5 * (r_w[1:] + r_w[:-1])

    grid = {'nr':nr, 'r_w':r_w, 'r':r}

    # aspect ratio
    h = 10 * (r / 100)**1.25
    hr = h / r
#    hr = 0.1

    # keplerian velocity
    vk = np.sqrt(gg * 1*ms / (r* au))

    # ==== determine which model ====
    fn, fnkwargs = model1, {'p':1, 'q':0.5}
#    fn, fnkwargs = model2, {'ringpar':np.array([
#        [50, 5]
#        ])}

    # ==== calculate pressure profile ====
    press = fn(r, **fnkwargs)
    beta = get_pressure_gradient(grid, fn, fnkwargs=fnkwargs)

    # ==== calculate the gas drift ====
    dvr = rel_v_r(St, alpha, beta, hr, vk)
    dvphi = rel_v_phi(St, alpha, beta, hr, vk)

    # angle
    ang = get_angle(dvr, dvphi)

    # ==== plotting ====
    plot1(grid, press, dvr, dvphi, ang)

    pdb.set_trace()

if __name__ == '__main__':
    main()
