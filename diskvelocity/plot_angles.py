"""
plot images of polarization angles given known drift directions
"""
import pdb
import numpy as np
import matplotlib.pyplot as plt

def sample_points(fn_r, rin=10, rout=100, unitlen=5):
    """
    sample points for x and y
    
    """
    # sample in radius
    nr = 5
    r = np.linspace(rin, rout, nr)

    # sample in phi
    nphi = 8
    phi = np.linspace(0, 2*np.pi, nphi, endpoint=False)

    # determine the x, y points
    rr, pp = np.meshgrid(r, phi, indexing='ij')
    xx = rr * np.cos(pp)
    yy = rr * np.sin(pp)

    # determine the vector directions
    ang = fn_r(rr)

    # determine the directions in cartesian coordinates
    vx = unitlen * np.cos(ang)
    vy = unitlen * np.sin(ang)

    return xx.flatten(), yy.flatten(), vx, vy

def add_plot(ax, fn_r):
    """
    add vectors
    fn_r : interp object
        this should give the angles as a function of radius in au
        the angle should be in radians
    """
    # plot the outer rim of the disk
    phi = np.linspace(0, 2*np.pi, 361)
    rout_x = rout * np.cos(phi)
    rout_y = rout * np.sin(phi)
    ax.plot(rout_x, rout_y, 'k')

    # plot a center point for the star
    ax.plot([0], [0], 'ko')

    # sample the points we want
    x, y, vx, vy = sample_points(fn_r)


