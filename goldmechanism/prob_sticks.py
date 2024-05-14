"""
just overplot transparent sticks to show that most of the projection will be along the gas drift

Assume the gas drift is in the plane
"""
import numpy as np
import matplotlib.pyplot as plt

def sdir(thet, phi):
    """
    vector onto the surface of the unit sphere
    """
    return np.array([
        np.sin(thet) * np.cos(phi),
        np.sin(thet) * np.sin(phi),
        np.cos(thet)
        ])

def get_stick(thet, phi, l=1):
    """
    calculate the coordinates of the stick using two points
    """
    vec = sdir(thet, phi)

    x = np.array([-1, 1]) * vec[0] * l/2
    y = np.array([-1, 1]) * vec[1] * l/2
    z = np.array([-1, 1]) * vec[2] * l/2
    return x, y, z

def main():
    # ==== settings ====
    npoints = 1000
    thet = np.random.rand(npoints) * np.pi
    phi = np.random.rand(npoints) * 2*np.pi

    # ==== make plot ====
    ax = plt.figure().add_subplot(projection='3d')
    for i in range(npoints):
        x, y, z = get_stick(thet[i], phi[i])
        ax.plot(x, y, z, color='k', alpha=0.025, linewidth=5)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    plt.show()

if __name__ == '__main__':
    main()

