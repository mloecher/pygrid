import numpy as np
import pylab as pl


def radial_2d(nproj, res):
    kx = np.zeros((nproj, res))
    ky = np.zeros((nproj, res))

    angles = np.linspace(0, np.pi, nproj + 1)
    angles = angles[:-1]

    rad = np.linspace(-.5, .5, res)
    for i in range(nproj):
        kx[i, :] = rad * np.cos(angles[i])
        ky[i, :] = rad * np.sin(angles[i])

    kxx = kx.flatten()
    kyy = ky.flatten()

    # pl.scatter(kxx, kyy)
    # pl.show()

    trajectory = np.zeros((2,) + kx.shape)
    trajectory[0, :] = kx
    trajectory[1, :] = ky

    density = np.zeros(kx.shape)
    density = np.sqrt(kx*kx + ky*ky)
    # density = 1/density

    return (trajectory, density)

if __name__ == "__main__":
    radial_2d(10, 19)
