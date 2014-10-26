""" This module contains the base gridder class."""
from kernels import kaiser_bessel2
from trajectories import radial_2d
from grid_methods import simple_igrid2d
from utils import *
import numpy as np
import pylab as pl

class GridParams:
    """ This is class to hold all of the important gridding parameters."""

    def __init__(self):
        pass

class GridMethod:
    """This is a parent class for all of the individual gridding methods"""

    def __init__(self):
        pass

class Gridder:

    """This is the main gridding class."""

    def __init__(self):
        """Init Gridder with default testing values."""
        self.krad = 1.5
        self.grid_mod = 100
        self.over_samp = 1.5

        self.kernx, self.kerny, fx, fy = kaiser_bessel2(
            self.krad*2.0, self.grid_mod, self.over_samp)

    def igrid_2d(self, kspace, trajectory):
        self.kshape = kspace.shape

        data = np.zeros(trajectory[0].shape, 'c16')

        for i in range(data.size):
            x = (trajectory[0].flat[i] + .5) * kspace.shape[0]
            y = (trajectory[1].flat[i] + .5) * kspace.shape[1]
            xmin = int(np.floor(x - self.krad))
            ymin = int(np.floor(y - self.krad))
            xmax = int(np.ceil(x + self.krad))
            ymax = int(np.ceil(y + self.krad))
            for iy in range(ymin, ymax + 1):
                if iy >= 0 and iy < kspace.shape[1]:
                    dy = y - iy
                    for ix in range(xmin, xmax + 1):
                        if ix >= 0 and ix < kspace.shape[0]:
                            dx = x - ix
                            dr = np.sqrt(dx * dx + dy * dy)
                            kval = self.get_kval(dr)
                            data.flat[i] += kval * kspace[ix, iy]
        return data

    def get_kval(self, dr):
        if dr >= self.krad:
            return 0.0
        else:
            i = 0
            while self.kernx[i] <= dr:
                i += 1
            y = self.kerny[i - 1] + (self.kerny[i] - self.kerny[i - 1]) * \
                (dr - self.kernx[i - 1]) / (self.kernx[i] - self.kernx[i - 1])
            return y

    def grid_2d(self, data, trajectory, density=0):
        kspace = np.zeros(self.kshape, 'c16')

        for i in range(data.size):
            x = (trajectory[0].flat[i] + .5) * kspace.shape[0]
            y = (trajectory[1].flat[i] + .5) * kspace.shape[1]
            xmin = int(np.floor(x - self.krad))
            ymin = int(np.floor(y - self.krad))
            xmax = int(np.ceil(x + self.krad))
            ymax = int(np.ceil(y + self.krad))
            for iy in range(ymin, ymax + 1):
                if iy >= 0 and iy < kspace.shape[1]:
                    dy = y - iy
                    for ix in range(xmin, xmax + 1):
                        if ix >= 0 and ix < kspace.shape[0]:
                            dx = x - ix
                            dr = np.sqrt(dx * dx + dy * dy)
                            kval = self.get_kval(dr)
                            kspace[ix, iy] += data.flat[i] * density.flat[i] * kval
        return kspace


if __name__ == "__main__":

    x = np.linspace(-.5, .5, 128)
    y = np.linspace(-.5, .5, 128)
    xx, yy = np.meshgrid(x, y)
    rr = np.sqrt(xx * xx + yy * yy)
    im = np.zeros((128, 128))
    im[rr < .4] = 2.0
    im[rr < .3] = 1.0

    im = zeropad_ratio(im, 1.5)

    kspace = np.fft.fftshift(im)
    kspace = np.fft.fft2(kspace)
    kspace = np.fft.ifftshift(kspace)

    (traj, dens) = radial_2d(160, 161)

    gridder = Gridder()
    data = gridder.igrid_2d(kspace, traj)
    k1 = gridder.grid_2d(data, traj, dens)

    im1 = np.fft.fftshift(k1)
    im1 = np.fft.ifft2(im1)
    im1 = np.fft.ifftshift(im1)

    # pl.plot(abs(data[0]))
    # pl.show()

    # pl.imshow(np.log(abs(k1)))
    # pl.imshow(abs(im1))
    # pl.show()

    demod = abs(im)/abs(im1)
    pl.plot(demod[:,96])
    pl.show()
