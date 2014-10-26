""" This module contains the base gridder class."""
from kernels import kaiser_bessel2
from trajectories import radial_2d
from grid_methods import simple_igrid2d
from grid_kernel import GridKernel
from python_forloops import *
from utils import *
import sys
import numpy as np
import pylab as pl

class GridParams:
    """ This is class to hold all of the important gridding parameters.

    It could probably just be a dictionary, but maybe we willl need some
    functions as we develop?
    """

    # Kernel defaults:
    kernel_type = 'kaiser_bessel'
    krad = 1.5
    grid_mod = 100
    over_samp = 1.25

    # Image defaults:
    imsize = (128, 128)


class GridMethod:
    """This is a parent class for all of the individual gridding methods"""

    def __init__(self):
        pass

    def setup(self):
        pass

class Gridder:

    """This is the main gridding class."""

    def __init__(self, trajectory, density):
        self.grid_params = GridParams()
        self.kernel = GridKernel(self.grid_params)
        self.traj = trajectory
        self.dens = density

        self.grid_method = GMPythonForLoops()
        self.grid_method.setup()

    def igrid_FT(self, im):
        im = zeropad_ratio(im, self.grid_params.over_samp)

        kspace = np.fft.fftshift(im)
        kspace = np.fft.fft2(kspace)
        kspace = np.fft.ifftshift(kspace)

        return self.grid_method.igrid_2d(kspace, self.grid_params, self.kernel, self.traj)

    def grid_FT(self, data):
        kspace = self.grid_method.grid_2d(data, self.grid_params, self.kernel, self.traj, self.dens)

        im = np.fft.fftshift(kspace)
        im = np.fft.ifft2(im)
        im = np.fft.ifftshift(im)

        im = crop_ratio(im, self.grid_params.over_samp)

        return im


if __name__ == "__main__":

    x = np.linspace(-.5, .5, 128)
    y = np.linspace(-.5, .5, 128)
    xx, yy = np.meshgrid(x, y)
    rr = np.sqrt(xx * xx + yy * yy)
    im = np.zeros((128, 128))
    im[rr < .4] = 2.0
    im[rr < .3] = 1.0

    (traj, dens) = radial_2d(20, 161)

    gridder = Gridder(traj, dens)

    data = gridder.igrid_FT(im)
    im1 = gridder.grid_FT(data)

    # pl.plot(abs(data[0]))
    pl.imshow(abs(im)/abs(im1))
    pl.show()

    sys.exit()

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
