""" This module contains the base gridder class."""
from kernels import kaiser_bessel2
from trajectories import radial_2d
from grid_methods import simple_igrid2d
from grid_kernel import GridKernel
from python_forloops import *
from python_vectors import *
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
    krad = 2.5
    grid_mod = 16
    over_samp = 2.5

    # Image defaults:
    imsize = (128, 128)
    imsize_os = tuple([x * over_samp for x in imsize])


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

        # self.grid_method = GMPythonForLoops()
        self.grid_method = GMPythonVectors()
        self.grid_method.setup()

    def igrid_FT(self, im):
        im = zeropad_ratio(im, self.grid_params.over_samp)

        # pl.imshow(abs(im))
        # pl.show()

        im = im / self.kernel.deapp

        kspace = np.fft.fftshift(im)
        kspace = np.fft.fft2(kspace)
        kspace = np.fft.ifftshift(kspace)

        # pl.imshow(np.log(abs(kspace)))
        # pl.show()

        return self.grid_method.igrid_2d(kspace, self.grid_params, self.kernel, self.traj)

    def grid_FT(self, data):
        kspace = self.grid_method.grid_2d(
            data, self.grid_params, self.kernel, self.traj, self.dens)

        # pl.imshow(np.log(abs(kspace)))
        # pl.show()

        im = np.fft.fftshift(kspace)
        im = np.fft.ifft2(im)
        im = np.fft.ifftshift(im)

        im = im / self.kernel.deapp

        im = crop_ratio(im, self.grid_params.over_samp)

        return im


if __name__ == "__main__":

    x = np.linspace(-.5, .5, 128)
    y = np.linspace(-.5, .5, 128)
    xx, yy = np.meshgrid(x, y)
    rr = np.sqrt(xx * xx + yy * yy)
    im = np.zeros((128, 128))
    im[rr < .4] = 1.0
    im[rr < .3] = .8
    im[rr < .2] = .6

    # im = np.zeros((128, 128))
    # im[10:118, 10:118] = 1.0

    (traj, dens) = radial_2d(24, 151)

    print traj.shape

    gridder = Gridder(traj, dens)

    data = gridder.igrid_FT(im)

    # pl.plot(abs(data[0]))
    # pl.plot(abs(data[10]))
    # pl.plot(dens[0])
    # pl.show()
    #
    #
    # print dens
    # print dens.shape
    # print abs(data[0]).max()

    im1 = gridder.grid_FT(data)

    diff = abs(im) / abs(im1)
    center = diff.shape[0] / 2

    pl.figure()
    pl.plot(diff[center, :])

    pl.figure()
    # pl.imshow(abs(im)/abs(im1), cmap=pl.gray())
    pl.imshow(np.abs(im1), cmap=pl.gray())
    pl.colorbar()
    pl.show()
