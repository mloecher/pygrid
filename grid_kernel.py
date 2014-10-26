import pylab as pl
import numpy as np

class GridKernel():

    def __init__(self, grid_params):
        self.calc_kernel_kb(grid_params)
        self.fourier_demod(grid_params)
        self.krad = grid_params.krad

    def get_kval(self, dr):
        if dr >= self.krad:
            return 0.0
        else:
            i = 0
            while self.kx[i] <= dr:
                i += 1
            y = self.ky[i - 1] + (self.ky[i] - self.ky[i - 1]) * \
                (dr - self.kx[i - 1]) / (self.kx[i] - self.kx[i - 1])
            return y

    def calc_kernel_kb(self, grid_params):

        kw0 = 2.0*grid_params.krad/grid_params.over_samp
        kr = grid_params.krad

        beta = np.pi*np.sqrt((kw0*(grid_params.over_samp-0.5))**2-0.8)

        x = np.linspace(0, kr, grid_params.grid_mod)

        x_bess = np.sqrt(1 - (x / kr) ** 2)

        y = np.i0(beta * x_bess)
        y = y / y[0]

        self.kx = x
        self.ky = y

    def fourier_demod(self, grid_params):
        xres = grid_params.imsize[0] * grid_params.over_samp
        Dx = np.arange(xres)
        Dx = Dx - xres/2.0
        Dy = np.zeros(Dx.size, 'c16')

        for i in range(1, self.kx.size):
            temp = self.ky[i] * 2*np.exp(2*1j*np.pi*Dx/xres*self.kx[i])
            Dy += temp

        Dy = Dy.real
        Dy = Dy - Dy.min()
        Dy = Dy/Dy.max()

    def calc_deapp(self, grid_params):
