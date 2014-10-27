import pylab as pl
import numpy as np

class GridKernel():

    def __init__(self, grid_params):
        self.calc_kernel_tri(grid_params)

        self.fourier_demod(grid_params)
        self.calc_deapp(grid_params)
        self.krad = grid_params.krad


    def get_kval(self, dr):
        if dr >= self.krad:
            return 0.0
        else:
            i = 0
            while self.kx[i] <= dr:
                i += 1
            dx = self.kx[i] - self.kx[i-1]
            ddr = dr - self.kx[i-1]
            # print ddr
            y = self.ky[i - 1] * (1 - ddr/dx) + self.ky[i] * ddr/dx
            # y = self.ky[i - 1] + (self.ky[i] - self.ky[i - 1]) * \
                # (dr - self.kx[i - 1]) / (self.kx[i] - self.kx[i - 1])
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

    def calc_kernel_tri(self, grid_params):

        kr = grid_params.krad

        x = np.linspace(0, kr, grid_params.grid_mod)

        y = 1.0 - x/kr

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
        Dy = Dy + 1e-3

        pl.plot(Dy)
        pl.show()

        self.Dx = Dx
        self.Dy = Dy

    def calc_deapp(self, grid_params):
        test = np.dot(self.Dy[np.newaxis,:].T, self.Dy[np.newaxis,:])

        # pl.imshow(test)
        # pl.colorbar()
        # pl.show()

        demod = self.Dy[self.Dy.size/2:]

        x = np.arange(grid_params.imsize_os[0]) - grid_params.imsize_os[0]/2
        y = np.arange(grid_params.imsize_os[1]) - grid_params.imsize_os[1]/2
        xx, yy = np.meshgrid(x, y)
        rr = np.sqrt(xx * xx + yy * yy)

        rr = np.sqrt(xx * xx + yy * yy)
        rr[rr>grid_params.imsize_os[0]/2-1] = grid_params.imsize_os[0]/2-1

        rr0 = np.floor(rr)
        rr1 = np.ceil(rr)

        y0 = demod[rr0.flatten().astype('int')]
        y1 = demod[rr1.flatten().astype('int')]

        drr = (rr - rr0).flatten()

        out = (1-drr)*y0 + drr*y1

        out = out.reshape(grid_params.imsize_os)

        self.deapp = test
