import pylab as pl
import numpy as np
import cmath


def kaiser_bessel(kwidth, grid_mod, over_samp = 1.375):
    # kernel = np.zeros(krad * grid_mod)
    beta = np.pi * np.sqrt((kwidth * kwidth) / (over_samp * over_samp)
                           * (over_samp - 0.5) * (over_samp - 0.5) - 0.8)

    print beta

    x = np.linspace(0, kwidth, grid_mod)

    x1 = np.sqrt(1 - (x / kwidth) ** 2)

    y = np.i0(beta * x1)
    y = y / y[0]

    # pl.plot(x, y)
    # pl.show()

    fx = np.arange(192)
    fx = fx - 192/2.0
    fy = np.zeros(192, 'c16')

    deapp = np.sqrt(((np.pi * fx / kwidth)**2 - beta**2).astype('c16'))

    fy = np.sin(deapp)/deapp

    return(x, y, fy)

def kaiser_bessel2(kwidth, grid_mod, over_samp, imsize=128):
    kw0 = kwidth/over_samp

    kr = kwidth/2.0

    beta = np.pi*np.sqrt((kw0*(over_samp-0.5))**2-0.8)

    # print beta

    kosf = np.floor(0.91/(over_samp*1e-3));
    # print kr
    # print kosf
    # print kosf*kr
    om = np.arange(kosf*kr)/(kosf*kr)
    # print om

    x = np.linspace(0, kr, grid_mod)

    x_bess = np.sqrt(1 - (x / kr) ** 2)
    # print x_bess
    y = np.i0(beta * x_bess)
    y = y / y[0]

    # pl.plot(y)
    # pl.show()

    fx = np.arange(imsize*over_samp)
    fx = fx - imsize*over_samp/2
    fx = fx / imsize * 1.5
    # print fx

    sqa = np.sqrt((np.pi*np.pi*kw0*kw0*fx*fx-beta*beta).astype('c16'))
    # print sqa
    fy = np.sin(sqa)/sqa
    fy = abs(fy)
    fy = fy / fy.max()
    # print fy

    # pl.plot(fy)
    # pl.show()

    return (x,y,fx,fy)

    #
    # x = [-o*n/2:osf*n/2-1]/(n);
    # sqa = sqrt(pi*pi*kw*kw*x.*x-beta*beta);
    # dax = sin(sqa)./(sqa);
    # % normalize by DC value
    # dax = dax/dax(osf*n/2);

    # % compute kernel, assume e1 is 0.001, assuming nearest neighbor
    # kosf = floor(0.91/(osf*1e-3));
    #
    # % compute kernel
    # om = [0:kosf*kwidth]/(kosf*kwidth); % 2GKx/W in Beatty paper
    # p = besseli(0,beta*sqrt(1-om.*om));
    # p = p./p(1);
    #
    #
    #
    # # kernel = np.zeros(krad * grid_mod)
    # beta = np.pi * np.sqrt((kwidth * kwidth) / (over_samp * over_samp)
    #                        * (over_samp - 0.5) * (over_samp - 0.5) - 0.8)
    #
    # print beta
    #
    # x = np.linspace(0, kwidth, grid_mod)
    #
    # x1 = np.sqrt(1 - (x / kwidth) ** 2)
    #
    # y = np.i0(beta * x1)
    # y = y / y[0]
    #
    # # pl.plot(x, y)
    # # pl.show()
    #
    # fx = np.arange(192)
    # fx = fx - 192/2.0
    # fy = np.zeros(192, 'c16')
    #
    # deapp = np.sqrt((np.pi * fx / kwidth)**2 - beta**2)
    #
    # fy = np.sin(deapp)/deapp
    #
    # return(x, y, fy)
