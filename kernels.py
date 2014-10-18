import pylab as pl
import numpy as np


def kaiser_bessel(kwidth, grid_mod, over_samp=1.375):
    # kernel = np.zeros(krad * grid_mod)
    beta = np.pi * np.sqrt((kwidth * kwidth) / (over_samp * over_samp)
                           * (over_samp - 0.5) * (over_samp - 0.5) - 0.8)

    x = np.linspace(0, kwidth, grid_mod)

    x1 = np.sqrt(1-(x/kwidth)**2);

    y = np.i0(beta*x1)
    y = y / y[0]

    pl.plot(x,y)
    pl.show()

    return(x,y)

kaiser_bessel(2.5, 24)
