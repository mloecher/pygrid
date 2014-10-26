""" This is a simle gridding method with a bunch of for loops in pure python.

It really shouldn't be used for anything becayse it is so slow, just testing and
such.
"""

import numpy as np
import pylab as pl

class GMPythonForLoops():

    def setup(self):
        pass

    def igrid_2d(self, kspace, grid_params, kernel, trajectory):

        data = np.zeros(trajectory[0].shape, 'c16')

        for i in range(data.size):
            x = (trajectory[0].flat[i] + .5) * kspace.shape[0]
            y = (trajectory[1].flat[i] + .5) * kspace.shape[1]
            xmin = int(np.floor(x - kernel.krad))
            ymin = int(np.floor(y - kernel.krad))
            xmax = int(np.ceil(x + kernel.krad))
            ymax = int(np.ceil(y + kernel.krad))
            for iy in range(ymin, ymax + 1):
                if iy >= 0 and iy < kspace.shape[1]:
                    dy = y - iy
                    for ix in range(xmin, xmax + 1):
                        if ix >= 0 and ix < kspace.shape[0]:
                            dx = x - ix
                            dr = np.sqrt(dx * dx + dy * dy)
                            kval = kernel.get_kval(dr)
                            data.flat[i] += kval * kspace[ix, iy]
        return data

    def grid_2d(self, data, grid_params, kernel, trajectory, density=1.0):

        ksize = tuple([x * grid_params.over_samp for x in grid_params.imsize])
        kspace = np.zeros(ksize, 'c16')

        for i in range(data.size):
            x = (trajectory[0].flat[i] + .5) * kspace.shape[0]
            y = (trajectory[1].flat[i] + .5) * kspace.shape[1]
            xmin = int(np.floor(x - kernel.krad))
            ymin = int(np.floor(y - kernel.krad))
            xmax = int(np.ceil(x + kernel.krad))
            ymax = int(np.ceil(y + kernel.krad))
            for iy in range(ymin, ymax + 1):
                if iy >= 0 and iy < kspace.shape[1]:
                    dy = y - iy
                    for ix in range(xmin, xmax + 1):
                        if ix >= 0 and ix < kspace.shape[0]:
                            dx = x - ix
                            dr = np.sqrt(dx * dx + dy * dy)
                            kval = kernel.get_kval(dr)
                            kspace[ix, iy] += data.flat[i] * density.flat[i] * kval
        return kspace
