""" This is a pure python/numpy implementation using vectors.

This is probably the fastest non-C, non-CUDA, non-sparse method of gridding.
Make sure you are using a proper math library linked with numpy though so that
the vector math is efficient.
"""

import numpy as np
import pylab as pl

class GMPythonVectors():

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
                    dy = abs(y - iy)
                    kvy = kernel.get_kval(dy)
                    for ix in range(xmin, xmax + 1):
                        if ix >= 0 and ix < kspace.shape[0]:
                            dx = abs(x - ix)
                            kvx = kernel.get_kval(dx)
                            # dr = np.sqrt(dx * dx + dy * dy)
                            # kval = kernel.get_kval(dr)
                            data.flat[i] += kvx * kvy * kspace[ix, iy]
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
                    dy = abs(y - iy)
                    kvy = kernel.get_kval(dy)
                    for ix in range(xmin, xmax + 1):
                        if ix >= 0 and ix < kspace.shape[0]:
                            dx = abs(x - ix)
                            kvx = kernel.get_kval(dx)
                            # dr = np.sqrt(dx * dx + dy * dy)
                            # kval = kernel.get_kval(dr)
                            kspace[ix, iy] += data.flat[i] * density.flat[i] * kvx * kvy
        return kspace
