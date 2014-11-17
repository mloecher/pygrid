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

        x = (trajectory[0].flatten() + .5) * kspace.shape[0]
        y = (trajectory[1].flatten() + .5) * kspace.shape[1]

        iwin = int(np.floor(kernel.krad))

        for ymod in range(-iwin, iwin+1):
            iy = np.round(y+ymod).astype('int')
            indy = np.where((iy >= 0) & (iy < kspace.shape[1]))[0]
            dy = np.abs(y-iy)
            kvy = kernel.get_kval_vec(dy)

            for xmod in range(-iwin, iwin+1):
                ix = np.round(x+xmod).astype('int')
                indx = np.where((ix >= 0) & (ix < kspace.shape[0]))[0]
                dx = np.abs(x-ix)
                kvx = kernel.get_kval_vec(dx)

                ind = np.intersect1d(indx, indy)

                data.flat[ind] += kvy[ind] * kvx[ind] * kspace[ix[ind], iy[ind]]

        return data

    def grid_2d(self, data, grid_params, kernel, trajectory, density=1.0):

        ksize = tuple([x * grid_params.over_samp for x in grid_params.imsize])
        kspace = np.zeros(ksize, 'c16')

        x = (trajectory[0].flatten() + .5) * kspace.shape[0]
        y = (trajectory[1].flatten() + .5) * kspace.shape[1]

        iwin = int(np.floor(kernel.krad))

        for ymod in range(-iwin, iwin+1):
            iy = np.round(y+ymod).astype('int')
            indy = np.where((iy >= 0) & (iy < kspace.shape[1]))[0]
            dy = np.abs(y-iy)
            kvy = kernel.get_kval_vec(dy)

            for xmod in range(-iwin, iwin+1):
                ix = np.round(x+xmod).astype('int')
                indx = np.where((ix >= 0) & (ix < kspace.shape[0]))[0]
                dx = np.abs(x-ix)
                kvx = kernel.get_kval_vec(dx)

                ind = np.intersect1d(indx, indy)

                lin_ind = np.ravel_multi_index((ix[ind], iy[ind]), kspace.shape)
                w = kvy[ind] * kvx[ind] * density.flat[ind]
                accum_r = np.bincount(lin_ind, w * data.flat[ind].real, kspace.size)
                accum_i = np.bincount(lin_ind, w * data.flat[ind].imag, kspace.size)

                kspace.flat += (accum_r + 1j*accum_i)
        return kspace
