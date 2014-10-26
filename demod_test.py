from kernels import kaiser_bessel
from trajectories import radial_2d
from grid_methods import simple_igrid2d
from utils import *
from pygrid import Gridder
import numpy as np
import pylab as pl

x, y, fy22 = kaiser_bessel(2.5, 100, 1.5)


pl.plot(fy22)
pl.show()

fx = np.arange(192)
fx = fx - 192/2.0
fy = np.zeros(192, 'c16')

print fx

for i in range(1, x.size):
    temp = y[i] * 2*np.cos(2*np.pi*(fx/np.sqrt(fx.size))*(x[i]/np.sqrt(x.size)))
    fy += temp
    # pl.plot(temp)
    # pl.show()


im = np.ones((128, 128))

im = zeropad_ratio(im, 1.5)

kspace = np.fft.fftshift(im)
kspace = np.fft.fft2(kspace)
kspace = np.fft.ifftshift(kspace)

(traj, dens) = radial_2d(30, 161)

gridder = Gridder()
data = gridder.igrid_2d(kspace, traj)
k1 = gridder.grid_2d(data, traj, dens)

im1 = np.fft.fftshift(k1)
im1 = np.fft.ifft2(im1)
im1 = np.fft.ifftshift(im1)

demod = abs(im1)/abs(im)
demod = demod/demod[96, 96]
fy = fy / fy.max()
diff = demod[:,96]/fy
# pl.plot(demod[:,96])
pl.plot(fy)
# pl.plot(diff)
# pl.plot(abs(fy2)/abs(fy2).max())
pl.show()

# pl.imshow(abs(demod))
# pl.show()

# pl.plot(abs(fy))
# pl.plot(abs(fy2))
# pl.show()
