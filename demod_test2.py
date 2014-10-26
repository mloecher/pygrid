from kernels import kaiser_bessel, kaiser_bessel2
from trajectories import radial_2d
from grid_methods import simple_igrid2d
from utils import *
from pygrid import Gridder
import numpy as np
import pylab as pl
import sys

def full_demod(demod, os):
    x = np.arange(128*os) - 128*os/2
    y = np.arange(128*os) - 128*os/2
    xx, yy = np.meshgrid(x, y)
    rr = np.sqrt(xx * xx + yy * yy)
    rr = np.round(rr)
    rr[rr>128/2] = 128/2

    out = demod[rr.flatten().astype('int')]

    out = out.reshape((128*os, 128*os))
    #
    # pl.figure()
    # pl.imshow(out)
    # pl.colorbar()

    rr = np.sqrt(xx * xx + yy * yy)
    rr[rr>128*os/2-10] = 128*os/2-10

    rr0 = np.floor(rr)
    rr1 = np.ceil(rr)

    y0 = demod[rr0.flatten().astype('int')]
    y1 = demod[rr1.flatten().astype('int')]

    drr = (rr - rr0).flatten()

    out = (1-drr)*y0 + drr*y1

    out = out.reshape((128*os, 128*os))

    print out.min()

    # out[out < .1] = .1

    # pl.figure()
    # pl.imshow(out)
    # pl.colorbar()
    # pl.show()

    return out


os = 1.5
x, y, fx, fy = kaiser_bessel2(3.0, 100, os)


bigsize = 128*os

Fx = np.arange(bigsize)
Fx = Fx - bigsize/2.0
Fy = np.zeros(bigsize, 'c16')

print y

for i in range(1, x.size):
    temp = y[i] * 2*np.exp(2*1j*np.pi*(Fx/128)*(x[i]))
    Fy += temp
Fy = Fy.real
print (Fy.min(), Fy.max())
Fy = Fy - Fy.min()
print (Fy.min(), Fy.max())
Fy = Fy/Fy.max()
print (Fy.min(), Fy.max())

Fx2 = np.arange(bigsize*2)
Fx2 = Fx2 - bigsize*2/2.0
Fy2 = np.zeros(bigsize*2, 'c16')

for i in range(1, x.size):
    temp = y[i] * 2*np.exp(2*1j*np.pi*(Fx2/128)*(x[i]))
    Fy2 += temp

Fy2 = Fy2.real / Fy2.real.max()

# pl.plot(Fx2,Fy2)
# pl.plot(Fx,Fy)
# pl.show()
#
# sys.exit()

demod = full_demod(Fy[Fy.size/2:], os)

# im = np.ones((128, 128))

x = np.linspace(-.5, .5, 128)
y = np.linspace(-.5, .5, 128)
xx, yy = np.meshgrid(x, y)
rr = np.sqrt(xx * xx + yy * yy)
im = np.zeros((128, 128))
im[rr < .4] = 2.0
im[rr < .3] = 1.0


im = zeropad_ratio(im, os)

kspace = np.fft.fftshift(im)
kspace = np.fft.fft2(kspace)
kspace = np.fft.ifftshift(kspace)

(traj, dens) = radial_2d(200, 161)

gridder = Gridder()
data = gridder.igrid_2d(kspace, traj)
k1 = gridder.grid_2d(data, traj, dens)

im1 = np.fft.fftshift(k1)
im1 = np.fft.ifft2(im1)
im1 = np.fft.ifftshift(im1)

im1 = im1/demod


im1 = crop_ratio(im1, os)

pl.imshow(abs(im1))
pl.show()

sys.exit()

demod = abs(im1)/abs(im)
center = demod.shape[0]/2
demod = demod/demod[center, center]

diff = demod[:,center]/fy

pl.plot(demod[:,center])
pl.plot(fy)
# pl.plot(diff)
pl.show()


sys.exit()

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
