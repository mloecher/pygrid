import numpy as np

def zeropad_ratio(im, ratio):
    shape0 = im.shape
    shape1 = tuple([x * ratio for x in shape0])

    out = np.zeros(shape1)

    d0_start = shape1[0]/2 - shape0[0]/2
    d0_stop = shape1[0]/2 + shape0[0]/2
    d1_start = shape1[1]/2 - shape0[1]/2
    d1_stop = shape1[1]/2 + shape0[1]/2


    out[d0_start:d0_stop, d1_start:d1_stop] = im
    return out

def crop_ratio(im, ratio):
    shape_big = im.shape
    shape_small = tuple([x / ratio for x in shape_big])

    out = np.zeros(shape_small)

    d0_start = shape_big[0]/2 - shape_small[0]/2
    d0_stop = shape_big[0]/2 + shape_small[0]/2
    d1_start = shape_big[1]/2 - shape_small[1]/2
    d1_stop = shape_big[1]/2 + shape_small[1]/2


    out = im[d0_start:d0_stop, d1_start:d1_stop]
    return out
