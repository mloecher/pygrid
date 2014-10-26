import numpy as np

def zeropad_ratio(im, ratio):
    shape0 = im.shape
    shape1 = tuple([x * ratio for x in shape0])

    print shape0
    print shape1

    out = np.zeros(shape1)

    d0_start = shape1[0]/2 - shape0[0]/2
    d0_stop = shape1[0]/2 + shape0[0]/2
    d1_start = shape1[1]/2 - shape0[1]/2
    d1_stop = shape1[1]/2 + shape0[1]/2


    out[d0_start:d0_stop, d1_start:d1_stop] = im
    return out
