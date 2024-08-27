import numpy as np
from numba import njit

##########################
# numba helper functions #
##########################


@njit
def get_shape_gradients_jit(glob_coords, _dN, _wts, _ipcount):
    wts = np.copy(_wts)
    dN = np.copy(_dN)

    for ip in range(_ipcount):
        dNip = dN[ip]
        J = dNip @ glob_coords
        wts[ip] *= np.linalg.det(J)
        dNip = np.linalg.inv(J) @ dNip
        dN[ip] = dNip

    return dN, wts
