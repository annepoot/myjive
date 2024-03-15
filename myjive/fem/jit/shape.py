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
        dNip = np.copy(dN[:, :, ip])
        J = glob_coords @ dNip
        wts[ip] *= np.linalg.det(J)
        dNip = dNip @ np.linalg.inv(J)
        dN[:, :, ip] = dNip

    return dN, wts
