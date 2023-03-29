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
        J = glob_coords @ dN[:, :, ip]
        wts[ip] *= np.linalg.det(J)
        dN[:, :, ip] = dN[:, :, ip] @ np.linalg.inv(J)

    return dN, wts
