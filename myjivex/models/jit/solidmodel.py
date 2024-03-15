import numpy as np
from numba import njit

##########################
# numba helper functions #
##########################


@njit
def get_B_matrix_jit(grads, _strcount, _dofcount, _nodecount, _rank):
    B = np.zeros((_strcount, _dofcount))
    if _rank == 1:
        B = grads.transpose()
    elif _rank == 2:
        for inode in range(_nodecount):
            i = 2 * inode
            gi = grads[inode, :]
            B[0, i + 0] = gi[0]
            B[1, i + 1] = gi[1]
            B[2, i + 0] = gi[1]
            B[2, i + 1] = gi[0]
    elif _rank == 3:
        for inode in range(_nodecount):
            i = 3 * inode
            gi = grads[inode, :]
            B[0, i + 0] = gi[0]
            B[1, i + 1] = gi[1]
            B[2, i + 2] = gi[2]
            B[3, i + 0] = gi[1]
            B[3, i + 1] = gi[0]
            B[4, i + 1] = gi[2]
            B[4, i + 2] = gi[1]
            B[5, i + 0] = gi[2]
            B[5, i + 2] = gi[0]
    return B


@njit
def get_N_matrix_jit(sfuncs, _dofcount, _rank):
    N = np.zeros((_rank, _dofcount))
    for i in range(_rank):
        N[i, i::_rank] = sfuncs.transpose()
    return N
