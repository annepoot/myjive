import numpy as np
import scipy.sparse as spsp
from numba import njit

#####################
# wrapper functions #
#####################

def solve_triangular(A, b, lower=True):
    if not spsp.isspmatrix_csr(A):
        raise ValueError('A has to be a sparse matrix in csr format')

    return solve_triangular_jit(A.data, A.indices, A.indptr, b, lower)


##########################
# numba helper functions #
##########################

@njit
def solve_triangular_jit(data, indices, indptr, b, lower=True):
    x = np.zeros_like(b)

    if lower:
        rng = range(len(b))
    else:
        rng = range(len(b)-1, -1, -1)

    for i in rng:
        if lower:
            start = indptr[i]
            end = indptr[i+1]-1
            L_ii = data[indptr[i+1]-1]
        else:
            start = indptr[i]+1
            end = indptr[i+1]
            L_ii = data[indptr[i]]

        cols = indices[start:end]
        vals = data[start:end]
        numerator = b[i]

        for j, L_ij in zip(cols, vals):
            numerator -= L_ij * x[j]
        x[i] = numerator / L_ii
    return x
