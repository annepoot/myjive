import numpy as np
import scipy.sparse as spsp
from numba import njit

from jive.solver.numba.sputil import rowcol2idx, idxs2rowscols

#####################
# wrapper functions #
#####################

def incomplete_cholesky(A):
    if not spsp.isspmatrix_csr(A):
        raise ValueError('A has to be a sparse matrix in csr format')

    L = spsp.tril(A, format='csr')
    L.data = incomplete_cholesky_jit(L.data, L.indices, L.indptr)

    return L


##########################
# numba helper functions #
##########################

@njit
def incomplete_cholesky_jit(data, indices, indptr):
    L_data = data.copy()
    rows, cols = idxs2rowscols(indices, indptr)

    # Get all row and column indices in pairs
    for idx, row, col in zip(range(len(indices)), rows, cols):
        if data[idx] == 0.0:
            continue

        A_ij = data[idx]

        # Get all entries belonging to row i and row j
        irowindices = indices[indptr[row]:indptr[row+1]]
        irowvalues = L_data[indptr[row]:indptr[row+1]]
        jrowindices = indices[indptr[col]:indptr[col+1]]
        jrowvalues = L_data[indptr[col]:indptr[col+1]]

        # Initialize rowsum computation
        rowsum = 0.0
        iidx = 0
        jidx = 0

        # Compute sum(L_ik * Ljk) for 0 <= k < j
        while iidx < len(irowindices) and jidx < len(jrowindices):
            icol = irowindices[iidx]
            jcol = jrowindices[jidx]

            if icol >= col or jcol >= col:
                break

            if icol < jcol:
                iidx += 1
            elif icol > jcol:
                jidx += 1
            else:
                rowsum += irowvalues[iidx] * jrowvalues[jidx]
                iidx += 1
                jidx += 1

        # Compute the next entry in the lower triangular matrix
        if row == col:
            L_ij = np.sqrt(A_ij - rowsum)
        else:
            idx_jj = rowcol2idx(indices, indptr, col, col)
            L_jj = L_data[idx_jj]
            L_ij = (A_ij - rowsum) / L_jj

        L_data[idx] = L_ij

    return L_data
