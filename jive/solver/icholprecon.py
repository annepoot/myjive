import numpy as np
import scipy.sparse as spsp
from numba import njit

class ICholPrecon:

    def __init__(self):
        super().__init__()

        self._sourcematrix = None
        self._L = None
        self._LT = None

    def update(self, sourcematrix):
        self._sourcematrix = sourcematrix
        self._L = incomplete_cholesky(self._sourcematrix)
        self._LT = self._L.T.tocsr()

    def dot(self, lhs):
        return self._L @ (self._L.T @ lhs)

    def solve(self, rhs):
        tmp = solve_triangular(self._L, rhs, lower=True)
        return solve_triangular(self._LT, tmp, lower=False)

    def get_matrix(self):
        return self._L @ self._L.T


def declare(factory):
    factory.declare_precon('ichol', ICholPrecon)


#####################
# wrapper functions #
#####################

def incomplete_cholesky(A):
    if not spsp.isspmatrix_csr(A):
        raise ValueError('A has to be a sparse matrix in csr format')

    L = spsp.tril(A, format='csr')
    L.data = incomplete_cholesky_jit(L.data, L.indices, L.indptr)

    return L


def solve_triangular(A, b, lower=True):
    if not spsp.isspmatrix_csr(A):
        raise ValueError('A has to be a sparse matrix in csr format')

    return solve_triangular_jit(A.data, A.indices, A.indptr, b, lower)


##########################
# numba helper functions #
##########################

@njit
def idx2rowcol(indices, indptr, idx):
    for r, ipt in enumerate(indptr):
        if ipt > idx:
            row = r - 1
            break
    col = indices[idx]
    return row, col


@njit
def rowcol2idx(indices, indptr, row, col):
    cols = indices[indptr[row]:indptr[row+1]]
    for i, c in enumerate(cols):
        if c == col:
            idx = i
            break
    idx = indptr[row] + idx
    return idx


@njit
def idxs2rowscols(indices, indptr):
    rows = np.zeros_like(indices)
    cols = np.zeros_like(indices)
    idx = 0
    for row in range(len(indptr)-1):
        for col in indices[indptr[row]:indptr[row+1]]:
            rows[idx] = row
            cols[idx] = col
            idx += 1
    assert idx == len(indices)
    return rows, cols


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
