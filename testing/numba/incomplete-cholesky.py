import sys
sys.path.append('../../')

import numpy as np
import scipy.sparse as spsp
from numba import njit
from time import perf_counter

import jive.util.proputils as pu
from jive.app import main
import gputils as gpu

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
def incomplete_cholesky(data, indices, indptr):

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

props = pu.parse_file('beam.pro')
props['init']['mesh']['file'] = 'meshes/beam_fine4.msh'
props['solver']['solver'] = 'direct'
globdat = main.jive(props)

K = globdat['matrix0']

time0 = perf_counter()
L = spsp.tril(K, format='csr')
time0 = perf_counter() - time0

time1 = perf_counter()
data = incomplete_cholesky(L.data, L.indices, L.indptr)
time1 = perf_counter() - time1

time2 = perf_counter()
data = incomplete_cholesky(L.data, L.indices, L.indptr)
time2 = perf_counter() - time2

time0 = time0 - perf_counter()
L.data = data
time0 = time0 + perf_counter()

time3 = perf_counter()
Lgpu = gpu.incomplete_cholesky(K)
time3 = perf_counter() - time3

print('Reference time (gpu.incomplete_cholesky):', time3)
print('Time including compilation:', time1)
print('Time excluding compilation:', time2)
print('Time for remaining tasks:', time0)
print('Speedup factor:', time3 / time2)

assert np.isclose(L.data, Lgpu.data).all()
