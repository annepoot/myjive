import numpy as np
from numba import njit


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
    cols = indices[indptr[row] : indptr[row + 1]]
    for i, c in enumerate(cols):
        if c > col:
            return -1
        elif c == col:
            idx = i
            break
    idx = indptr[row] + idx
    return idx


@njit
def idxs2rowscols(indices, indptr):
    rows = np.zeros_like(indices)
    cols = np.zeros_like(indices)
    idx = 0
    for row in range(len(indptr) - 1):
        for col in indices[indptr[row] : indptr[row + 1]]:
            rows[idx] = row
            cols[idx] = col
            idx += 1
    assert idx == len(indices)
    return rows, cols
