import numpy as np
import scipy.sparse as spsp
from numba import njit

from .sputil import rowcol2idx, idxs2rowscols

#####################
# wrapper functions #
#####################


def incomplete_cholesky(A):
    if not spsp.isspmatrix_csr(A):
        raise ValueError("A has to be a sparse matrix in csr format")

    L = spsp.tril(A, format="csr")
    L.data = incomplete_cholesky_jit(L.data, L.indices, L.indptr)

    return L


def sparse_cholesky(A, mode="LDL"):
    if not spsp.isspmatrix_csr(A):
        raise ValueError("A has to be a sparse matrix in csr format")

    A.sort_indices()

    if mode == "LL":
        L = sparse_cholesky_LL(A)
    elif mode == "LDL":
        L = sparse_cholesky_LDL(A)

    L.sort_indices()

    return L


def sparse_cholesky_LL(A):
    Ldata, Lindices, Lindptr = sparse_cholesky_jit(A.data, A.indices, A.indptr)
    L = spsp.csr_array((Ldata, Lindices, Lindptr), dtype=A.dtype)
    return L


def sparse_cholesky_LDL(A):
    Ldata, Lindices, Lindptr, Ddata = sparse_LDL_jit(A.data, A.indices, A.indptr)
    L = spsp.csr_array((Ldata, Lindices, Lindptr), dtype=A.dtype)
    sqrtD = spsp.diags(np.sqrt(Ddata), shape=L.shape)
    return L @ sqrtD


##########################
# numba helper functions #
##########################


@njit
def incomplete_cholesky_jit(data, indices, indptr):
    L_data = data.copy()
    rows, cols = idxs2rowscols(indices, indptr)

    # Get all row and column indices in pairs
    for idx, (row, col) in enumerate(zip(rows, cols)):
        if data[idx] == 0.0:
            continue

        A_ij = data[idx]

        # Get all entries belonging to row i and row j
        irowindices = indices[indptr[row] : indptr[row + 1]]
        irowvalues = L_data[indptr[row] : indptr[row + 1]]
        jrowindices = indices[indptr[col] : indptr[col + 1]]
        jrowvalues = L_data[indptr[col] : indptr[col + 1]]

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
            if A_ij - rowsum <= 0:
                raise ValueError("Matrix is not positive definite")

            L_ij = np.sqrt(A_ij - rowsum)
        else:
            idx_jj = indptr[col + 1] - 1
            L_jj = L_data[idx_jj]
            L_ij = (A_ij - rowsum) / L_jj

        L_data[idx] = L_ij

    return L_data


@njit
def sparse_cholesky_jit(data, indices, indptr):
    Ldata = []
    Lindices = []
    Lindptr = [0]

    # Go over all rows, and all relevant columns
    # (starting from the first non-zero column in that row)
    for row in range(len(indptr) - 1):
        for col in range(indices[indptr[row]], row + 1):
            idx = rowcol2idx(indices, indptr, row, col)
            if idx >= 0:
                A_ij = data[idx]
            else:
                A_ij = 0

            # Get all entries belonging to row i and row j
            irowindices = Lindices[Lindptr[row] :]
            irowvalues = Ldata[Lindptr[row] :]

            if row == col:
                jrowindices = irowindices
                jrowvalues = irowvalues
            else:
                jrowindices = Lindices[Lindptr[col] : Lindptr[col + 1]]
                jrowvalues = Ldata[Lindptr[col] : Lindptr[col + 1]]

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
                if A_ij - rowsum <= 0:
                    raise ValueError("Matrix is not positive definite")

                Lij = np.sqrt(A_ij - rowsum)

            else:
                idx_jj = Lindptr[col + 1] - 1
                Ljj = Ldata[idx_jj]
                Lij = (A_ij - rowsum) / Ljj

            if Lij != 0:
                Ldata.append(Lij)
                Lindices.append(col)

        Lindptr.append(len(Lindices))

    return Ldata, Lindices, Lindptr


@njit
def sparse_LDL_jit(data, indices, indptr):
    Ldata = []
    Lindices = []
    Lindptr = [0]
    Ddata = np.zeros(len(indptr) - 1)

    # Go over all rows, and all relevant columns
    # (starting from the first non-zero column in that row)
    for row in range(len(indptr) - 1):
        for col in range(indices[indptr[row]], row + 1):
            idx = rowcol2idx(indices, indptr, row, col)
            if idx >= 0:
                A_ij = data[idx]
            else:
                A_ij = 0

            # Get all entries belonging to row i and row j
            irowindices = Lindices[Lindptr[row] :]
            irowvalues = Ldata[Lindptr[row] :]

            if row == col:
                jrowindices = irowindices
                jrowvalues = irowvalues
            else:
                jrowindices = Lindices[Lindptr[col] : Lindptr[col + 1]]
                jrowvalues = Ldata[Lindptr[col] : Lindptr[col + 1]]

            # Initialize rowsum computation
            rowsum = 0.0
            iidx = 0
            jidx = 0

            # Compute sum(L_ik * Ljk * Dk) for 0 <= k < j
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
                    rowsum += irowvalues[iidx] * jrowvalues[jidx] * Ddata[icol]
                    iidx += 1
                    jidx += 1

            # Compute the next entry in the lower triangular matrix
            if row == col:
                Lij = 1.0
                Ddata[row] = A_ij - rowsum

            else:
                Dj = Ddata[col]
                Lij = (A_ij - rowsum) / Dj

            if Lij != 0:
                Ldata.append(Lij)
                Lindices.append(col)

        if Ddata[len(Lindptr) - 1] <= 0:
            raise ValueError("Matrix is not positive definite")

        Lindptr.append(len(Lindices))

    return Ldata, Lindices, Lindptr, Ddata
