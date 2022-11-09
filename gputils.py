import numpy as np
import scipy.sparse as spsp

def cholesky(A, get_nops=False):
    L = np.zeros_like(A)
    n_ops = 0

    for i in range(A.shape[0]):
        for j in range(i+1):
            rowsum = 0

            for k in range(j):
                n_ops += 1
                rowsum += L[i,k] * L[j,k]

            if i == j:
                L[i,j] = np.sqrt(A[i,i] - rowsum)
            else:
                L[i,j] = 1 / L[j,j] * (A[i,j] - rowsum)

    if get_nops:
        return L, n_ops
    else:
        return L

def incomplete_cholesky(A, get_nops=False):
    if not spsp.isspmatrix_csr(A):
        raise TypeError('A must be a scipy csr_array')

    # Get the lower triangle of A (for the sparsity structure)
    L = spsp.tril(A, format='csr')
    n_ops = 0

    # Get the data of L
    indptr = L.indptr
    indices = L.indices
    data = L.data

    # Get all row and column indices in pairs
    for row, col in zip(*L.nonzero()):

        # Get all entries belonging to row i and row j
        irowindices = indices[indptr[row]:indptr[row+1]]
        irowvalues = data[indptr[row]:indptr[row+1]]
        jrowindices = indices[indptr[col]:indptr[col+1]]
        jrowvalues = data[indptr[col]:indptr[col+1]]

        # Initialize rowsum computation
        rowsum = 0
        iidx = 0
        jidx = 0

        # Compute sum(L_ik * Ljk) for 0 <= k < j
        while iidx < len(irowindices) and jidx < len(jrowindices):
            n_ops += 1

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
        n_ops += 1
        if row == col:
            L[row,col] = np.sqrt(A[row,col] - rowsum)
        else:
            L[row,col] = (A[row,col] - rowsum) / L[col,col]

    if get_nops:
        return L, n_ops
    else:
        return L

def sparse_cholesky(A, get_nops=False):
    if not spsp.isspmatrix_csr(A):
        raise TypeError('A must be a scipy csr_array')

    # Get the lower triangle of A (for the sparsity structure)
    Ldata = []
    Lindices = []
    Lindptr = [0]
    n_ops = 0

    # Go over all rows, and all relevant columns
    # (starting from the first non-zero column in that row)
    for row in range(A.shape[0]):
        for col in range(A.indices[A.indptr[row]], row+1):

            # Get all entries belonging to row i and row j
            irowindices = Lindices[Lindptr[row]:]
            irowvalues = Ldata[Lindptr[row]:]

            if row == col:
                jrowindices = irowindices
                jrowvalues = irowvalues
            else:
                jrowindices = Lindices[Lindptr[col]:Lindptr[col+1]]
                jrowvalues = Ldata[Lindptr[col]:Lindptr[col+1]]

            # Initialize rowsum computation
            rowsum = 0
            iidx = 0
            jidx = 0

            # Compute sum(L_ik * Ljk) for 0 <= k < j
            while iidx < len(irowindices) and jidx < len(jrowindices):
                n_ops += 1

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
            n_ops += 1
            if row == col:
                Lij = np.sqrt(A[row,col] - rowsum)
            else:
                Ljj = Ldata[Lindptr[col+1]-1]
                Lij = (A[row,col] - rowsum) / Ljj

            if not np.isclose(Lij, 0):
                Ldata.append(Lij)
                Lindices.append(col)

        Lindptr.append(len(Lindices))

    L = spsp.csr_array((Ldata, Lindices, Lindptr), dtype=float)

    if get_nops:
        return L, n_ops
    else:
        return L
