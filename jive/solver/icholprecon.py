import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spspla

class ICholPrecon:

    def __init__(self):
        super().__init__()

        self._sourcematrix = None
        self._matrix = None

    def update(self, sourcematrix):
        self._sourcematrix = sourcematrix
        self._L = self.incomplete_cholesky(self._sourcematrix)

    def dot(self, lhs):
        return self._L @ (self._L.T @ lhs)

    def solve(self, rhs):
        tmp = spspla.spsolve_triangular(self._L, rhs, lower=True)
        return spspla.spsolve_triangular(self._L.T, tmp, lower=False)

    def get_matrix(self):
        return self._L @ self._L.T

    def incomplete_cholesky(self, K):

        # Get the lower triangle of A (for the sparsity structure)
        L = spsp.tril(K, format='csr')

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
                L[row,col] = np.sqrt(K[row,col] - rowsum)
            else:
                L[row,col] = (K[row,col] - rowsum) / L[col,col]

        return L

def declare(factory):
    factory.declare_precon('ichol', ICholPrecon)
