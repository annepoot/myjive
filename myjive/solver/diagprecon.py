import scipy.sparse as spsp

from .preconditioner import Preconditioner

__all__ = ["DiagPrecon"]


class DiagPrecon(Preconditioner):
    def __init__(self):
        super().__init__()

        self._sourcematrix = None
        self._M = None
        self._M_inv = None

    def update(self, sourcematrix):
        self._sourcematrix = sourcematrix
        diag = self._sourcematrix.diagonal()
        self._M = spsp.diags(diag, format="csr")
        self._M_inv = spsp.diags(1 / diag, format="csr")

    def dot(self, lhs):
        return self._M @ lhs

    def solve(self, rhs):
        return self._M_inv @ rhs

    def get_matrix(self):
        return self._M

    def get_inv_matrix(self):
        return self._M_inv
