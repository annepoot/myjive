import scipy.sparse as spsp

from .preconditioner import Preconditioner

__all__ = ["IdPrecon"]


class IdPrecon(Preconditioner):
    def __init__(self):
        super().__init__()

        self._sourcematrix = None
        self._M = None

    def update(self, sourcematrix):
        self._sourcematrix = sourcematrix
        self._M = spsp.identity(self._sourcematrix.shape[0], format="csr")

    def dot(self, lhs):
        return lhs

    def solve(self, rhs):
        return rhs

    def get_matrix(self):
        return self._M

    def get_inv_matrix(self):
        return self._M
