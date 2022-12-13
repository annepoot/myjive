from jive.solver.preconditioner import Preconditioner

from jive.solver.numba.cholesky import incomplete_cholesky
from jive.solver.numba.spsolve import solve_triangular

class ICholPrecon(Preconditioner):

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
