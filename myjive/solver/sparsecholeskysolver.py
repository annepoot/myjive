from .directsolver import DirectSolver
from .constrainer import Constrainer
from .util import reorder as reord
from .jit.cholesky import sparse_cholesky
from .jit.spsolve import solve_triangular

__all__ = ["SparseCholeskySolver"]


class SparseCholeskySolver(DirectSolver):
    def __init__(self):
        super().__init__()

        self._P = None
        self._L = None
        self._LT = None

    def update(self, matrix, constraints, preconditioner=None):
        self._cons = constraints
        self._conman = Constrainer(self._cons, matrix)
        self._matrix = self._conman.get_output_matrix()

        self._P = reord.get_reorder(self._matrix)
        self._matrix_r = reord.reorder_matrix(self._matrix, self._P)
        self._L = sparse_cholesky(self._matrix_r)
        self._LT = self._L.T.tocsr()

    def improve(self, lhs, rhs):
        if self.precon_mode:
            f = rhs
        else:
            f = self._conman.get_rhs(rhs)

        f_r = reord.reorder_vector(f, self._P)

        tmp = solve_triangular(self._L, f_r, lower=True)
        lhs_r = solve_triangular(self._LT, tmp, lower=False)

        lhs = reord.rev_reorder_vector(lhs_r, self._P)

        return lhs

    def get_matrix(self):
        return self._matrix

    def get_constraints(self):
        return self._cons
