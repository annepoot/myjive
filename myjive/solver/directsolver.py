import scipy.sparse.linalg as spspla

from .solver import Solver
from .constrainer import Constrainer

__all__ = ["DirectSolver"]


class DirectSolver(Solver):
    def __init__(self):
        super().__init__()

        self._matrix = None
        self._cons = None
        self._conman = None

    def update(self, matrix, constraints, preconditioner=None):
        self._cons = constraints
        self._conman = Constrainer(self._cons, matrix)
        self._matrix = self._conman.get_output_matrix()

    def improve(self, lhs, rhs):
        if self.precon_mode:
            f = rhs
        else:
            f = self._conman.get_rhs(rhs)

        lhs = spspla.spsolve(self._matrix, f)
        return lhs

    def get_matrix(self):
        return self._matrix

    def get_constraints(self):
        return self._cons
