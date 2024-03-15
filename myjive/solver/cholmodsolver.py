from sksparse import cholmod as cm

from .sparsecholeskysolver import SparseCholeskySolver
from .constrainer import Constrainer

__all__ = ["CholmodSolver"]


class CholmodSolver(SparseCholeskySolver):
    def update(self, matrix, constraints, preconditioner=None):
        self._cons = constraints
        self._conman = Constrainer(self._cons, matrix)
        self._matrix = self._conman.get_output_matrix()

        # Transpose to convert the matrix from csr to csc
        self._chol = cm.cholesky(
            self._matrix.T, mode="simplicial", ordering_method="amd"
        )

    def improve(self, lhs, rhs):
        if self.precon_mode:
            f = rhs
        else:
            f = self._conman.get_rhs(rhs)

        lhs = self._chol.solve_A(f)

        return lhs

    def get_matrix(self):
        return self._matrix

    def get_constraints(self):
        return self._cons
