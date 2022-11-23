import numpy as np
import scipy.sparse.linalg as spspla

from jive.solver.solver import Solver
from jive.solver.constrainer import Constrainer

MAXITER = 'maxIter'

class StdIterativeSolver(Solver):

    def __init__(self):
        super().__init__()

        self._matrix = None
        self._cons = None
        self._conman = None
        self._precon = None

        self._maxiter = 1000

    def configure(self, props):
        self._maxiter = props.get(MAXITER, self._maxiter)

    def update(self, matrix, constraints, preconditioner=None):
        self._cons = constraints
        self._conman = Constrainer(self._cons, matrix)
        self._matrix = self._conman.get_output_matrix()
        self._precon = preconditioner

    def start(self):
        if self._precon is not None:
            self._precon.start()

    def finish(self):
        if self._precon is not None:
            self._precon.finish()

    def improve(self, lhs, rhs):
        f = self._conman.get_rhs(rhs)
        u = self._conman.get_lhs(lhs)

        for iter in range(self._maxiter):
            res = self.get_residual(u, f)
            error = np.linalg.norm(res)

            if error < self._precision:
                break

            du = self._solve(res)
            u += du
        else:
            raise RuntimeError('maximum number of iterations {} exceeded'.format(self._maxiter))

        return u

    def _solve(self, res):
        return spspla.spsolve(self._matrix, -res)

    def get_residual(self, lhs, rhs):
        return self._matrix @ lhs - rhs

    def get_matrix(self):
        return self._matrix

    def get_constraints(self):
        return self._cons


def declare(factory):
    factory.declare_solver('StdIterativeSolver', StdIterativeSolver)
