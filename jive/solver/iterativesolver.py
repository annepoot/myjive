import numpy as np
import warnings

from jive.solver.solver import Solver
from jive.solver.constrainer import Constrainer

MAXITER = 'maxIter'
ALLOWMAXITER = 'allowMaxIter'
NOTIMPLEMENTEDMSG = 'this function needs to be implemented in an derived class'

class IterativeSolver(Solver):

    def __init__(self):
        super().__init__()

        self._matrix = None
        self._cons = None
        self._conman = None
        self._precon = None

        self._maxiter = 2000
        self._allowmaxiter = False

    def configure(self, props, globdat):
        super().configure(props, globdat)
        self._maxiter = int(props.get(MAXITER, self._maxiter))
        if ALLOWMAXITER in props:
            self._allowmaxiter = bool(eval(props[ALLOWMAXITER]))

    def update(self, matrix, constraints, preconditioner=None):
        self._cons = constraints
        self._conman = Constrainer(self._cons, matrix)
        self._matrix = self._conman.get_output_matrix()
        self._precon = preconditioner
        if self._precon is not None:
            self._precon.update(self._matrix)

    def improve(self, lhs, rhs):
        if self.precon_mode:
            f = rhs
            u = lhs
        else:
            f = self._conman.get_rhs(rhs)
            u = self._conman.get_lhs(lhs)

        self.start()

        for _ in range(self._maxiter):
            res = self.get_residual(u, f)
            error = np.linalg.norm(res)

            if error < self._precision:
                break

            du = self.iterate(res)
            u += du
        else:
            if self._allowmaxiter:
                warnings.warn('maximum number of iterations {} exceeded'.format(self._maxiter), RuntimeWarning)
            else:
                raise RuntimeError('maximum number of iterations {} exceeded'.format(self._maxiter))

        self.finish()

        return u

    def iterate(self, res):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def start(self):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def finish(self):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def get_residual(self, lhs, rhs):
        return self._matrix @ lhs - rhs

    def get_matrix(self):
        return self._matrix

    def get_constraints(self):
        return self._cons


def declare(factory):
    factory.declare_solver('iterative', IterativeSolver)
