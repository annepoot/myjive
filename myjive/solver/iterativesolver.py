import numpy as np
import warnings

from myjive.util.proputils import optional_argument
from .solver import Solver
from .constrainer import Constrainer

NOTIMPLEMENTEDMSG = "this function needs to be implemented in an derived class"

__all__ = ["IterativeSolver"]


class IterativeSolver(Solver):
    def __init__(self):
        super().__init__()

        self._matrix = None
        self._cons = None
        self._conman = None
        self._precon = None

        self._init_guess = None

    def configure(self, globdat, **props):
        super().configure(globdat, **props)

        self._maxiter = optional_argument(self, props, "maxIter", default=2000)
        self._allow_max_iter = optional_argument(
            self, props, "allowMaxIter", default=False
        )

    def update(self, matrix, constraints, preconditioner=None):
        self._cons = constraints
        self._conman = Constrainer(self._cons, matrix)
        self._matrix = self._conman.get_output_matrix()
        self._precon = preconditioner
        if self._precon is not None:
            self._precon.update(self._matrix)

    def solve(self, rhs):
        if self._init_guess is None:
            lhs = np.zeros_like(rhs)
        else:
            lhs = self._init_guess

        lhs = self.improve(lhs, rhs)

        return lhs

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
            if self._allow_max_iter:
                warnings.warn(
                    "maximum number of iterations {} exceeded".format(self._maxiter),
                    RuntimeWarning,
                )
            else:
                raise RuntimeError(
                    "maximum number of iterations {} exceeded".format(self._maxiter)
                )

        self.finish()

        return u

    def iterate(self, res):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def start(self):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def finish(self):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def set_init_guess(self, init_guess):
        self._init_guess = init_guess

    def get_init_guess(self):
        return self._init_guess

    def get_residual(self, lhs, rhs):
        return self._matrix @ lhs - rhs

    def get_matrix(self):
        return self._matrix

    def get_constraints(self):
        return self._cons
