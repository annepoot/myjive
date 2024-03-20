from myjive.util.proputils import optional_argument
import numpy as np

NOTIMPLEMENTEDMSG = "this function needs to be implemented in an derived class"

__all__ = ["Solver", "SolverFactory"]


class SolverFactory:
    def __init__(self):
        self._solvers = {}

    def declare_solver(self, typ, solver):
        self._solvers[typ] = solver

    def get_solver(self, typ):
        solver = self._solvers.get(typ)
        if not solver:
            raise ValueError(typ)
        return solver()

    def is_solver(self, typ):
        return typ in self._solvers


class Solver:
    def __init__(self):
        self.precon_mode = False

    @classmethod
    def declare(cls, factory):
        name = cls.__name__
        if len(name) > 6 and name[-6:] == "Solver":
            name = name[:-6]
        factory.declare_solver(name, cls)

    def configure(self, globdat, **props):
        self._precision = optional_argument(self, props, "precision", default=1e-8)

    def start(self):
        pass

    def finish(self):
        pass

    def solve(self, rhs):
        lhs = np.zeros_like(rhs)

        lhs = self.improve(lhs, rhs)

        return lhs

    def multisolve(self, rhs):
        if hasattr(rhs, "toarray"):
            rhs_mat = rhs.toarray()
        else:
            rhs_mat = rhs

        lhs = np.zeros_like(rhs_mat)

        for j in range(rhs_mat.shape[1]):
            lhs[:, j] = self.solve(rhs_mat[:, j])

        return lhs

    def improve(self, lhs, rhs):
        return NotImplementedError(NOTIMPLEMENTEDMSG)

    def get_matrix(self):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def get_constraints(self):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)
