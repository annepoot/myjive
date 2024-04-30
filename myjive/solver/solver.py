from ..util import saveconfig as sg
import numpy as np

NOTIMPLEMENTEDMSG = "this function needs to be implemented in an derived class"

__all__ = ["Solver", "SolverFactory"]


class SolverFactory:
    def __init__(self):
        self._solvers = {}

    def declare_solver(self, typ, solver):
        self._solvers[typ] = solver

    def get_solver(self, typ, name):
        solver = self._solvers.get(typ)
        if not solver:
            raise ValueError(typ)
        return solver(name)

    def is_solver(self, typ):
        return typ in self._solvers


class Solver:
    def __init__(self, name):
        self._name = name
        self._config = {}
        self.precon_mode = False

    @classmethod
    def get_type(cls):
        typ = cls.__name__
        if typ[-6:] == "Solver":
            typ = typ[:-6]
        return typ

    @classmethod
    def declare(cls, factory):
        typ = cls.get_type()
        factory.declare_solver(typ, cls)

    @sg.save_config
    def configure(self, globdat, *, precision=1e-8):
        self._precision = precision

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

    def get_config(self):
        if len(self._config) == 0:
            raise NotImplementedError("Empty solver get_config")
        else:
            return self._config

    def get_name(self):
        return self._name

    @staticmethod
    def save_config(configure):
        return sg.save_config(configure)
