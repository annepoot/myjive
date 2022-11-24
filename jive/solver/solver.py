import numpy as np

PRECISION = 'precision'
NOTIMPLEMENTEDMSG = 'this function needs to be implemented in an derived class'

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
        self._precision = 1e-5

    def configure(self, props):
        self._precision = props.get(PRECISION, self._precision)

    def start(self):
        pass

    def finish(self):
        pass

    def solve(self, rhs):
        lhs = np.zeros_like(rhs)
        return self.improve(lhs, rhs)

    def improve(self, lhs, rhs):
        return lhs

    def get_matrix(self):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def get_constraints(self):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)
