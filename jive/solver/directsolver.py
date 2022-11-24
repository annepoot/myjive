import scipy.sparse.linalg as spspla

from jive.solver.solver import Solver
from jive.solver.constrainer import Constrainer

class DirectSolver(Solver):

    def __init__(self):
        super().__init__()

        self._matrix = None
        self._cons = None
        self._conman = None

    def configure(self, props):
        super().configure(props)

    def update(self, matrix, constraints):
        self._cons = constraints
        self._conman = Constrainer(self._cons, matrix)
        self._matrix = self._conman.get_output_matrix()

    def solve(self, rhs):
        f = self._conman.get_rhs(rhs)
        u = spspla.spsolve(self._matrix, f)
        return u

    def get_matrix(self):
        return self._matrix

    def get_constraints(self):
        return self._cons


def declare(factory):
    factory.declare_solver('direct', DirectSolver)
