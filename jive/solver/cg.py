from jive.solver.iterativesolver import IterativeSolver

class CG(IterativeSolver):

    def __init__(self):
        super().__init__()

        self._beta = None
        self._p = None

    def iterate(self, res):
        r = -res

        if self._p is None:
            self._p = r

        alpha = (r @ r) / (self._p @ self._matrix @ self._p)
        du = alpha * self._p

        r_new = r - alpha * self._matrix @ self._p

        self._beta = (r_new @ r_new) / (r @ r)
        self._p = r_new + self._beta * self._p

        return du

def declare(factory):
    factory.declare_solver('CG', CG)
