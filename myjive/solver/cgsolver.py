from .iterativesolver import IterativeSolver

__all__ = ["CGSolver"]


class CGSolver(IterativeSolver):
    def __init__(self):
        super().__init__()

        self._p = None

    def start(self):
        self._p = None

    def finish(self):
        self._p = None

    def iterate(self, res):
        r = -res

        if self._precon is None:
            z = r
        else:
            z = self._precon.solve(r)

        if self._p is None:
            self._p = z

        alpha = (r @ z) / (self._p @ self._matrix @ self._p)
        du = alpha * self._p

        r_new = r - alpha * self._matrix @ self._p

        if self._precon is None:
            z_new = r_new
        else:
            z_new = self._precon.solve(r_new)

        beta = (r_new @ z_new) / (r @ z)
        self._p = z_new + beta * self._p

        return du
