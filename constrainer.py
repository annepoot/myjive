import numpy as np

class Constrainer:
    def __init__(self):
        self._ddofs = []
        self._dvals = []
        self._ndofs = []
        self._nvals = []

    def add_constraint(self, dof, val):
        self.add_dirichlet(dof, val)

    def add_dirichlet(self, dof, val):
        self._ddofs.append(dof)
        self._dvals.append(val)

    def add_neumann(self, dof, val):
        self._ndofs.append(dof)
        self._nvals.append(val)

    def constrain(self, k, f):
        kc = k.copy()
        fc = f.copy()

        for dof, val in zip(self._ddofs, self._dvals):
            for i in range(kc.shape[0]):
                if i == dof:
                    fc[i] = val
                else:
                    fc[i] -= kc[i, dof] * val

            kc[:, dof] = kc[dof, :] = 0.0
            kc[dof, dof] = 1.0

        return kc, fc

    def get_constraints(self):
        return self.get_dirichlet()

    def get_dirichlet(self):
        return self._ddofs, self._dvals

    def get_neumann(self):
        return self._ndofs, self._nvals
