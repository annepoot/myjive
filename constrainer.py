import numpy as np


class Constrainer:
    def __init__(self):
        self._dofs = []
        self._vals = []

    def add_constraint(self, dof, val):
        self._dofs.append(dof)
        self._vals.append(val)

    def constrain(self, k, f):
        kc = np.copy(k)
        fc = np.copy(f)

        for dof, val in zip(self._dofs, self._vals):
            for i in range(kc.shape[0]):
                if i == dof:
                    fc[i] = val
                else:
                    fc[i] -= kc[i, dof] * val

            kc[:, dof] = kc[dof, :] = 0.0
            kc[dof, dof] = 1.0

        return kc, fc

    def get_constraints(self):
        return self._dofs, self._vals
