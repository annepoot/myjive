import numpy as np
import scipy.sparse as spsp

__all__ = ["Constrainer"]


class Constrainer:
    def __init__(self, constraints, inputmatrix):
        self.update(constraints, inputmatrix)

    def get_lhs(self, u):
        uc = u.copy()

        for dof, val in zip(*self._cons.get_constraints()):
            uc[dof] = val

        return uc

    def get_rhs(self, f):
        fc = f.copy()
        fc += self._rhs

        for dof, val in zip(*self._cons.get_constraints()):
            fc[dof] = val

        return fc

    def get_input_matrix(self):
        return self._input

    def get_output_matrix(self):
        return self._output

    def update(self, constraints, inputmatrix):
        self._cons = constraints
        self._input = inputmatrix
        self._output = spsp.csr_array(self._input.copy())
        self._rhs = np.zeros(self._output.shape[0])

        dofs, vals = self._cons.get_constraints()

        self._rhs -= self._output[:, dofs] @ vals
        self._rhs[dofs] = vals

        self._output[:, dofs] *= 0.0
        self._output[dofs, :] *= 0.0
        self._output[dofs, dofs] = 1.0

    def constrain(self, k, f):
        return self.apply_dirichlet(k, f)

    def apply_dirichlet(self, k, f):
        assert k is self._input

        return self.get_output_matrix(), self.get_rhs(f)

    def apply_neumann(self, f):
        fc = f.copy()

        for dof, val in zip(*self._cons.get_neumann()):
            fc[dof] += val

        return fc

    def new_constrainer(self, inputmatrix):
        c_new = Constrainer(inputmatrix)
        c_new.add_constraints(*self.get_constraints())
        return c_new
