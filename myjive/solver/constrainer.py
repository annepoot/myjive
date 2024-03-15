import numpy as np
import scipy.sparse as spsp

__all__ = ["Constrainer"]


class Constrainer:
    def __init__(self, constraints, inputmatrix):
        self._cons = constraints
        self._input = inputmatrix
        self._output = spsp.csr_array(self._input.copy())
        self._rhs = np.zeros(self._output.shape[0])

        self.update()

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

    def update(self, constraints=None, inputmatrix=None):
        if constraints is not None:
            self._cons = constraints

        if inputmatrix is not None:
            self._input = inputmatrix

        for dof, val in zip(*self._cons.get_constraints()):
            for i in range(self._output.shape[0]):
                if i == dof:
                    self._rhs[i] = val
                else:
                    self._rhs[i] -= self._output[i, dof] * val

            self._output[:, [dof]] *= 0.0
            self._output[[dof], :] *= 0.0
            self._output[dof, dof] = 1.0

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
