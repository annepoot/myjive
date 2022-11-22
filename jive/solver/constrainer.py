import numpy as np
import scipy.sparse as spsp

class Constrainer:
    def __init__(self, inputmatrix):
        self._ddofs = []
        self._dvals = []
        self._ndofs = []
        self._nvals = []

        self._input = inputmatrix
        self._output = spsp.csr_array(self._input.copy())

        self._rhs = np.zeros(self._output.shape[0])

    def add_constraints(self, dofs, vals):
        for dof, val in zip(dofs, vals):
            self.add_constraint(dof, val)

    def add_constraint(self, dof, val):
        self.add_dirichlet(dof, val)

    def add_dirichlet(self, dof, val):
        self._ddofs.append(dof)
        self._dvals.append(val)

        for i in range(self._output.shape[0]):
            if i == dof:
                self._rhs[i] = val
            else:
                self._rhs[i] -= self._output[i, dof] * val

        self._output[:,[dof]] *= 0.0
        self._output[[dof], :] *= 0.0
        self._output[dof, dof] = 1.0

    def add_neumann(self, dof, val):
        self._ndofs.append(dof)
        self._nvals.append(val)

    def get_lhs(self, u):
        uc = u.copy()

        if len(self._ddofs) == 0:
            return u
        else:
            for dof, val in zip(self._ddofs, self._dvals):
                uc[dof] = val

        return uc

    def get_rhs(self, f):
        fc = f.copy()
        fc += self._rhs

        for dof, val in zip(self._ddofs, self._dvals):
            fc[dof] = val

        return fc

    def get_input_matrix(self):
        return self._input

    def get_output_matrix(self):
        return self._output

    def constrain(self, k, f):
        assert k is self._input

        return self.get_output_matrix(), self.get_rhs(f)

    def apply_dirichlet(self, k, f):
        # assert k is self._input

        return self.get_output_matrix(), self.get_rhs(f)

    def apply_neumann(self, f):
        fc = f.copy()

        for dof, val in zip(self._ndofs, self._nvals):
            fc[dof] += val

        return fc

    def get_constraints(self):
        return self.get_dirichlet()

    def get_dirichlet(self):
        return self._ddofs, self._dvals

    def get_neumann(self):
        return self._ndofs, self._nvals

    def new_constrainer(self, inputmatrix):
        c_new = Constrainer(inputmatrix)
        c_new.add_constraints(*self.get_constraints())
        return c_new
