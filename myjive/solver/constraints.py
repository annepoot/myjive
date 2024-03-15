__all__ = ["Constraints"]


class Constraints:
    def __init__(self):
        self._ddofs = []
        self._dvals = []
        self._ndofs = []
        self._nvals = []

    def add_constraints(self, dofs, vals):
        for dof, val in zip(dofs, vals):
            self.add_constraint(dof, val)

    def add_constraint(self, dof, val):
        self.add_dirichlet(dof, val)

    def add_dirichlet(self, dof, val):
        self._ddofs.append(dof)
        self._dvals.append(val)

    def add_neumann(self, dof, val):
        self._ndofs.append(dof)
        self._nvals.append(val)

    def get_constraints(self):
        return self.get_dirichlet()

    def get_dirichlet(self):
        return self._ddofs, self._dvals

    def get_neumann(self):
        return self._ndofs, self._nvals
