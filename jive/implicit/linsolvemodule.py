import numpy as np
import scipy.sparse as spsp

from jive.fem.names import GlobNames as gn
from jive.fem.names import ParamNames as pn
from jive.fem.names import Actions as act

import jive.util.proputils as pu
from jive.implicit.solvermodule import SolverModule
from jive.solver.constraints import Constraints
from jive.util.table import Table
from jive.util.xtable import to_xtable

SOLVER = "solver"
PRECONDITIONER = "preconditioner"
TYPE = "type"
GETMASSMATRIX = "getMassMatrix"
GETSTRAINMATRIX = "getStrainMatrix"
TABLES = "tables"

__all__ = ["LinsolveModule"]


class LinsolveModule(SolverModule):
    def init(self, props, globdat):
        super().init(props, globdat)

        myprops = props[self._name]
        self._get_mass_matrix = bool(eval(myprops.get(GETMASSMATRIX, "False")))
        self._get_strain_matrix = bool(eval(myprops.get(GETSTRAINMATRIX, "False")))
        self._tnames = pu.parse_list(myprops.get(TABLES, "[]"))

        self._model = globdat[gn.MODEL]
        self._dc = globdat[gn.DOFSPACE].dof_count()

        solverprops = myprops.get(SOLVER, {})
        solver = solverprops.get(TYPE, "Cholmod")
        self._solver = globdat[gn.SOLVERFACTORY].get_solver(solver)
        self._solver.configure(solverprops, globdat)

        preconprops = myprops.get(PRECONDITIONER, {})
        self._precon = None
        precon = preconprops.get(TYPE)
        if precon is not None:
            self._precon = globdat[gn.PRECONFACTORY].get_precon(precon)
            self._precon.configure(preconprops, globdat)

    def solve(self, globdat):
        model = globdat[gn.MODEL]

        print("Running LinsolverModule")
        globdat[gn.TIMESTEP] = 1

        K, _ = self.update_matrix(globdat)
        f = self.get_ext_vector(globdat)
        c = self.update_constraints(K, globdat)

        # Optionally get the mass matrix
        if self._get_mass_matrix:
            M = self.update_mass_matrix(globdat)

        # Optionally get the strain matrix
        if self._get_strain_matrix:
            B = self.update_strain_matrix(globdat)

        # Update the solver
        self._solver.update(K, c, self._precon)

        # Solve the system
        u = self._solver.solve(f)

        # Store rhs and solution in Globdat
        globdat[gn.EXTFORCE] = f
        globdat[gn.STATE0] = u

        # Store stiffness matrix and constraints in Globdat
        globdat[gn.MATRIX0] = K
        globdat[gn.CONSTRAINTS] = c

        # Optionally store mass matrix in Globdat
        if self._get_mass_matrix:
            globdat[gn.MATRIX2] = M

        # Optionally store strain matrix in Globdat
        if self._get_strain_matrix:
            globdat[gn.MATRIXB] = B

        # Compute stresses, strains, etc.
        if gn.TABLES not in globdat:
            globdat[gn.TABLES] = {}

        for name in self._tnames:
            nodecount = len(globdat[gn.NSET])

            params = {}
            params[pn.TABLE] = Table(size=nodecount)
            params[pn.TABLENAME] = name
            params[pn.TABLEWEIGHTS] = np.zeros(nodecount)

            model.take_action(act.GETTABLE, params, globdat)

            to_xtable(params[pn.TABLE])

            for jcol in range(params[pn.TABLE].column_count()):
                values = params[pn.TABLE].get_col_values(None, jcol)
                params[pn.TABLE].set_col_values(
                    None, jcol, values / params[pn.TABLEWEIGHTS]
                )

            params[pn.TABLE].to_table()
            globdat[gn.TABLES][name] = params[pn.TABLE]

        return "ok"

    def configure(self, props, globdat):
        pass

    def get_ext_vector(self, globdat):
        f_ext = np.zeros(self._dc)
        params = {pn.EXTFORCE: f_ext}

        self._model.take_action(act.GETEXTFORCE, params, globdat)

        return f_ext

    def get_neumann_vector(self, globdat):
        f_neum = np.zeros(self._dc)
        params = {pn.NEUMANNFORCE: f_neum}

        self._model.take_action(act.GETNEUMANNFORCE, params, globdat)

        return f_neum

    def update_matrix(self, globdat):
        params = {}
        params[pn.MATRIX0] = self._get_empty_matrix(globdat)
        params[pn.INTFORCE] = np.zeros(self._dc)

        self._model.take_action(act.GETMATRIX0, params, globdat)

        return params[pn.MATRIX0], params[pn.INTFORCE]

    def update_mass_matrix(self, globdat):
        params = {}
        params[pn.MATRIX2] = self._get_empty_matrix(globdat)

        self._model.take_action(act.GETMATRIX2, params, globdat)

        return params[pn.MATRIX2]

    def update_strain_matrix(self, globdat):
        params = {}
        params[pn.MATRIXB] = self._get_empty_bmatrix(globdat)
        params[pn.TABLEWEIGHTS] = np.zeros(params[pn.MATRIXB].shape[0])

        self._model.take_action(act.GETMATRIXB, params, globdat)

        # Divide non-zero entries by weights
        str_indices, dof_indices = params[pn.MATRIXB].nonzero()
        for i, str_idx in enumerate(str_indices):
            dof_idx = dof_indices[i]
            params[pn.MATRIXB][str_idx, dof_idx] /= params[pn.TABLEWEIGHTS][str_idx]

        return params[pn.MATRIXB]

    def update_constraints(self, K, globdat):
        c = Constraints()
        params = {pn.CONSTRAINTS: c}

        self._model.take_action(act.GETCONSTRAINTS, params, globdat)

        return c

    def advance(self, globat):
        # globdat[gn.MODEL].take_action(act.ADVANCE)
        pass

    def cancel(self, globdat):
        # globdat[gn.MODEL].take_action(act.CANCEL)
        pass

    def commit(self, globdat):
        # globdat[gn.MODEL].take_action(act.COMMIT)
        return True

    def _get_empty_matrix(self, globdat):
        rowindices = []
        colindices = []

        doftypes = globdat[gn.DOFSPACE].get_types()
        dc = globdat[gn.DOFSPACE].dof_count()

        for elem in globdat[gn.ESET]:
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, doftypes)

            for row in idofs:
                for col in idofs:
                    rowindices.append(row)
                    colindices.append(col)

        assert len(rowindices) == len(colindices)
        values = np.zeros(len(rowindices))

        K_empty = spsp.csr_array(
            (values, (rowindices, colindices)), shape=(dc, dc), dtype=float
        )
        return K_empty

    def _get_empty_bmatrix(self, globdat):
        rowindices = []
        colindices = []

        doftypes = globdat[gn.DOFSPACE].get_types()
        dc = globdat[gn.DOFSPACE].dof_count()
        nc = globdat[gn.NSET].size()
        rank = globdat[gn.DOFSPACE].type_count()
        strcount = rank * (rank + 1) // 2

        for elem in globdat[gn.ESET]:
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, doftypes)
            node_count = len(inodes)

            # Get the node index vector
            node_idx = np.zeros(node_count * strcount, dtype=int)
            for i in range(strcount):
                node_idx[i * node_count : (i + 1) * node_count] = inodes + nc * i

            for row in node_idx:
                for col in idofs:
                    rowindices.append(row)
                    colindices.append(col)

        assert len(rowindices) == len(colindices)
        values = np.zeros(len(rowindices))

        B_empty = spsp.csr_array(
            (values, (rowindices, colindices)), shape=(nc * strcount, dc), dtype=float
        )
        return B_empty
