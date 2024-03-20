import numpy as np
import scipy.sparse as spsp

from ..names import GlobNames as gn

from .solvermodule import SolverModule
from ..solver import Constraints
from ..util import Table, to_xtable
from ..util.proputils import optional_argument

TYPE = "type"

__all__ = ["LinsolveModule"]


class LinsolveModule(SolverModule):
    def init(self, globdat, **props):

        super().init(globdat, **props)

        # Get props
        solverprops = optional_argument(
            self, props, "solver", default={TYPE: "Cholmod"}
        )
        preconprops = optional_argument(self, props, "preconditioner", default={})
        self._get_mass_matrix = optional_argument(
            self, props, "getMassMatrix", default=False
        )
        self._get_strain_matrix = optional_argument(
            self, props, "getStrainMatrix", default=False
        )
        self._tnames = optional_argument(self, props, "tables", default=[])

        self._models = globdat[gn.MODELS]
        self._dc = globdat[gn.DOFSPACE].dof_count()

        solvertype = solverprops[TYPE]
        self._solver = globdat[gn.SOLVERFACTORY].get_solver(solvertype)
        self._solver.configure(globdat, **solverprops)

        self._precon = None

        if len(preconprops) > 0:
            precontype = preconprops[TYPE]
            self._precon = globdat[gn.PRECONFACTORY].get_precon(precontype)
            self._precon.configure(globdat, **preconprops)

    def solve(self, globdat):
        print("Running LinsolverModule")
        globdat[gn.TIMESTEP] = 1

        K = self.update_matrix(globdat)
        f = self.get_ext_vector(globdat)
        c = self.update_constraints(globdat)

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
            table = Table(size=nodecount)
            tbwts = np.zeros(nodecount)

            for model in self.get_relevant_models("GETTABLE", self._models):
                table, tbwts = model.GETTABLE(name, table, tbwts, globdat)

            to_xtable(table)

            for jcol in range(table.column_count()):
                values = table.get_col_values(None, jcol)
                table.set_col_values(None, jcol, values / tbwts)

            table.to_table()
            globdat[gn.TABLES][name] = table

        return "ok"

    def get_ext_vector(self, globdat):
        f_ext = np.zeros(globdat[gn.DOFSPACE].dof_count())

        for model in self.get_relevant_models("GETEXTFORCE", self._models):
            f_ext = model.GETEXTFORCE(f_ext, globdat)

        return f_ext

    def get_neumann_vector(self, globdat):
        f_neum = np.zeros(globdat[gn.DOFSPACE].dof_count())

        for model in self.get_relevant_models("GETNEUMANNFORCE", self._models):
            f_neum = model.GETNEUMANNFORCE(f_neum, globdat)

        return f_neum

    def update_matrix(self, globdat):
        K = self._get_empty_matrix(globdat)

        for model in self.get_relevant_models("GETMATRIX0", self._models):
            K = model.GETMATRIX0(K, globdat)

        return K

    def update_mass_matrix(self, globdat):
        M = self._get_empty_matrix(globdat)

        for model in self.get_relevant_models("GETMATRIX2", self._models):
            M = model.GETMATRIX2(M, globdat)

        return M

    def update_strain_matrix(self, globdat):
        B = self._get_empty_bmatrix(globdat)
        wts = np.zeros(B.shape[0])

        for model in self.get_relevant_models("GETMATRIXB", self._models):
            B, wts = model.GETMATRIXB(B, wts, globdat)

        # Divide non-zero entries by weights
        str_indices, dof_indices = B.nonzero()
        for i, str_idx in enumerate(str_indices):
            dof_idx = dof_indices[i]
            B[str_idx, dof_idx] /= wts[str_idx]

        return B

    def update_constraints(self, globdat):
        c = Constraints()

        for model in self.get_relevant_models("GETCONSTRAINTS", self._models):
            c = model.GETCONSTRAINTS(c, globdat)

        return c

    def advance(self, globat):
        pass

    def cancel(self, globdat):
        pass

    def commit(self, globdat):
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
