import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spspla

from jive.fem.names import GlobNames as gn
from jive.fem.names import ParamNames as pn
from jive.fem.names import Actions as act

import jive.util.proputils as pu
from jive.implicit.solvermodule import SolverModule
from jive.solver.constrainer import Constrainer
from jive.util.table import Table

STOREMATRIX = 'storeMatrix'
STORECONSTRAINTS = 'storeConstraints'
GETMASSMATRIX = 'getMassMatrix'
GETUNITMASSMATRIX = 'getUnitMassMatrix'
TABLES = 'tables'


class LinsolveModule(SolverModule):

    def init(self, props, globdat):

        myprops = props[self._name]
        self._store_matrix = bool(eval(myprops.get(STOREMATRIX,'False')))
        self._store_constraints = bool(eval(myprops.get(STORECONSTRAINTS,'False')))
        self._get_mass_matrix = bool(eval(myprops.get(GETMASSMATRIX,'False')))
        self._tnames = pu.parse_list(myprops.get(TABLES, '[]'))

        self._model = globdat[gn.MODEL]
        self._dc = globdat[gn.DOFSPACE].dof_count()

    def solve(self, globdat):

        model = globdat[gn.MODEL]

        print('Running LinsolverModule')
        globdat[gn.TIMESTEP] = 1

        K, _ = self.update_matrix(globdat)
        f = self.get_ext_vector(globdat)
        c = self.update_constraints(globdat)

        # Constrain K and f
        Kc, fc = c.constrain(K, f)

        # Optionally get the mass matrix
        if self._get_mass_matrix:
            M = self._get_empty_matrix(globdat)
            params = {pn.MATRIX2: M}
            model.take_action(act.GETMATRIX2, params, globdat)

        # Solve the system
        u = spspla.spsolve(Kc, fc)

        # Store rhs and solution in Globdat
        globdat[gn.EXTFORCE] = f
        globdat[gn.STATE0] = u

        # Optionally store stiffness matrix in Globdat
        if self._store_matrix:
            globdat[gn.MATRIX0] = K

            # Optionally store mass matrix in Globdat
            if self._get_mass_matrix:
                globdat[gn.MATRIX2] = M

        # Optionally store the constrainer in Globdat
        if self._store_constraints:
            globdat[gn.CONSTRAINTS] = c

        # Compute stresses, strains, etc.
        globdat[gn.TABLES] = {}

        for name in self._tnames:
            params = {}
            params[pn.TABLE] = Table()
            params[pn.TABLENAME] = name
            params[pn.TABLEWEIGHTS] = np.zeros(len(globdat[gn.NSET]))

            model.take_action(act.GETTABLE, params, globdat)

            globdat[gn.TABLES][name] = params[pn.TABLE]

        return 'ok'

    def configure(self, props, globdat):
        pass

    def get_ext_vector(self, globdat):
        f_ext = np.zeros(self._dc)
        params = {pn.EXTFORCE: f_ext}

        self._model.take_action(act.GETEXTFORCE, params, globdat)

        return f_ext

    def update_matrix(self, globdat):
        K = self._get_empty_matrix(globdat)
        f_int = np.zeros(self._dc)
        params = {pn.MATRIX0: K, pn.INTFORCE:f_int}

        self._model.take_action(act.GETMATRIX0, params, globdat)

        return K, f_int

    def update_constraints(self, globdat):
        c = Constrainer()
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

        K_empty = spsp.csr_array((values, (rowindices, colindices)), shape=(dc,dc), dtype=float)
        return K_empty


def declare(factory):
    factory.declare_module('Linsolve', LinsolveModule)
