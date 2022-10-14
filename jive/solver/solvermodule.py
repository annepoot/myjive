import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spspla

from jive.fem.names import GlobNames as gn
from jive.fem.names import ParamNames as pn
from jive.fem.names import Actions as act

import jive.util.proputils as pu
from jive.app.module import Module
from jive.solver.constrainer import Constrainer
from jive.util.table import Table

NSTEPS = 'nsteps'
STOREMATRIX = 'storeMatrix'
STORECONSTRAINTS = 'storeConstraints'
GETMASSMATRIX = 'getMassMatrix'
GETUNITMASSMATRIX = 'getUnitMassMatrix'
TABLES = 'tables'

class SolverModule(Module):

    def init(self, props, globdat):

        myprops = props[self._name]
        self._step = 0
        self._nsteps = int(myprops.get(NSTEPS,1))
        self._store_matrix = bool(eval(myprops.get(STOREMATRIX,'False')))
        self._store_constraints = bool(eval(myprops.get(STORECONSTRAINTS,'False')))
        self._get_mass_matrix = bool(eval(myprops.get(GETMASSMATRIX,'False')))
        self._tnames = pu.parse_list(myprops.get(TABLES, '[]'))

    def run(self, globdat):

        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[gn.MODEL]

        self._step += 1
        print('Running time step', self._step)
        globdat[gn.TIMESTEP] = self._step

        K = spsp.csr_array((dc, dc))
        f = np.zeros(dc)
        c = Constrainer()

        params = {pn.MATRIX0: K, pn.EXTFORCE: f, pn.CONSTRAINTS: c}

        # Assemble K
        model.take_action(act.GETMATRIX0, params, globdat)

        # Assemble f
        model.take_action(act.GETEXTFORCE, params, globdat)

        # Get constraints
        model.take_action(act.GETCONSTRAINTS, params, globdat)

        # Constrain K and f
        Kc, fc = c.constrain(K, f)

        # Optionally get the mass matrix
        if self._get_mass_matrix:
            M = spsp.csr_array((dc, dc))
            params[pn.MATRIX2] = M
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

        if self._step >= self._nsteps:
            return 'exit'
        else:
            return 'ok'

    def shutdown(self, globdat):
        pass

    def __solve(self, globdat):
        pass


def declare(factory):
    factory.declare_module('Solver', SolverModule)
