import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

from module import *
from constrainer import Constrainer

NSTEPS = 'nsteps'


class SolverModule(Module):

    def init(self, props, globdat):
        myprops = props[self._name]

        self._nsteps = int(myprops[NSTEPS])
        self._step = 0

    def run(self, globdat):
        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[gn.MODEL]

        self._step += 1
        print('Running time step', self._step)
        globdat[gn.TIMESTEP] = self._step

        K = np.zeros((dc, dc))
        f = np.zeros(dc)
        c = Constrainer()

        params = {}
        params[pn.MATRIX0] = K
        params[pn.EXTFORCE] = f
        params[pn.CONSTRAINTS] = c

        # Assemble K
        model.take_action(act.GETMATRIX0, params, globdat)

        # Assemble f
        model.take_action(act.GETEXTFORCE, params, globdat)

        # Get constraints
        model.take_action(act.GETCONSTRAINTS, params, globdat)

        # Constrain K and f
        Kc, fc = c.constrain(K, f)

        # Sparsify and solve
        smat = sparse.csr_matrix(Kc)
        u = linalg.spsolve(smat, fc)

        # Store solution in Globdat
        globdat[gn.STATE0] = u

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
