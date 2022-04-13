import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

# from module import *
from module import Module
from constrainer import Constrainer
import gaussianhelper as gh

NSTEPS = 'nsteps'
NOBS = 'nobs'
STOREMATRIX = 'storeMatrix'
RANDOMOBS = 'randomObs'

class GaussianModule(Module):

    def init(self, props, globdat):
        self._step = 0
        myprops = props[self._name]
        self._nsteps = int(myprops.get(NSTEPS,1))
        self._store_matrix = bool(eval(myprops.get(STOREMATRIX,'False')))

        self._nobs = int(myprops.get(NOBS,globdat[gn.DOFSPACE].dof_count()))

    def run(self, globdat):
        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[gn.MODEL]

        self._step += 1
        print('Running time step', self._step)
        globdat[gn.TIMESTEP] = self._step

        K = np.zeros((dc, dc))
        M = np.zeros((dc, dc))
        f = np.zeros(dc)
        c = Constrainer()
        phi, phi_sub = gh.get_phis(N_obs=self._nobs, N_mesh=dc)

        params = {}
        params[pn.MATRIX0] = K
        params[pn.MATRIX2] = M
        params[pn.EXTFORCE] = f
        params[pn.CONSTRAINTS] = c

        # Assemble K
        model.take_action(act.GETMATRIX0, params, globdat)

        # Assmemble M
        model.take_action(act.GETMATRIX2, params, globdat)

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


        # Optionally store stiffness matrix in Globdat
        if ( self._store_matrix ):
          globdat[gn.MATRIX0] = K
          globdat[gn.MATRIX2] = M
          globdat['phi'] = phi
          globdat['phi_sub'] = phi_sub

        if self._step >= self._nsteps:
            return 'exit'
        else:
            return 'ok'

    def shutdown(self, globdat):
        pass

    def __solve(self, globdat):
        pass


def declare(factory):
    factory.declare_module('Gaussian', GaussianModule)
