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
OBSNOISE = 'obsNoise'
ALPHA = 'alpha'
STOREMATRIX = 'storeMatrix'
RANDOMOBS = 'randomObs'

class GaussianModule(Module):

    def init(self, props, globdat):
        self._step = 0
        myprops = props[self._name]
        self._nsteps = int(myprops.get(NSTEPS,1))
        self._store_matrix = bool(eval(myprops.get(STOREMATRIX,'False')))

        self._nobs = int(myprops.get(NOBS,globdat[gn.DOFSPACE].dof_count()))
        self._noise2 = float(myprops.get(OBSNOISE))**2
        self._alpha2 = float(myprops.get(ALPHA))**2

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

        # Get the observed nodal forces
        f_obs = phi.T @ fc

        #############################
        # POSTERIOR MEAN ON u AND f #
        #############################

        # Get the observation covariance matrix
        obs_cov = self._alpha2 * phi.T @ M @ phi + np.identity(self._nobs) * self._noise2

        # Use cholesky, because we need L1 later on again
        L1 = np.linalg.cholesky(obs_cov)
        v0 = np.linalg.solve(L1, f_obs)

        # Get the posterior of both the nodal forces and nodal displacements
        f_post = self._alpha2 * M @ phi @ np.linalg.solve(L1.T, v0)
        u_post = linalg.spsolve(smat, f_post)

        #########################
        # PRIOR COVARIANCE ON f #
        #########################

        # Sigma = alpha2 * M
        sigma_f_prior = self._alpha2 * M.diagonal()

        #############################
        # POSTERIOR COVARIANCE ON f #
        #############################

        # Sigma = alpha2 * M - alpha4 * M * phi * inv(alpha2 * phi.T * M * phi + Sigma_e) * phi.T * M
        V1 = np.linalg.solve(L1, phi.T @ M)
        sigma_f_post = sigma_f_prior.copy()
        for i in range(dc):
            sigma_f_post[i] -= self._alpha2**2 * V1[:,i].T @ V1[:,i]

        #########################
        # PRIOR COVARIANCE ON u #
        #########################

        # Sigma = alpha2 * inv(K) * M * inv(K)
        L2 = np.linalg.cholesky(M)
        V2 = linalg.spsolve(smat, L2)
        sigma_u_prior = np.zeros_like(sigma_f_prior)
        for i in range(dc):
            sigma_u_prior[i] = self._alpha2 * V2[i,:] @ V2[i,:].T

        #############################
        # POSTERIOR COVARIANCE ON u #
        #############################

        # Sigma = alpha2 * inv(K) * M * inv(K) - alpha4 * inv(K) * M * phi * inv(alpha2 * phi.T * M * phi.T + Sigma_e) * phi.T * M * inv(K)
        V3 = linalg.spsolve(smat, V1.T)
        sigma_u_post = sigma_u_prior.copy()
        for i in range(dc):
            sigma_u_post[i] -= self._alpha2**2 * V3[i,:] @ V3[i,:].T

        # Optionally store stiffness matrix in Globdat
        if ( self._store_matrix ):
          globdat[gn.MATRIX0] = K
          globdat[gn.MATRIX2] = M
          globdat['phi'] = phi
          globdat['phi_sub'] = phi_sub
          globdat['u_post'] = u_post
          globdat['sigma_u_prior'] = sigma_u_prior
          globdat['sigma_u_post'] = sigma_u_post

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
