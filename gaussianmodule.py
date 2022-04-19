import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

from names import GlobNames as gn
from names import ParamNames as pn
from names import PropNames as prn
from names import Actions as act

from module import Module
from constrainer import Constrainer
import gaussianhelper as gh
import testutils as tu
import matplotlib.pyplot as plt

NOBS = 'nobs'
OBSNOISE = 'obsNoise'
ALPHA = 'alpha'
STOREMATRIX = 'storeMatrix'
RANDOMOBS = 'randomObs'

class GaussianModule(Module):

    def init(self, props, globdat):
        self._step = 0
        myprops = props[self._name]
        self._store_matrix = bool(eval(myprops.get(STOREMATRIX,'False')))
        self._dc = globdat[gn.DOFSPACE].dof_count()

        # Get the appropriate model for this module
        self._modelname = myprops.get(gn.MODEL, gn.MODEL)
        modelprops = props[self._modelname]
        modelfac = globdat[gn.MODELFACTORY]

        # Get the number of observations either as a ratio or as a fixed integer
        self._nobs = float(myprops.get(NOBS,self._dc))
        if self._nobs < 1:
            self._nobs = self._dc * self._nobs
        self._nobs = int(np.round(self._nobs))

        self._noise2 = float(myprops.get(OBSNOISE))**2

        # Get the alpha parameter, or use the optimal alpha
        self._alpha2 = myprops.get(ALPHA, 'opt')
        if self._alpha2.isnumeric():
            self._alpha2 = float(self._alpha2)**2

        # Initialize model
        print('GaussianModule: Creating model...')
        m = modelfac.get_model(modelprops[prn.TYPE], self._modelname)
        m.configure(modelprops, globdat)
        globdat[self._modelname] = m

    def run(self, globdat):
        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[self._modelname]

        self._step += 1
        print('Running time step', self._step)
        globdat[gn.TIMESTEP] = self._step

        # Get the observation operators
        phi, phi_sub = gh.get_phis(N_obs=self._nobs, N_mesh=dc)

        # Get the relevant FEM results from Globdat
        K = globdat[gn.MATRIX0]
        M = globdat[gn.MATRIX2]
        u = globdat[gn.STATE0]
        f = globdat[gn.EXTFORCE]

        # Constrain K and f
        c = Constrainer()
        Kc, fc = c.constrain(K, f)

        # Sparsify the stiffness matrix
        smat = sparse.csr_matrix(Kc)

        # Get the observed nodal forces
        f_obs = phi.T @ fc

        #################
        # OPTIMAL ALPHA #
        #################

        # Check if alpha should be optimized
        if self._alpha2 == 'opt':

            # If so, determine the optimal value of alpha
            L = np.linalg.cholesky(phi.T @ M @ phi)
            v = np.linalg.solve(L, f_obs)
            self._alpha2 = v.T @ v / f_obs.shape[0]

            print('Setting alpha to optimal value: {:.4f}'.format(self._alpha2))

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
        V2 = np.linalg.solve(Kc, L2)
        sigma_u_prior = np.zeros_like(sigma_f_prior)
        for i in range(dc):
            sigma_u_prior[i] = self._alpha2 * V2[i,:] @ V2[i,:].T

        #############################
        # POSTERIOR COVARIANCE ON u #
        #############################

        # Sigma = alpha2 * inv(K) * M * inv(K) - alpha4 * inv(K) * M * phi * inv(alpha2 * phi.T * M * phi.T + Sigma_e) * phi.T * M * inv(K)
        V3 = np.linalg.solve(Kc, V1.T)
        sigma_u_post = sigma_u_prior.copy()
        for i in range(dc):
            sigma_u_post[i] -= self._alpha2**2 * V3[i,:] @ V3[i,:].T

        # Optionally store stiffness matrix in Globdat
        if ( self._store_matrix ):
            globdat['phi'] = phi
            globdat['phi_sub'] = phi_sub
            globdat['f_post'] = f_post
            globdat['u_post'] = u_post
            globdat['sigma_f_prior'] = sigma_f_prior
            globdat['sigma_f_post'] = sigma_f_post
            globdat['sigma_u_prior'] = sigma_u_prior
            globdat['sigma_u_post'] = sigma_u_post

        return 'exit'

    def shutdown(self, globdat):
        pass

    def __solve(self, globdat):
        pass


def declare(factory):
    factory.declare_module('Gaussian', GaussianModule)
