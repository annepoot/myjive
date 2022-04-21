import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

from names import GlobNames as gn
from names import ParamNames as pn
from names import GPParamNames as gppn
from names import PropNames as prn
from names import Actions as act
from names import GPActions as gpact

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

        # Initialize model
        print('GaussianModule: Creating model...')
        m = modelfac.get_model(modelprops[prn.TYPE], self._modelname)
        m.configure(modelprops, globdat)
        globdat[self._modelname] = m

    def run(self, globdat):
        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[self._modelname]

        # Configure the model again, to make sure K, M and f are stored there as well
        model.configure()

        self._step += 1
        print('Running time step', self._step)
        globdat[gn.TIMESTEP] = self._step

        # Get the observation operators
        phi, phi_sub = gh.get_phis(N_obs=self._nobs, N_mesh=dc)

        # Store the observation operator in Globdat
        globdat['phi'] = phi

        # Get the relevant FEM results from Globdat
        K = globdat[gn.MATRIX0]
        M = globdat[gn.MATRIX2]
        u = globdat[gn.STATE0]
        f = globdat[gn.EXTFORCE]

        # Constrain K and f
        c = Constrainer()
        Kc, fc = c.constrain(K, f)

        u_prior = np.zeros(dc)
        u_post = np.zeros(dc)
        sigma_u_prior = np.zeros((dc, dc))
        sigma_u_post = np.zeros((dc, dc))

        u_params = {}
        u_params[gppn.FIELD] = 'u'
        u_params[gppn.FULLCOVARIANCE] = False
        u_params[gppn.PRIORMEAN] = u_prior
        u_params[gppn.POSTERIORMEAN] = u_post
        u_params[gppn.PRIORCOVARIANCE] = sigma_u_prior
        u_params[gppn.POSTERIORCOVARIANCE] = sigma_u_post

        f_prior = np.zeros(dc)
        f_post = np.zeros(dc)
        sigma_f_prior = np.zeros((dc, dc))
        sigma_f_post = np.zeros((dc, dc))

        f_params = {}
        f_params[gppn.FIELD] = 'f'
        f_params[gppn.FULLCOVARIANCE] = False
        f_params[gppn.PRIORMEAN] = f_prior
        f_params[gppn.POSTERIORMEAN] = f_post
        f_params[gppn.PRIORCOVARIANCE] = sigma_f_prior
        f_params[gppn.POSTERIORCOVARIANCE] = sigma_f_post

        model.take_action(gpact.GETPRIORMEAN, u_params, globdat)
        model.take_action(gpact.GETPOSTERIORMEAN, u_params, globdat)
        model.take_action(gpact.GETPRIORCOVARIANCE, u_params, globdat)
        model.take_action(gpact.GETPOSTERIORCOVARIANCE, u_params, globdat)

        model.take_action(gpact.GETPRIORMEAN, f_params, globdat)
        model.take_action(gpact.GETPOSTERIORMEAN, f_params, globdat)
        model.take_action(gpact.GETPRIORCOVARIANCE, f_params, globdat)
        model.take_action(gpact.GETPOSTERIORCOVARIANCE, f_params, globdat)

        # # Sparsify the stiffness matrix
        # smat = sparse.csr_matrix(Kc)

        # # Get the observed nodal forces
        # f_obs = phi.T @ fc

        # #############################
        # # POSTERIOR MEAN ON u AND f #
        # #############################

        # # Get the observation covariance matrix
        # obs_cov = self._alpha2 * phi.T @ M @ phi + np.identity(self._nobs) * self._noise2

        # # Use cholesky, because we need L1 later on again
        # L1 = np.linalg.cholesky(obs_cov)
        # v0 = np.linalg.solve(L1, f_obs)

        # # Get the posterior of both the nodal forces and nodal displacements
        # f_post = self._alpha2 * M @ phi @ np.linalg.solve(L1.T, v0)
        # u_post = linalg.spsolve(smat, f_post)

        # #############################
        # # POSTERIOR COVARIANCE ON f #
        # #############################

        # # Sigma = alpha2 * M - alpha4 * M * phi * inv(alpha2 * phi.T * M * phi + Sigma_e) * phi.T * M
        # V1 = np.linalg.solve(L1, phi.T @ M)
        # sigma_f_post = sigma_f_prior.copy()
        # for i in range(dc):
        #     sigma_f_post[i] -= self._alpha2**2 * V1[:,i].T @ V1[:,i]

        # #########################
        # # PRIOR COVARIANCE ON u #
        # #########################

        # # Sigma = alpha2 * inv(K) * M * inv(K)
        # L2 = np.linalg.cholesky(M)
        # V2 = np.linalg.solve(Kc, L2)
        # sigma_u_prior = np.zeros_like(sigma_f_prior)
        # for i in range(dc):
        #     sigma_u_prior[i] = self._alpha2 * V2[i,:] @ V2[i,:].T

        # #############################
        # # POSTERIOR COVARIANCE ON u #
        # #############################

        # # Sigma = alpha2 * inv(K) * M * inv(K) - alpha4 * inv(K) * M * phi * inv(alpha2 * phi.T * M * phi.T + Sigma_e) * phi.T * M * inv(K)
        # V3 = np.linalg.solve(Kc, V1.T)
        # sigma_u_post = sigma_u_prior.copy()
        # for i in range(dc):
        #     sigma_u_post[i] -= self._alpha2**2 * V3[i,:] @ V3[i,:].T

        # Optionally store stiffness matrix in Globdat
        if ( self._store_matrix ):
            globdat['phi'] = phi
            globdat['phi_sub'] = phi_sub
            globdat['f_post'] = f_params[gppn.POSTERIORMEAN]
            globdat['u_post'] = u_params[gppn.POSTERIORMEAN]
            globdat['sigma_f_prior'] = f_params[gppn.PRIORCOVARIANCE]
            globdat['sigma_f_post'] = f_params[gppn.POSTERIORCOVARIANCE]
            globdat['sigma_u_prior'] = u_params[gppn.PRIORCOVARIANCE]
            globdat['sigma_u_post'] = u_params[gppn.POSTERIORCOVARIANCE]

        return 'exit'

    def shutdown(self, globdat):
        pass

    def __solve(self, globdat):
        pass


def declare(factory):
    factory.declare_module('Gaussian', GaussianModule)
