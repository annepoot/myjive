import numpy as np
from scipy.linalg import solve, solve_triangular

from names import Actions as act
from names import GPActions as gpact
from names import ParamNames as pn
from names import GPParamNames as gppn
from names import GlobNames as gn
from names import PropNames as prn
from gpmodel import GPModel

NOBS = 'nobs'
OBSNOISE = 'obsNoise'
ALPHA = 'alpha'
BETA = 'beta'
PRIOR = 'prior'
TYPE = 'type'
FUNC = 'func'
HYPERPARAMS = 'hyperparams'
RANDOMOBS = 'randomObs'
SHAPE = 'shape'
INTSCHEME = 'intScheme'
PDNOISE = 'pdNoise'
ENSEMBLE = 'ensemble'


class GPEnKFModel(GPModel):

    def configure(self, props, globdat):

        super().configure(props, globdat)

        self._nens = props.get(ENSEMBLE, 100)


    def _configure_prior(self, params, globdat):

        # !!! Note: here, Sigma is still being assembled.
        # How to get the ensemble without sampling from the prior?
        super()._configure_prior(params, globdat)

        # u = m + sqrt(Sigma) * z

        # Get the relevant matrices
        # self._get_sqrtSigma()
        sqrtSigma = 0.1 * np.diag(np.sqrt(self._Mc.diagonal()))

        rng = params.get(gppn.RNG, np.random.default_rng())

        # Define the array that stores the samples
        self._X = np.zeros((self._dc, self._nens))

        # Get a loop to create all samples
        for i in range(self._nens):

            # Get a sample from the standard normal distribution
            z = rng.standard_normal(self._dc)

            # Get the sample of the force field
            f = sqrtSigma @ z + self._m

            self._X[:,i] = solve(self._Kc, f)

        # Get the mean of X, and the deviation from the mean
        EX = np.mean(self._X, axis=1)
        self._A = self._X - np.tile(EX, (self._nens,1)).T

        # Compute the product of H and A (which is nobs x nens)
        self._HA = self._H @ self._A

        # Delete Sigma altogether, to make sure it's never used
        del self._Sigma


    def _get_prior_covariance(self, params, globdat):

        fullSigma = params.get(gppn.FULLCOVARIANCE, False)

        ####################
        # PRIOR COVARIANCE #
        ####################

        # Check if the full covariance matrix should be returned
        if fullSigma:

            # If so, compute the full covariance matrix
            Sigma_prior = self._A @ self._A.T / (self._nens - 1)
            params[gppn.PRIORCOVARIANCE] = Sigma_prior

        else:

            # If not, compute only the diagonal of the covariance matrix
            var_prior = np.zeros(self._dc)
            for i in range(self._dc):
                var_prior[i] = self._A[i,:] @ self._A[i,:].T / (self._nens - 1)

            params[gppn.PRIORCOVARIANCE] = var_prior


    def _get_Sigma_obs(self):

        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._HA @ self._HA.T / (self._nens-1) + np.identity(self._nobs) * self._noise2

    def _get_sqrtObs(self):

        if not '_sqrtObs' in vars(self):
            self._get_Sigma_obs()
            self._sqrtObs = np.linalg.cholesky(self._Sigma_obs)

    def _get_sqrtSigma(self):

        if not '_sqrtSigma' in vars(self):
            self._sqrtSigma = self._A / np.sqrt(self._nens-1)

    def _premul_Sigma(self, X):

        return self._A @ (self._A.T @ X) / (self._nens-1)

    def _postmul_Sigma(self, X):

        return (X @ self._A) @ self._A.T / (self._nens-1)


def declare(factory):
    factory.declare_model('GPEnKF', GPEnKFModel)
