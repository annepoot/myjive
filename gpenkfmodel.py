import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spspla

from jive.solver.jit.cholesky import sparse_cholesky
from jive.fem.names import GPActions as gpact
from jive.fem.names import GPParamNames as gppn
from gpmodel import GPModel

ENSEMBLE = 'ensemble'
SEED = 'seed'
PRIOR = 'prior'
PREMULTIPLIER = 'premultiplier'
DIAGONALIZED = 'diagonalized'


class GPEnKFModel(GPModel):

    def configure(self, props, globdat):

        super().configure(props, globdat)

        self._nens = int(props.get(ENSEMBLE, 100))
        self._seed = eval(props.get(SEED,'None'))
        if not self._seed is None:
            self._seed = int(self._seed)

        self._premultiplier = props[PRIOR].get(PREMULTIPLIER, None)
        self._diagonalized = props[PRIOR].get(DIAGONALIZED, False)

    def _configure_prior(self, params, globdat):

        super()._configure_prior(params, globdat)

        # Define a dictionary with relevant functions
        eval_dict = self._get_eval_dict()

        # Get the premultiplier matrix (if necessary)
        if not self._premultiplier is None:
            preK = eval(self._premultiplier, eval_dict)

        # Check if we have a diagonalizable prior
        if self._diagonalized:

            # If so, get the cholesky root from the diagonal
            self._sqrtSigma = spsp.diags(np.sqrt(self._Sigma.sum(axis=1)), format='csr')

        else:

            # If not, do a full Cholesky decomposition (assuming Sigma is a sparse matrix)
            self._sqrtSigma = sparse_cholesky(self._Sigma)

        # Set the params for the ensemble samples
        ensemble_params = {gppn.NSAMPLE:self._nens, gppn.RNG:np.random.default_rng(self._seed)}

        # Take the GETPRIORSAMPLES action to build the ensemble
        self.take_action(gpact.GETPRIORSAMPLES, ensemble_params, globdat)

        # Obtain the ensemble
        self._X = ensemble_params[gppn.PRIORSAMPLES]

        # Apply the premultiplier to the samples
        for i in range(self._nens):
            self._X[:,i] = spspla.spsolve(preK, self._X[:,i])

        # Now, delete self._Sigma and self._sqrtSigma. We don't need them any more!
        del self._Sigma
        del self._sqrtSigma

        # Get the mean of X, and the deviation from the mean
        EX = np.mean(self._X, axis=1)
        self._A = self._X - np.tile(EX, (self._nens,1)).T

        # Compute the product of H and A (which is nobs x nens)
        self._HA = self._H @ self._A

    def _get_prior_covariance(self, params, globdat):

        ####################
        # PRIOR COVARIANCE #
        ####################

        # Compute the full covariance matrix
        Sigma_prior = self._A @ self._A.T / (self._nens - 1)
        params[gppn.PRIORCOVARIANCE] = Sigma_prior

    def _apply_covariance_bcs(self, Sigma):
        Sigmac = Sigma.copy()

        # Set the covariance of the DBCs to 0
        Sigmac[self._cdofs,:] *= 0.0
        Sigmac[:,self._cdofs] *= 0.0

        # Add a tiny noise to ensure Sigma is positive definite rather than semidefinite
        Sigmac += self._pdnoise2 * spsp.identity(self._dc)

        return Sigmac

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
