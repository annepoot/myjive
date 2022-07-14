import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spspla

from names import GPParamNames as gppn
from gpmodel import GPModel

ENSEMBLE = 'ensemble'
PRIOR = 'prior'
PREMULTIPLIER = 'premultiplier'
DIAGONALIZED = 'diagonalized'


class GPEnKFModel(GPModel):

    def configure(self, props, globdat):

        super().configure(props, globdat)

        self._nens = int(props.get(ENSEMBLE, 100))
        self._premultiplier = props[PRIOR].get(PREMULTIPLIER, None)
        self._diagonalized = props[PRIOR].get(DIAGONALIZED, not self._premultiplier is None)


    def _configure_prior(self, params, globdat):

        # Define a dictionary with relevant functions
        eval_dict = {'inv':spspla.inv, 'exp':np.exp, 'norm':np.linalg.norm, 'np':np}
        eval_dict.update(self._hyperparams)

        # Add the mass and stiffness matrices to the dictionary
        eval_dict['M'] = self._Mc
        eval_dict['K'] = self._Kc

        # Get the covariance matrix
        Sigma = eval(self._covariance, eval_dict)

        # Set the covariance of the DBCs to 0
        Sigma[self._cdofs,:] = Sigma[:,self._cdofs] = 0.0

        # Add a tiny noise to ensure Sigma is positive definite rather than semidefinite
        Sigma += self._pdnoise2 * spsp.identity(self._dc)

        # Get the premultiplier matrix (if necessary)
        if not self._premultiplier is None:
            preK = eval(self._premultiplier, eval_dict)

        # Check if we have a diagonalizable prior
        if self._diagonalized:

            # If so, get the cholesky root from the diagonal
            sqrtSigma = spsp.diags(np.sqrt(Sigma.diagonal()), format='csr')

        else:

            # If not, do a full Cholesky decomposition
            sqrtSigma = np.linalg.cholesky(Sigma)

        rng = params.get(gppn.RNG, np.random.default_rng())

        # Define the array that stores the samples
        self._X = np.zeros((self._dc, self._nens))

        # Get a loop to create all samples
        for i in range(self._nens):

            # Get a sample from the standard normal distribution
            z = rng.standard_normal(self._dc)

            # Get the sample of the force field
            u = sqrtSigma @ z + self._m

            # Apply the premultiplier matrix to each sample (if necessary)
            if not self._premultiplier is None:
                u = spspla.spsolve(preK, u)

            # Add the sample to the ensemble
            self._X[:,i] = u

        # Get the mean of X, and the deviation from the mean
        EX = np.mean(self._X, axis=1)
        self._A = self._X - np.tile(EX, (self._nens,1)).T

        # Compute the product of H and A (which is nobs x nens)
        self._HA = self._H @ self._A


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
