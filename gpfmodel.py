import numpy as np
import scipy.sparse as spsp

from jive.solver.numba.cholesky import sparse_cholesky
from jive.solver.numba.spsolve import solve_triangular
from jive.fem.names import GPParamNames as gppn
from jive.fem.names import GlobNames as gn
from gpmodel import GPModel

PRIOR = 'prior'
DIAGONALIZED = 'diagonalized'


class GPfModel(GPModel):

    def configure(self, props, globdat):

        super().configure(props, globdat)

        self._diagonalized = props[PRIOR].get(DIAGONALIZED, False)

    def _configure_prior(self, params, globdat):

        # Get the observation operator
        self._Phic = spsp.csr_array(self._Phic)
        self._H = self._Phic.T
        self._H = self._H.tocsr()

        # Define the mean in terms of the force vector as well
        self._m = self._mf

        # Get the centered observation vector (as a deviation from the prior mean)
        self._y = self._g - self._Phic.T @ self._m

        # Check if alpha or beta should be optimized
        if len(self._hyperparams) == 1:
            key, value = list(self._hyperparams.items())[0]

            if value == 'opt':

                if key == 'alpha':
                    Sigma = globdat[gn.MATRIX2]
                elif key == 'beta':
                    Sigma = globdat[gn.MATRIX0]
                else:
                    raise ValueError('cannot find optimal value for ' + key)

                # Apply boundary conditions to the prior
                Sigma = self._apply_covariance_bcs(Sigma)

                self._hyperparams[key] = self._get_param_opt(Sigma)

        super()._configure_prior(params, globdat)


    def _get_prior_mean(self, params, globdat):

        # Get the prior mean on f
        super()._get_prior_mean(params, globdat)

        # Inform GPSolverModule that force-related info is returned
        params[gppn.FIELD] = 'f'

    def _get_posterior_mean(self, params, globdat):

        # Get the posterior mean on f
        super()._get_posterior_mean(params, globdat)

        # Inform GPSolverModule that force-related info is returned
        params[gppn.FIELD] = 'f'

    def _get_prior_covariance(self, params, globdat):

        # Get the prior covariance on f first
        super()._get_prior_covariance(params, globdat)

        # Inform GPSolverModule that displacement-related info is returned
        params[gppn.FIELD] = 'f'

    def _get_posterior_covariance(self, params, globdat):

        # Get the posterior covariance on f first
        super()._get_posterior_covariance(params, globdat)

        # Inform GPSolverModule that displacement-related info is returned
        params[gppn.FIELD] = 'f'

    def _get_prior_samples(self, params, globdat):

        # Get the prior samples on f
        super()._get_prior_samples(params, globdat)

        # Inform GPSolverModule that force-related info is returned
        params[gppn.FIELD] = 'f'

    def _get_posterior_samples(self, params, globdat):

        # Get the posterior samples on f
        super()._get_posterior_samples(params, globdat)

        # Inform GPSolverModule that force-related info is returned
        params[gppn.FIELD] = 'f'

    def _kalman_update(self, params, globdat):

        # Perform the Kalman update on f
        super()._kalman_update(params, globdat)

        # Inform GPSolverModule that force-related info is returned
        params[gppn.FIELD] = 'f'

    def _get_param_opt(self, Sigma_fc):

        # Determine the optimal value of alpha
        L = sparse_cholesky(self._H @ Sigma_fc @ self._H.T)
        v = self._solve_triangular(L, self._y, lower=True)
        alpha2 = v.T @ v / self._nobs

        return np.sqrt(alpha2)

    def _apply_covariance_bcs(self, Sigma):
        Sigmac = Sigma.copy()

        # Set the covariance of the DBCs to 0
        Sigmac[self._cdofs,:] *= 0.0
        Sigmac[:,self._cdofs] *= 0.0

        # Add a tiny noise to ensure Sigma is positive definite rather than semidefinite
        Sigmac += self._pdnoise2 * spsp.identity(self._dc)

        return Sigmac

    def _solve_triangular(self, A, b, lower):
        return solve_triangular(A.tocsr(), b, lower=lower)

    def _get_Sigma_obs(self):

        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._H @ self._Sigma @ self._H.T + spsp.identity(self._nobs) * self._noise2

    def _get_sqrtObs(self):

        if not '_sqrtObs' in vars(self):
            self._get_Sigma_obs()
            self._sqrtObs = sparse_cholesky(self._Sigma_obs)

    def _get_sqrtSigma(self):

        if not '_sqrtSigma' in vars(self):
            if self._diagonalized:
                self._sqrtSigma = spsp.diags(np.sqrt(self._Sigma.sum(axis=1)), format='csr')
            else:
                self._sqrtSigma = sparse_cholesky(self._Sigma)

    def _get_sqrtNoise(self):

        if not '_sqrtNoise' in vars(self):
            self._sqrtNoise = np.sqrt(self._noise2) * spsp.identity(self._nobs)



def declare(factory):
    factory.declare_model('GPf', GPfModel)
