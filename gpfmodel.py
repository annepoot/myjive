import numpy as np
import scipy.sparse as spsp
import scipy.linalg as spla
import scipy.sparse.linalg as spspla

from jive.fem.names import GPActions as gpact
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
        self._H = self._Phic.T

        # Define the mean in terms of the force vector as well
        self._m = self._mf

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

        # Get the relevant matrices
        self._get_sqrtSigma()

        # Solve the system for each dof
        if not '_V2' in vars(self):
            self._V2 = spspla.spsolve(self._Kc, self._sqrtSigma)

        # Convert to dense if necessary
        if hasattr(self._V2, 'todense'):
            self._V2 = self._V2.todense()

        # Compute the full covariance matrix
        Sigma_prior = self._V2 @ self._V2.T
        params[gppn.PRIORCOVARIANCE] = Sigma_prior

        # Inform GPSolverModule that displacement-related info is returned
        params[gppn.FIELD] = 'u'

    def _get_posterior_covariance(self, params, globdat):

        # Get the posterior covariance on f first
        super()._get_posterior_covariance(params, globdat)

        # Check if the prior variance is given in the params, otherwise, take an action to obtain it
        if not gppn.PRIORCOVARIANCE in params:
            self.take_action(gpact.GETPRIORCOVARIANCE, params, globdat)

        # Get the prior variance from the params
        var_prior = params[gppn.PRIORCOVARIANCE]

        # Get the relevant matrices
        self._get_V1()

        # Solve the system for each coarse dof
        if not '_V3' in vars(self):
            self._V3 = spspla.spsolve(self._Kc, self._V1.T)

        # Convert to dense if necessary
        if hasattr(self._V3, 'todense'):
            self._V3 = self._V3.todense()

        # If so, compute the full covariance matrix
        Sigma_post = var_prior.copy()
        Sigma_post -= self._V3 @ self._V3.T
        params[gppn.POSTERIORCOVARIANCE] = Sigma_post

        # Inform GPSolverModule that displacement-related info is returned
        params[gppn.FIELD] = 'u'

    def _get_prior_samples(self, params, globdat):

        # Get the prior samples on f
        super()._get_prior_samples(params, globdat)

        # Inform GPSolverModule that force-related info is returned
        params[gppn.FIELD] = 'f'

    def _get_posterior_samples(self, params, globdat):

        # Get the prior samples on f
        super()._get_posterior_samples(params, globdat)

        # Inform GPSolverModule that force-related info is returned
        params[gppn.FIELD] = 'f'

    def _get_param_opt(self, Sigma_fc):

        # Determine the optimal value of alpha
        L = np.linalg.cholesky(self._H @ Sigma_fc @ self._H.T)
        v = spla.solve_triangular(L, self._y)
        alpha2 = v.T @ v / self._nobs

        return np.sqrt(alpha2)

    def _get_sqrtSigma(self):

        if not '_sqrtSigma' in vars(self):
            if self._diagonalized:
                self._sqrtSigma = spsp.diags(np.sqrt(self._Sigma.sum(axis=1)), format='csr')
            else:
                if hasattr(self._Sigma, 'todense'):
                    self._sqrtSigma = spsp.csr_array(np.linalg.cholesky(self._Sigma.todense()))
                else:
                    self._sqrtSigma = np.linalg.cholesky(self._Sigma)


def declare(factory):
    factory.declare_model('GPf', GPfModel)
