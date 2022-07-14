import numpy as np
import scipy.sparse as spsp
import scipy.linalg as spla
import scipy.sparse.linalg as spspla

from names import GPActions as gpact
from names import GPParamNames as gppn
from names import GlobNames as gn
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
        self._m = self._Kc @ self._m

        # Check if alpha or beta should be optimized
        if len(self._hyperparams) == 1:
            key, value = list(self._hyperparams.items())[0]

            if value == 'opt':

                if key == 'alpha':
                    Sigma = globdat[gn.MATRIX2]
                elif key == 'beta':
                    Sigma == globdat[gn.MATRIX0]
                else:
                    raise ValueError('cannot find optimal value for ' + key)

                # Apply boundary conditions to the prior
                Sigma = self._apply_covariance_bcs(Sigma)

                self._hyperparams[key] = self._get_param_opt(Sigma)

        super()._configure_prior(params, globdat)


    def _get_prior_mean(self, params, globdat):

        # Check if the prior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')

        if field == 'u':

            # Return the prior of the displacement field
            params[gppn.PRIORMEAN] = spspla.spsolve(self._Kc, self._m)

        elif field == 'f':

            # Return the prior of the force field
            params[gppn.PRIORMEAN] = self._m

        else:

            raise ValueError(field)


    def _get_posterior_mean(self, params, globdat):

        # Get the posterior mean on f first
        super()._get_posterior_mean(params, globdat)

        # Check if the prior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')

        if field == 'u':

            # Solve to get the posterior of the displacement field
            params[gppn.POSTERIORMEAN] = spspla.spsolve(self._Kc, params[gppn.POSTERIORMEAN])

        elif field == 'f':

            pass

        else:

            raise ValueError(field)


    def _get_prior_covariance(self, params, globdat):

        # Get the prior covariance on f first
        super()._get_prior_covariance(params, globdat)

        # Check if the prior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')
        fullSigma = params.get(gppn.FULLCOVARIANCE, False)

        if field == 'u':

            # Get the relevant matrices
            self._get_sqrtSigma()

            # Solve the system for each dof
            if not '_V2' in vars(self):
                self._V2 = spspla.spsolve(self._Kc, self._sqrtSigma)

            # Convert to dense if necessary
            if hasattr(self._V2, 'todense'):
                self._V2 = self._V2.todense()

            # Check if the full covariance matrix should be returned
            if fullSigma:

                # If so, compute the full covariance matrix
                Sigma_prior = self._V2 @ self._V2.T
                params[gppn.PRIORCOVARIANCE] = Sigma_prior

            else:

                # Compute only the diagonal of the covariance matrix
                var_prior = np.zeros(self._dc)
                for i in range(self._dc):
                    var_prior[i] = self._V2[i,:] @ self._V2[i,:].T

                params[gppn.PRIORCOVARIANCE] = var_prior

        elif field == 'f':

            pass

        else:

            raise ValueError(field)


    def _get_posterior_covariance(self, params, globdat):

        # Get the posterior covariance on f first
        super()._get_posterior_covariance(params, globdat)

        # Check if the posterior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')
        fullSigma = params.get(gppn.FULLCOVARIANCE, False)

        # Check if the prior variance is given in the params, otherwise, take an action to obtain it
        if not gppn.PRIORCOVARIANCE in params:
            self.take_action(gpact.GETPRIORCOVARIANCE, params, globdat)

        # Get the prior variance from the params
        var_prior = params[gppn.PRIORCOVARIANCE]

        if field == 'u':

            # Get the relevant matrices
            self._get_V1()

            # Solve the system for each coarse dof
            if not '_V3' in vars(self):
                self._V3 = spspla.spsolve(self._Kc, self._V1.T)

            # Convert to dense if necessary
            if hasattr(self._V3, 'todense'):
                self._V3 = self._V3.todense()

            # Check if the full covariance matrix should be returned
            if fullSigma:

                # If so, compute the full covariance matrix
                Sigma_post = var_prior.copy()
                Sigma_post -= self._V3 @ self._V3.T
                params[gppn.POSTERIORCOVARIANCE] = Sigma_post

            else:

                # If not, compute only the diagonal of the covariance matrix
                var_post = var_prior.copy()
                for i in range(self._dc):
                    var_post[i] -= self._V3[i,:] @ self._V3[i,:].T

                params[gppn.POSTERIORCOVARIANCE] = var_post

        elif field == 'f':

            pass

        else:

            raise ValueError(field)


    def _get_prior_samples(self, params, globdat):

        # Get the prior samples on f first
        super()._get_prior_samples(params, globdat)

        # Check if the prior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')

        if field == 'u':

            # Compute the corresponding displacement field for each sample
            params[gppn.PRIORSAMPLES] = spspla.spsolve(self._Kc, params[gppn.PRIORSAMPLES])

        elif field == 'f':

            pass

        else:

            raise ValueError(field)


    def _get_posterior_samples(self, params, globdat):

        # Get the prior samples on f first
        super()._get_posterior_samples(params, globdat)

        # Check if the prior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')

        if field == 'u':

            # Compute the corresponding displacement field for each sample
            params[gppn.POSTERIORSAMPLES] = spspla.spsolve(self._Kc, params[gppn.POSTERIORSAMPLES])

        elif field == 'f':

            pass

        else:

            raise ValueError(field)


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
