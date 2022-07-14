import numpy as np
from scipy.linalg import solve_triangular
from scipy.sparse.linalg import spsolve

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
RANDOMOBS = 'randomObs'
SHAPE = 'shape'
INTSCHEME = 'intScheme'
PDNOISE = 'pdNoise'


class GPfModel(GPModel):


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

                Sigma[self._cdofs,:] = Sigma[:,self._cdofs] = 0.0
                Sigma += self._pdnoise2 * np.identity(self._dc)

                self._hyperparams[key] = self._get_param_opt(Sigma)

        super()._configure_prior(params, globdat)


    def _get_prior_mean(self, params, globdat):

        # Check if the prior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')

        if field == 'u':

            # Return the prior of the displacement field
            params[gppn.PRIORMEAN] = spsolve(self._Kc, self._m)

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
            params[gppn.POSTERIORMEAN] = spsolve(self._Kc, params[gppn.POSTERIORMEAN])

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
                self._V2 = spsolve(self._Kc, self._sqrtSigma)

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
                self._V3 = spsolve(self._Kc, self._V1.T)

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
            params[gppn.PRIORSAMPLES] = spsolve(self._Kc, params[gppn.PRIORSAMPLES])

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
            params[gppn.POSTERIORSAMPLES] = spsolve(self._Kc, params[gppn.POSTERIORSAMPLES])

        elif field == 'f':

            pass

        else:

            raise ValueError(field)


    def _get_param_opt(self, Sigma_fc):

        # Determine the optimal value of alpha
        L = np.linalg.cholesky(self._H @ Sigma_fc @ self._H.T)
        v = solve_triangular(L, self._y)
        alpha2 = v.T @ v / self._nobs

        return np.sqrt(alpha2)


def declare(factory):
    factory.declare_model('GPf', GPfModel)
