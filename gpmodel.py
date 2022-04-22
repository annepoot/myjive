import numpy as np

from names import GPActions as gpact
from names import GPParamNames as gppn
from names import GlobNames as gn
from model import Model, ModelFactory
from constrainer import Constrainer

NOBS = 'nobs'
OBSNOISE = 'obsNoise'
ALPHA = 'alpha'
RANDOMOBS = 'randomObs'

class GPModel(Model):

    def take_action(self, action, params, globdat):
        print('BarModel taking action', action)

        if action == gpact.GETPRIORMEAN:
            self._get_prior_mean(params, globdat)
        elif action == gpact.GETPRIORCOVARIANCE:
            self._get_prior_covariance(params, globdat)
        elif action == gpact.GETPOSTERIORMEAN:
            self._get_posterior_mean(params, globdat)
        elif action == gpact.GETPOSTERIORCOVARIANCE:
            self._get_posterior_covariance(params, globdat)


    def configure(self, props, globdat):

        self._dc = globdat[gn.DOFSPACE].dof_count()

        # Get the number of observations either as a ratio or as a fixed integer
        self._nobs = float(props.get(NOBS,self._dc))
        if self._nobs < 1:
            self._nobs = self._dc * self._nobs
        self._nobs = int(np.round(self._nobs))

        # Get the observational noise
        self._noise2 = float(props.get(OBSNOISE))**2

        # Get the alpha parameter, or use the optimal alpha
        self._alpha2 = props.get(ALPHA, 'opt')
        if self._alpha2 != 'opt':
            self._alpha2 = float(self._alpha2)**2


    def configure_fem(self, globdat):

        # Get K, M and f from globdat
        K = globdat.get(gn.MATRIX0)
        M = globdat.get(gn.MATRIX2)
        f = globdat.get(gn.EXTFORCE)
        c = globdat.get(gn.CONSTRAINTS)

        # Constrain K, M and f
        self._Mc = c.constrain(M, f)[0]
        self._Kc, self._fc = c.constrain(K, f)

        # Get phi from Globdat
        self._phi = self._get_phis()[0]

        # Compute the observation vector
        self._y = self._phi.T @ self._fc

        # Set alpha to the optimal value if necessary
        if self._alpha2 == 'opt':
            self._alpha2 = self._get_alpha2_opt()
            print('Setting alpha^2 to optimal value: {:.4f}'.format(self._alpha2))


    def _get_prior_mean(self, params, globdat):

        #########################
        # PRIOR MEAN ON u and f #
        #########################

        # Only a 0-mean prior has been implemented
        params[gppn.PRIORMEAN] = np.zeros(self._dc)


    def _get_posterior_mean(self, params, globdat):

        # Get the observation covariance matrix
        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._alpha2 * self._phi.T @ self._Mc @ self._phi + np.identity(self._nobs) * self._noise2

        # Do the cholesky decomposition only if necessary
        if not '_sqrtObs' in vars(self):
            self._sqrtObs = np.linalg.cholesky(self._Sigma_obs)

        # Solve the system for the observed forces
        if not '_v0' in vars(self):
            self._v0 = np.linalg.solve(self._sqrtObs, self._y)

        # Check if the prior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')

        if field == 'u':

            #######################
            # POSTERIOR MEAN ON u #
            #######################

            if not '_u_post' in vars(self):

                # Get the posterior of the force field
                if not '_f_post' in vars(self):
                    self._f_post = self._alpha2 * self._Mc @ self._phi @ np.linalg.solve(self._sqrtObs.T, self._v0)

                # Solve to get the posterior of the displacement field
                self._u_post = np.linalg.solve(self._Kc, self._f_post)

            # Return the posterior of the displacement field
            params[gppn.POSTERIORMEAN] = self._u_post

        elif field == 'f':

            #######################
            # POSTERIOR MEAN ON f #
            #######################

            # Get the posterior of the force field
            if not '_f_post' in vars(self):
                self._f_post = self._alpha2 * self._Mc @ self._phi @ np.linalg.solve(self._sqrtObs.T, self._v0)

            # Return the posterior of the force field
            params[gppn.POSTERIORMEAN] = self._f_post

        else:

            raise ValueError(field)


    def _get_prior_covariance(self, params, globdat):

        # Check if the prior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')
        fullSigma = params.get(gppn.FULLCOVARIANCE, False)

        if field == 'u':

            #########################
            # PRIOR COVARIANCE ON u #
            #########################

            # Sigma = alpha2 * inv(K) * M * inv(K)

            # Do the cholesky decomposition only if necessary
            if not '_sqrtM' in vars(self):
                self._sqrtM = np.linalg.cholesky(self._Mc)

            # Solve the system for each dof
            if not '_V2' in vars(self):
                self._V2 = np.linalg.solve(self._Kc, self._sqrtM)

            # Check if the full covariance matrix should be returned
            if fullSigma:

                # If so, compute the full covariance matrix
                Sigma_prior = self._alpha2 * self._V2 @ self._V2.T

            else:

                # Compute only the diagonal of the covariance matrix
                sigma_prior = np.zeros(self._dc)
                for i in range(self._dc):
                    sigma_prior[i] = self._alpha2 * self._V2[i,:] @ self._V2[i,:].T

        elif field == 'f':

            #########################
            # PRIOR COVARIANCE ON f #
            #########################

            # Sigma = alpha2 * M

            # Check if the full covariance matrix should be returned
            if fullSigma:

                # If so, compute the full covariance matrix
                Sigma_prior = self._alpha2 * self._Mc

            else:

                # If not, compute only the diagonal of the covariance matrix
                sigma_prior = self._alpha2 * self._Mc.diagonal()

        else:

            raise ValueError(field)

        # Return either the full covariance matrix, or only its diagonal
        if fullSigma:
            params[gppn.PRIORCOVARIANCE] = Sigma_prior
        else:
            params[gppn.PRIORCOVARIANCE] = sigma_prior


    def _get_posterior_covariance(self, params, globdat):

        # Check if the posterior on u, eps or f should be obtained
        field = params.get(gppn.FIELD, 'u')
        fullSigma = params.get(gppn.FULLCOVARIANCE, False)

        # Check if the prior variance is given in the params, otherwise, take an action to obtain it
        if not gppn.PRIORCOVARIANCE in params:
            self.take_action(gpact.GETPRIORCOVARIANCE, params, globdat)

        # Get the prior variance from the params
        sigma_prior = params[gppn.PRIORCOVARIANCE]

        # Get the observation covariance matrix
        if not '_Sigma_obs' in vars(self):
            self._Sigma_obs = self._alpha2 * self._phi.T @ self._Mc @ self._phi + np.identity(self._nobs) * self._noise2

        # Do the cholesky decomposition only if necessary
        if not '_sqrtObs' in vars(self):
            self._sqrtObs = np.linalg.cholesky(self._Sigma_obs)

        # Solve the inverse observation covariance once for each observation
        if not '_V1' in vars(self):
            self._V1 = np.linalg.solve(self._sqrtObs, self._phi.T @ self._Mc)

        if field == 'u':

            #############################
            # POSTERIOR COVARIANCE ON u #
            #############################

            # Sigma = alpha2 * inv(K) * M * inv(K) - alpha4 * inv(K) * M * phi * inv(alpha2 * phi.T * M * phi.T + Sigma_e) * phi.T * M * inv(K)

            # Solve the system for each dof
            if not '_V3' in vars(self):
                self._V3 = np.linalg.solve(self._Kc, self._V1.T)

            # Check if the full covariance matrix should be returned
            if fullSigma:

                # If so, compute the full covariance matrix
                Sigma_post = sigma_prior.copy()
                Sigma_post -= self._alpha2**2 * self._V3 @ self._V3.T

            else:

                # If not, compute only the diagonal of the covariance matrix
                sigma_post = sigma_prior.copy()
                for i in range(self._dc):
                    sigma_post[i] -= self._alpha2**2 * self._V3[i,:] @ self._V3[i,:].T

        elif field == 'f':

            #############################
            # POSTERIOR COVARIANCE ON f #
            #############################

            # Sigma = alpha2 * M - alpha4 * M * phi * inv(alpha2 * phi.T * M * phi + Sigma_e) * phi.T * M

            # Check if the full covariance matrix should be returned
            if fullSigma:

                # If so, compute the full covariance matrix
                Sigma_post = sigma_prior.copy()
                Sigma_post -= self._alpha2**2 * self._V1.T @ self._V1

            else:

                # If not, compute only the diagonal of the covariance matrix
                sigma_post = sigma_prior.copy()
                for i in range(self._dc):
                    sigma_post[i] -= self._alpha2**2 * self._V1[:,i].T @ self._V1[:,i]

        else:

            raise ValueError(field)

        # Return either the full covariance matrix, or only its diagonal
        if fullSigma:
            params[gppn.POSTERIORCOVARIANCE] = Sigma_post
        else:
            params[gppn.POSTERIORCOVARIANCE] = sigma_post


    def _get_alpha2_opt(self):

        #################
        # OPTIMAL ALPHA #
        #################

        # If so, determine the optimal value of alpha
        L = np.linalg.cholesky(self._phi.T @ self._Mc @ self._phi)
        v = np.linalg.solve(L, self._y)
        alpha2 = v.T @ v / self._nobs

        return alpha2


    def _get_phis(self):

        # Store the number of dofs, and number of observatiosn
        dc = self._dc
        nobs = self._nobs

        # Create an empty array
        phi = np.zeros((dc, nobs))
        phi_sub = np.zeros((dc, nobs))

        # Get the indices of the nodes that will be observed
        vec_obs = np.rint(np.linspace(0, dc-1, nobs))
        vec_obs = vec_obs.astype(int)

        for i in range(nobs):

            phi[vec_obs[i], i] = 1
            phi_sub[vec_obs[i], i] = 1

        j = 0
        j_old = 0

        for i in range(nobs-1):
            while j < dc:
                if np.isclose(phi[j,i+1], 1):
                    phi[j_old:j+1,i] = np.linspace(1, 0, j-j_old+1)
                    phi[j_old:j+1,i+1] = np.linspace(0, 1, j-j_old+1)
                    j_old = j
                    break

                j = j+1

        return phi, phi_sub


def declare(factory):
    factory.declare_model('GP', GPModel)
