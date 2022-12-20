import numpy as np
import scipy.sparse as spsp

from jive.solver.numba.cholesky import sparse_cholesky
from jive.solver.numba.spsolve import solve_triangular
from jive.fem.names import GPActions as gpact
from jive.fem.names import GPParamNames as gppn
from jive.fem.names import GlobNames as gn
from gpmodel import GPModel

PRIOR = 'prior'
EXPLICITINVERSE = 'explicitInverse'
DIAGONALIZED = 'diagonalized'
SOLVER = 'solver'
PRECONDITIONER = 'preconditioner'
TYPE = 'type'

class GPfModel(GPModel):

    def configure(self, props, globdat):

        super().configure(props, globdat)

        self._explicit_inverse = bool(eval(props.get(EXPLICITINVERSE, 'True')))
        self._diagonalized = props[PRIOR].get(DIAGONALIZED, False)

        solverprops = props.get(SOLVER, {})
        solver = solverprops.get(TYPE, 'cholmod')
        self._solver = globdat[gn.SOLVERFACTORY].get_solver(solver)
        self._solver.configure(solverprops, globdat)

        preconprops = props.get(PRECONDITIONER, {})
        self._precon = None
        precon = preconprops.get(TYPE)
        if precon is not None:
            self._precon = globdat[gn.PRECONFACTORY].get_precon(precon)
            self._precon.configure(preconprops, globdat)

    def _configure_fem(self, params, globdat):

        super()._configure_fem(params, globdat)

        if self._explicit_inverse:
            # Check if K^-1 should just be precomputed once explictly
            self._Kinv = np.linalg.inv(self._Kc.toarray())
        else:
            # Otherwise, configure the solver
            K = globdat.get(gn.MATRIX0)
            c = globdat.get(gn.CONSTRAINTS)

            # Update the solver
            self._solver.update(K, c, self._precon)

            # Make sure that the solver is in precon mode
            # (so it does not add the DBC contribution to the rhs)
            self._solver.precon_mode = True

    def _configure_prior(self, params, globdat):

        # Get the observation operator
        self._Phic = spsp.csr_array(self._Phic)
        self._H = self._Phic.T
        self._H = self._H.tocsr()

        # Define the mean in terms of the force vector as well
        self._m = self._Kc @ self._m

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

        # Solve to get the prior mean on u
        if self._explicit_inverse:
            params[gppn.PRIORMEAN] = self._Kinv @ params[gppn.PRIORMEAN]
        else:
            params[gppn.PRIORMEAN] = self._solver.solve(params[gppn.PRIORMEAN])

    def _get_posterior_mean(self, params, globdat):

        # Get the posterior mean on f
        super()._get_posterior_mean(params, globdat)

        # Solve to get the posterior mean on u
        if self._explicit_inverse:
            params[gppn.POSTERIORMEAN] = self._Kinv @ params[gppn.POSTERIORMEAN]
        else:
            params[gppn.POSTERIORMEAN] = self._solver.solve(params[gppn.POSTERIORMEAN])

    def _get_prior_covariance(self, params, globdat):

        # Get the prior covariance on f
        super()._get_prior_covariance(params, globdat)

        # Solve to get the prior covariance on u
        if self._explicit_inverse:
            params[gppn.PRIORCOVARIANCE] = self._Kinv @ params[gppn.PRIORCOVARIANCE] @ self._Kinv
        else:
            params[gppn.PRIORCOVARIANCE] = self._solver.solve(self._solver.solve(params[gppn.PRIORCOVARIANCE]).T)

    def _get_posterior_covariance(self, params, globdat):

        # Unsolve to get the prior covariance on f
        if not gppn.PRIORCOVARIANCE in params:
            self.take_action(gpact.GETPRIORCOVARIANCE, params, globdat)
        tmp = params[gppn.PRIORCOVARIANCE]
        params[gppn.PRIORCOVARIANCE] = self._Kc @ params[gppn.PRIORCOVARIANCE] @ self._Kc

        # Get the posterior covariance on f
        super()._get_posterior_covariance(params, globdat)

        # Solve to get the posterior covariance on u
        if self._explicit_inverse:
            params[gppn.POSTERIORCOVARIANCE] = self._Kinv @ params[gppn.POSTERIORCOVARIANCE] @ self._Kinv
        else:
            params[gppn.POSTERIORCOVARIANCE] = self._solver.solve(self._solver.solve(params[gppn.POSTERIORCOVARIANCE]).T)

        # Store the original prior covariance back
        params[gppn.PRIORCOVARIANCE] = tmp

    def _get_prior_samples(self, params, globdat):

        # Get the prior samples on f
        super()._get_prior_samples(params, globdat)

        # Solve to get the prior samples on u
        if self._explicit_inverse:
            params[gppn.PRIORSAMPLES] = self._Kinv @ params[gppn.PRIORSAMPLES]
        else:
            params[gppn.PRIORSAMPLES] = self._solver.solve(params[gppn.PRIORSAMPLES])

    def _get_posterior_samples(self, params, globdat):

        # Get the posterior samples on f
        super()._get_posterior_samples(params, globdat)

        # Solve to get the posterior samples on u
        if self._explicit_inverse:
            params[gppn.POSTERIORSAMPLES] = self._Kinv @ params[gppn.POSTERIORSAMPLES]
        else:
            params[gppn.POSTERIORSAMPLES] = self._solver.solve(params[gppn.POSTERIORSAMPLES])

    def _kalman_update(self, params, globdat):

        # Unsolve to get the prior samples on f
        tmp = params[gppn.PRIORSAMPLES]
        params[gppn.PRIORSAMPLES] = self._Kc @ params[gppn.PRIORSAMPLES]

        # Perform the Kalman update on f
        super()._kalman_update(params, globdat)

        # Solve to get the posterior samples on u
        if self._explicit_inverse:
            params[gppn.POSTERIORSAMPLES] = self._Kinv @ params[gppn.POSTERIORSAMPLES]
        else:
            params[gppn.POSTERIORSAMPLES] = self._solver.solve(params[gppn.POSTERIORSAMPLES])

        # Store the original prior samples back
        params[gppn.PRIORSAMPLES] = tmp

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
