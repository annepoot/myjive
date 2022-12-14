import numpy as np
import scipy.sparse as spsp

from jive.fem.names import GlobNames as gn
from jive.fem.names import ParamNames as pn
from jive.fem.names import Actions as act
from jive.fem.names import GPActions as gpact
from jive.fem.names import GPParamNames as gppn

from jive.implicit.linsolvemodule import LinsolveModule

GETUNITMASSMATRIX = 'getUnitMassMatrix'
EXPLICITINVERSE = 'explicitInverse'
NSAMPLE = 'nsample'
SEED = 'seed'

class GPExactModule(LinsolveModule):

    def init(self, props, globdat):

        # Initialize solvermodule first
        super().init(props, globdat)

        myprops = props[self._name]
        self._get_unit_mass_matrix = bool(eval(myprops.get(GETUNITMASSMATRIX, 'True')))
        self._explicit_inverse = bool(eval(myprops.get(EXPLICITINVERSE, 'True')))

        self._nsample = int(myprops.get(NSAMPLE,1))
        self._seed = eval(myprops.get(SEED,'None'))
        if not self._seed is None:
            self._seed = int(self._seed)

    def run(self, globdat):

        # Run solvermodule first
        output = super().run(globdat)

        model = globdat[gn.MODEL]

        # Optionally get the mass matrix
        if self._get_unit_mass_matrix:
            M = self._get_empty_matrix(globdat)
            params = {
                pn.MATRIX2: M,
                pn.UNITMATRIX: True
            }
            model.take_action(act.GETMATRIX2, params, globdat)

            # Store mass matrix in Globdat
            globdat[gn.MATRIX2] = M

        # Configure the GP based on the fine FEM results
        model.take_action(gpact.CONFIGUREFEM, params, globdat)
        model.take_action(gpact.CONFIGUREPRIOR, params, globdat)

        # Define a dictionary for the output params
        params = {}

        self._solver.precon_mode = True

        # Get K_inv if necessary
        if self._explicit_inverse:
            I = spsp.identity(self._solver.get_matrix().shape[0])
            K_inv = self._solver.solve(I)

        # Get the prior mean
        model.take_action(gpact.GETPRIORMEAN, params, globdat)

        # Store the prior mean in globdat
        if params[gppn.FIELD] == 'f':
            globdat['f_prior'] = params[gppn.PRIORMEAN]
            if self._explicit_inverse:
                globdat['u_prior'] = K_inv @ params[gppn.PRIORMEAN]
            else:
                globdat['u_prior'] = self._solver.solve(params[gppn.PRIORMEAN])
        else:
            globdat['u_prior'] = params[gppn.PRIORMEAN]
            globdat['f_prior'] = self._solver.get_matrix() @ params[gppn.PRIORMEAN]

        # Get the posterior mean
        model.take_action(gpact.GETPOSTERIORMEAN, params, globdat)

        # Store the posterior mean in globdat
        if params[gppn.FIELD] == 'f':
            globdat['f_post'] = params[gppn.POSTERIORMEAN]
            if self._explicit_inverse:
                globdat['u_post'] = K_inv @ params[gppn.POSTERIORMEAN]
            else:
                globdat['u_post'] = self._solver.solve(params[gppn.POSTERIORMEAN])
        else:
            globdat['u_post'] = params[gppn.POSTERIORMEAN]
            globdat['f_post'] = self._solver.get_matrix() @ params[gppn.POSTERIORMEAN]

        # Get the prior covariance
        model.take_action(gpact.GETPRIORCOVARIANCE, params, globdat)

        # Store the prior covariance in globdat
        if params[gppn.FIELD] == 'f':
            globdat['var_f_prior'] = params[gppn.PRIORCOVARIANCE]
            if self._explicit_inverse:
                globdat['var_u_prior'] = K_inv @ params[gppn.PRIORCOVARIANCE] @ K_inv
            else:
                globdat['var_u_prior'] = self._solver.solve(self._solver.solve(params[gppn.PRIORCOVARIANCE]).T)
        else:
            globdat['var_u_prior'] = params[gppn.PRIORCOVARIANCE]
            globdat['var_f_prior'] = self._solver.get_matrix() @ params[gppn.PRIORCOVARIANCE] @ self._solver.get_matrix()

        # Get the posterior covariance
        model.take_action(gpact.GETPOSTERIORCOVARIANCE, params, globdat)

        # Store the posterior covariance in globdat
        if params[gppn.FIELD] == 'f':
            globdat['var_f_post'] = params[gppn.POSTERIORCOVARIANCE]
            if self._explicit_inverse:
                globdat['var_u_post'] = K_inv @ params[gppn.POSTERIORCOVARIANCE] @ K_inv
            else:
                globdat['var_u_post'] = self._solver.solve(self._solver.solve(params[gppn.POSTERIORCOVARIANCE]).T)
        else:
            globdat['var_u_post'] = params[gppn.POSTERIORCOVARIANCE]
            globdat['var_f_post'] = self._solver.get_matrix() @ params[gppn.POSTERIORCOVARIANCE] @ self._solver.get_matrix()

        # Get the log likelihood and store it in globdat
        model.take_action(gpact.GETLOGLIKELIHOOD, params, globdat)
        globdat['logLikelihood'] = params[gppn.LOGLIKELIHOOD]

        # Set the params for the samples
        params[gppn.NSAMPLE] = self._nsample
        params[gppn.RNG] = np.random.default_rng(self._seed)

        # Get the prior samples
        model.take_action(gpact.GETPRIORSAMPLES, params, globdat)

        # Store the prior samples in globdat
        if params[gppn.FIELD] == 'f':
            globdat['samples_f_prior'] = params[gppn.PRIORSAMPLES]
            if self._explicit_inverse:
                globdat['samples_u_prior'] = K_inv @ params[gppn.PRIORSAMPLES]
            else:
                globdat['samples_u_prior'] = self._solver.solve(params[gppn.PRIORSAMPLES])
        else:
            globdat['samples_u_prior'] = params[gppn.PRIORSAMPLES]
            globdat['samples_f_prior'] = self._solver.get_matrix() @ params[gppn.PRIORSAMPLES]

        # Update the prior samples to posterior samples
        model.take_action(gpact.KALMANUPDATE, params, globdat)

        # Store the updated samples in globdat
        if params[gppn.FIELD] == 'f':
            globdat['samples_f_post'] = params[gppn.POSTERIORSAMPLES]
            if self._explicit_inverse:
                globdat['samples_u_post'] = K_inv @ params[gppn.POSTERIORSAMPLES]
            else:
                globdat['samples_u_post'] = self._solver.solve(params[gppn.POSTERIORSAMPLES])
        else:
            globdat['samples_u_post'] = params[gppn.POSTERIORSAMPLES]
            globdat['samples_f_post'] = self._solver.get_matrix() @ params[gppn.POSTERIORSAMPLES]

        self._solver.precon_mode = False

        return output

    def shutdown(self, globdat):
        pass


def declare(factory):
    factory.declare_module('GPSolver', GPExactModule)
