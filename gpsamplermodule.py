import numpy as np

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

class GPSamplerModule(LinsolveModule):

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

        # Get the rhs prior mean
        mf = self.get_neumann_vector(globdat)
        m = self._solver.solve(mf)
        params[gppn.PRIORMEAN] = m

        # Configure the GP based on the fine FEM results
        model.take_action(gpact.CONFIGUREFEM, params, globdat)
        model.take_action(gpact.CONFIGUREPRIOR, params, globdat)

        # Define a dictionary for the output params
        params = {}
        params[gppn.NSAMPLE] = self._nsample
        params[gppn.RNG] = np.random.default_rng(self._seed)

        # Get the prior samples
        model.take_action(gpact.GETPRIORSAMPLES, params, globdat)

        # Store the prior samples in globdat
        globdat['samples_u_prior'] = params[gppn.PRIORSAMPLES]
        globdat['samples_f_prior'] = self._solver.get_matrix() @ params[gppn.PRIORSAMPLES]

        # Update the prior samples to posterior samples
        model.take_action(gpact.KALMANUPDATE, params, globdat)

        # Store the updated samples in globdat
        globdat['samples_u_post'] = params[gppn.POSTERIORSAMPLES]
        globdat['samples_f_post'] = self._solver.get_matrix() @ params[gppn.POSTERIORSAMPLES]

        # Compute the prior and posterior mean
        globdat['f_prior'] = np.mean(globdat['samples_f_prior'], axis=1)
        globdat['u_prior'] = np.mean(globdat['samples_u_prior'], axis=1)
        globdat['f_post'] = np.mean(globdat['samples_f_post'], axis=1)
        globdat['u_post'] = np.mean(globdat['samples_u_post'], axis=1)

        # Compute the prior and posterior variance
        globdat['std_f_prior'] = np.std(globdat['samples_f_prior'], axis=1)
        globdat['std_u_prior'] = np.std(globdat['samples_u_prior'], axis=1)
        globdat['std_f_post'] = np.std(globdat['samples_f_post'], axis=1)
        globdat['std_u_post'] = np.std(globdat['samples_u_post'], axis=1)

        return output

    def shutdown(self, globdat):
        pass


def declare(factory):
    factory.declare_module('GPSampler', GPSamplerModule)
