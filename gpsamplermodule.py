import numpy as np

from jive.fem.names import GlobNames as gn
from jive.fem.names import ParamNames as pn
from jive.fem.names import Actions as act
from jive.fem.names import GPActions as gpact
from jive.fem.names import GPParamNames as gppn

from jive.implicit.linsolvemodule import LinsolveModule

GETUNITMASSMATRIX = 'getUnitMassMatrix'
EXPLICITINVERSE = 'explicitInverse'
POSTPROJECT = 'postproject'
PRIORMEAN = 'priorMean'
NSAMPLE = 'nsample'
SEED = 'seed'

class GPSamplerModule(LinsolveModule):

    def init(self, props, globdat):

        # Initialize solvermodule first
        super().init(props, globdat)

        myprops = props[self._name]
        self._get_unit_mass_matrix = bool(eval(myprops.get(GETUNITMASSMATRIX, 'True')))
        self._explicit_inverse = bool(eval(myprops.get(EXPLICITINVERSE, 'True')))

        self._postproject = bool(eval(props.get(POSTPROJECT, 'False')))

        self._priormean = myprops.get(PRIORMEAN, 'zero')
        if self._priormean not in ['zero', 'dirichlet', 'neumann']:
            raise ValueError('priorMean has to be "zero", "dirichlet" or "neumann".')

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
        if self._priormean == 'zero':
            m = np.zeros_like(globdat[gn.STATE0])
        else:
            if self._priormean == 'dirichlet':
                mf = np.zeros_like(globdat[gn.EXTFORCE])
            elif self._priormean == 'neumann':
                mf = self.get_neumann_vector(globdat)
            else:
                raise ValueError('priorMean has to be "zero", "dirichlet" or "neumann".')

            m = self._solver.solve(mf)

        params[gppn.PRIORMEAN] = m

        # Configure the GP based on the fine FEM results
        model.take_action(gpact.CONFIGUREFEM, params, globdat)
        model.take_action(gpact.CONFIGUREPRIOR, params, globdat)

        # Create a dictionary for the gp output
        globdat[gn.GP] = {}
        globdat[gn.GP][gn.MEAN] = {}
        globdat[gn.GP][gn.MEAN][gn.PRIOR] = {}
        globdat[gn.GP][gn.MEAN][gn.POSTERIOR] = {}
        globdat[gn.GP][gn.STD] = {}
        globdat[gn.GP][gn.STD][gn.PRIOR] = {}
        globdat[gn.GP][gn.STD][gn.POSTERIOR] = {}
        globdat[gn.GP][gn.SAMPLES] = {}
        globdat[gn.GP][gn.SAMPLES][gn.PRIOR] = {}
        globdat[gn.GP][gn.SAMPLES][gn.POSTERIOR] = {}

        # Define a dictionary for the output params
        params = {}
        params[gppn.NSAMPLE] = self._nsample
        params[gppn.RNG] = np.random.default_rng(self._seed)

        # Get the prior samples
        model.take_action(gpact.GETPRIORSAMPLES, params, globdat)

        # Store the prior samples in globdat
        globdat[gn.GP][gn.SAMPLES][gn.PRIOR][gn.STATE0]   = params[gppn.PRIORSAMPLES]
        globdat[gn.GP][gn.SAMPLES][gn.PRIOR][gn.EXTFORCE] = self._solver.get_matrix() @ params[gppn.PRIORSAMPLES]

        # Update the prior samples to posterior samples
        model.take_action(gpact.KALMANUPDATE, params, globdat)

        # Store the updated samples in globdat
        globdat[gn.GP][gn.SAMPLES][gn.POSTERIOR][gn.STATE0]   = params[gppn.POSTERIORSAMPLES]
        globdat[gn.GP][gn.SAMPLES][gn.POSTERIOR][gn.EXTFORCE] = self._solver.get_matrix() @ params[gppn.POSTERIORSAMPLES]

        # Project the samples over the coarse space if needed
        if self._postproject:
            model.take_action(gpact.PROJECTSAMPLES, params, globdat)

            # Store the projected samples in globdat
            globdat[gn.GP][gn.SAMPLES][gn.PRIOR][gn.STATE0]     = params[gppn.PRIORSAMPLES]
            globdat[gn.GP][gn.SAMPLES][gn.POSTERIOR][gn.STATE0] = params[gppn.POSTERIORSAMPLES]

        # Compute the prior and posterior mean
        globdat[gn.GP][gn.MEAN][gn.PRIOR][gn.STATE0]       = np.mean(globdat[gn.GP][gn.SAMPLES][gn.PRIOR][gn.STATE0], axis=1)
        globdat[gn.GP][gn.MEAN][gn.PRIOR][gn.EXTFORCE]     = np.mean(globdat[gn.GP][gn.SAMPLES][gn.PRIOR][gn.EXTFORCE], axis=1)
        globdat[gn.GP][gn.MEAN][gn.POSTERIOR][gn.STATE0]   = np.mean(globdat[gn.GP][gn.SAMPLES][gn.POSTERIOR][gn.STATE0], axis=1)
        globdat[gn.GP][gn.MEAN][gn.POSTERIOR][gn.EXTFORCE] = np.mean(globdat[gn.GP][gn.SAMPLES][gn.POSTERIOR][gn.EXTFORCE], axis=1)

        # Compute the prior and posterior variance
        globdat[gn.GP][gn.STD][gn.PRIOR][gn.STATE0]       = np.std(globdat[gn.GP][gn.SAMPLES][gn.PRIOR][gn.STATE0], axis=1)
        globdat[gn.GP][gn.STD][gn.PRIOR][gn.EXTFORCE]     = np.std(globdat[gn.GP][gn.SAMPLES][gn.PRIOR][gn.EXTFORCE], axis=1)
        globdat[gn.GP][gn.STD][gn.POSTERIOR][gn.STATE0]   = np.std(globdat[gn.GP][gn.SAMPLES][gn.POSTERIOR][gn.STATE0], axis=1)
        globdat[gn.GP][gn.STD][gn.POSTERIOR][gn.EXTFORCE] = np.std(globdat[gn.GP][gn.SAMPLES][gn.POSTERIOR][gn.EXTFORCE], axis=1)

        return output

    def shutdown(self, globdat):
        pass


def declare(factory):
    factory.declare_module('GPSampler', GPSamplerModule)
