import numpy as np

from jive.fem.names import GlobNames as gn
from jive.fem.names import ParamNames as pn
from jive.fem.names import Actions as act
from jive.fem.names import GPActions as gpact
from jive.fem.names import GPParamNames as gppn

from jive.implicit.linsolvemodule import LinsolveModule
from jive.util.table import Table

GETUNITMASSMATRIX = 'getUnitMassMatrix'
POSTPROJECT = 'postproject'
PRIORMEAN = 'priorMean'
NSAMPLE = 'nsample'
SEED = 'seed'

class GPModule(LinsolveModule):

    def init(self, props, globdat):

        # Initialize solvermodule first
        super().init(props, globdat)

        myprops = props[self._name]
        self._get_unit_mass_matrix = bool(eval(myprops.get(GETUNITMASSMATRIX, 'True')))

        self._postproject = bool(eval(myprops.get(POSTPROJECT, 'False')))

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
        globdat[gn.GP][gn.COVARIANCE] = {}
        globdat[gn.GP][gn.COVARIANCE][gn.PRIOR] = {}
        globdat[gn.GP][gn.COVARIANCE][gn.POSTERIOR] = {}
        globdat[gn.GP][gn.STD] = {}
        globdat[gn.GP][gn.STD][gn.PRIOR] = {}
        globdat[gn.GP][gn.STD][gn.POSTERIOR] = {}
        globdat[gn.GP][gn.SAMPLES] = {}
        globdat[gn.GP][gn.SAMPLES][gn.PRIOR] = {}
        globdat[gn.GP][gn.SAMPLES][gn.POSTERIOR] = {}

        return output

    def shutdown(self, globdat):
        pass

    def take_sample_action(self, samples, globdat):

        model = globdat[gn.MODEL]

        # Define a dictionary for the output params
        params = {}
        params[gppn.NSAMPLE] = self._nsample
        params[gppn.RNG] = np.random.default_rng(self._seed)
        params[gppn.PRIORSAMPLES] = np.zeros((self._dc, self._nsample))
        params[gppn.POSTERIORSAMPLES] = np.zeros((self._dc, self._nsample))

        # Get the prior samples
        model.take_action(gpact.GETPRIORSAMPLES, params, globdat)

        # Store the prior samples in globdat
        samples[gn.PRIOR][gn.STATE0]   = params[gppn.PRIORSAMPLES]
        samples[gn.PRIOR][gn.EXTFORCE] = self._solver.get_matrix() @ params[gppn.PRIORSAMPLES]

        # Update the prior samples to posterior samples
        model.take_action(gpact.KALMANUPDATE, params, globdat)

        # Store the updated samples in globdat
        samples[gn.POSTERIOR][gn.STATE0]   = params[gppn.POSTERIORSAMPLES]
        samples[gn.POSTERIOR][gn.EXTFORCE] = self._solver.get_matrix() @ params[gppn.POSTERIORSAMPLES]

        # Project the samples over the coarse space if needed
        if self._postproject:
            model.take_action(gpact.PROJECTSAMPLES, params, globdat)

            # Store the projected samples in globdat
            samples[gn.PRIOR][gn.STATE0]     = params[gppn.PRIORSAMPLES]
            samples[gn.POSTERIOR][gn.STATE0] = params[gppn.POSTERIORSAMPLES]

    def take_sample_field_action(self, samples, globdat):

        model = globdat[gn.MODEL]

        # Loop over the different fields
        for name in self._tnames:

            # Get the fine solution table for that field
            globtable = globdat[gn.TABLES][name]

            # Get the tables for both prior and posterior
            for distribution in [gn.PRIOR, gn.POSTERIOR]:

                # Create a new sub-dictionary
                samples[distribution][name] = {}

                # Add a sample matrix for each component
                for comp in globtable.get_column_names():
                    samples[distribution][name][comp] = np.zeros((globdat[gn.NSET].size(), self._nsample))

                # Go over each state0 sample
                for i, sample in enumerate(samples[distribution][gn.STATE0].T):

                    # Define a dictionary for the output params
                    params = {}
                    params[pn.TABLE] = Table()
                    params[pn.TABLENAME] = name
                    params[pn.TABLEWEIGHTS] = np.zeros(len(globdat[gn.NSET]))
                    params[pn.SOLUTION] = sample
                    params[pn.VERBOSE] = False

                    # Get the relevant fields
                    model.take_action(act.GETTABLE, params, globdat)

                    # Add the field to the sample matrices for each component
                    for comp in params[pn.TABLE].get_column_names():
                        samples[distribution][name][comp][:,i] = params[pn.TABLE][comp]

    def take_mean_action(self, mean, globdat):

        model = globdat[gn.MODEL]

        # Define a dictionary for the output params
        params = {}

        # Get the prior and posterior mean
        model.take_action(gpact.GETPRIORMEAN, params, globdat)
        model.take_action(gpact.GETPOSTERIORMEAN, params, globdat)

        # Store the prior and posterior mean in globdat
        mean[gn.PRIOR][gn.STATE0]       = params[gppn.PRIORMEAN]
        mean[gn.PRIOR][gn.EXTFORCE]     = self._solver.get_matrix() @ params[gppn.PRIORMEAN]
        mean[gn.POSTERIOR][gn.STATE0]   = params[gppn.POSTERIORMEAN]
        mean[gn.POSTERIOR][gn.EXTFORCE] = self._solver.get_matrix() @ params[gppn.POSTERIORMEAN]

    def take_variance_action(self, cov, globdat):

        model = globdat[gn.MODEL]

        # Define a dictionary for the output params
        params = {}

        # Get the prior and posterior covariance
        model.take_action(gpact.GETPRIORCOVARIANCE, params, globdat)
        model.take_action(gpact.GETPOSTERIORCOVARIANCE, params, globdat)

        # Store the prior and posterior covariance in globdat
        cov[gn.PRIOR][gn.STATE0]       = params[gppn.PRIORCOVARIANCE]
        cov[gn.PRIOR][gn.EXTFORCE]     = self._solver.get_matrix() @ params[gppn.PRIORCOVARIANCE] @ self._solver.get_matrix()
        cov[gn.POSTERIOR][gn.STATE0]   = params[gppn.POSTERIORCOVARIANCE]
        cov[gn.POSTERIOR][gn.EXTFORCE] = self._solver.get_matrix() @ params[gppn.POSTERIORCOVARIANCE] @ self._solver.get_matrix()


def declare(factory):
    factory.declare_module('GP', GPModule)
