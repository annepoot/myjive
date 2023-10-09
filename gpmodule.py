import numpy as np

from jive.fem.names import GlobNames as gn
from jive.fem.names import ParamNames as pn
from jive.fem.names import Actions as act
from jive.fem.names import GPActions as gpact
from jive.fem.names import GPParamNames as gppn

from jive.implicit.linsolvemodule import LinsolveModule
from jive.util.table import Table
from jive.util.xtable import to_xtable

from jive.implicit.linsolvemodule import GETSTRAINMATRIX

GETUNITMASSMATRIX = 'getUnitMassMatrix'
EXPLICITINVERSE = 'explicitInverse'
POSTPROJECT = 'postproject'
NSAMPLE = 'nsample'
SEED = 'seed'

class GPModule(LinsolveModule):

    def init(self, props, globdat):

        # Initialize solvermodule first
        super().init(props, globdat)

        myprops = props[self._name]
        self._get_unit_mass_matrix = bool(eval(myprops.get(GETUNITMASSMATRIX, 'True')))
        self._get_strain_matrix = bool(eval(myprops.get(GETSTRAINMATRIX, 'True')))
        self._explicit_inverse = bool(eval(props.get(EXPLICITINVERSE, 'True')))
        self._postproject = bool(eval(myprops.get(POSTPROJECT, 'False')))

        self._nsample = int(myprops.get(NSAMPLE,1))
        self._seed = eval(myprops.get(SEED,'None'))
        if not self._seed is None:
            self._seed = int(self._seed)

    def run(self, globdat):

        # Run solvermodule first
        output = super().run(globdat)

        # Check if we should invert K explicitly
        if self._explicit_inverse:
            self._Kinv = np.linalg.inv(self._solver.get_matrix().toarray())

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

        # Get the prior and posterior samples
        model.take_action(gpact.GETPRIORSAMPLES, params, globdat)
        model.take_action(gpact.KALMANUPDATE, params, globdat)

        # Store the prior and posterior samples in globdat
        field = params[gppn.FIELD]

        samples[gn.PRIOR][field]     = params[gppn.PRIORSAMPLES]
        samples[gn.POSTERIOR][field] = params[gppn.POSTERIORSAMPLES]

        # Convert state0 to extForce or vice versa
        self._state0_extforce_conversion(samples[gn.PRIOR], field, 'samples')
        self._state0_extforce_conversion(samples[gn.POSTERIOR], field, 'samples')

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
                    nodecount = len(globdat[gn.NSET])

                    # Define a dictionary for the output params
                    params = {}
                    params[pn.TABLE] = Table(size=nodecount)
                    params[pn.TABLENAME] = name
                    params[pn.TABLEWEIGHTS] = np.zeros(nodecount)
                    params[pn.SOLUTION] = sample
                    params[pn.VERBOSE] = False

                    # Get the relevant fields
                    model.take_action(act.GETTABLE, params, globdat)

                    to_xtable(params[pn.TABLE])

                    for jcol in range(params[pn.TABLE].column_count()):
                        values = params[pn.TABLE].get_col_values(None, jcol)
                        params[pn.TABLE].set_col_values(None, jcol, values / params[pn.TABLEWEIGHTS])

                    params[pn.TABLE].to_table()

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
        field = params[gppn.FIELD]

        mean[gn.PRIOR][field]     = params[gppn.PRIORMEAN]
        mean[gn.POSTERIOR][field] = params[gppn.POSTERIORMEAN]

        # Convert state0 to extForce or vice versa
        self._state0_extforce_conversion(mean[gn.PRIOR], field, 'vector')
        self._state0_extforce_conversion(mean[gn.POSTERIOR], field, 'vector')

    def take_variance_action(self, cov, globdat):

        model = globdat[gn.MODEL]

        # Define a dictionary for the output params
        params = {}

        # Get the prior and posterior covariance
        model.take_action(gpact.GETPRIORCOVARIANCE, params, globdat)
        model.take_action(gpact.GETPOSTERIORCOVARIANCE, params, globdat)

        # Store the prior and posterior covariance in globdat
        field = params[gppn.FIELD]

        cov[gn.PRIOR][field]     = params[gppn.PRIORCOVARIANCE]
        cov[gn.POSTERIOR][field] = params[gppn.POSTERIORCOVARIANCE]

        # Convert state0 to extForce or vice versa
        self._state0_extforce_conversion(cov[gn.PRIOR], field, 'matrix')
        self._state0_extforce_conversion(cov[gn.POSTERIOR], field, 'matrix')

    def _state0_extforce_conversion(self, target, field, shape):
        assert shape in ['vector', 'matrix', 'samples']

        if field == gn.STATE0:
            if shape == 'vector' or shape == 'samples':
                target[gn.EXTFORCE] = self._solver.get_matrix() @ target[gn.STATE0]
            elif shape == 'matrix':
                target[gn.EXTFORCE] = self._solver.get_matrix() @ target[gn.STATE0] @ self._solver.get_matrix()
            else:
                raise ValueError('shape should be "vector", "matrix" or "samples"')
        elif field == gn.EXTFORCE:
            if self._explicit_inverse:
                if shape == 'vector' or shape == 'samples':
                    target[gn.STATE0] = self._Kinv @ target[gn.EXTFORCE]
                elif shape == 'matrix':
                    target[gn.STATE0] = self._Kinv @ target[gn.EXTFORCE] @ self._Kinv
                else:
                    raise ValueError('shape should be "vector", "matrix" or "samples"')
            else:
                self._solver.precon_mode = True
                if shape == 'vector':
                    target[gn.STATE0] = self._solver.solve(target[gn.EXTFORCE])
                elif shape == 'samples':
                    target[gn.STATE0] = self._solver.multisolve(target[gn.EXTFORCE])
                elif shape == 'matrix':
                    target[gn.STATE0] = self._solver.multisolve(self._solver.multisolve(target[gn.EXTFORCE]).T)
                else:
                    raise ValueError('shape should be "vector", "matrix" or "samples"')
                self._solver.precon_mode = False
        else:
            raise ValueError('field should be either STATE0 or EXTFORCE')


def declare(factory):
    factory.declare_module('GP', GPModule)
