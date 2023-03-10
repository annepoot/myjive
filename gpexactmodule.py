import numpy as np

from jive.fem.names import GlobNames as gn
from jive.fem.names import GPActions as gpact
from jive.fem.names import GPParamNames as gppn

from gpmodule import GPModule

GETUNITMASSMATRIX = 'getUnitMassMatrix'
POSTPROJECT = 'postproject'
PRIORMEAN = 'priorMean'
NSAMPLE = 'nsample'
SEED = 'seed'

class GPExactModule(GPModule):

    def run(self, globdat):

        # Run gpmodule (and solvermodule) first
        output = super().run(globdat)

        model = globdat[gn.MODEL]

        # Define a dictionary for the output params
        params = {}

        # Get the sub-dictionaries of the mean, covariance, std and samples
        mean = globdat[gn.GP][gn.MEAN]
        cov = globdat[gn.GP][gn.COVARIANCE]
        std = globdat[gn.GP][gn.STD]
        samples = globdat[gn.GP][gn.SAMPLES]

        # Get the mean
        self.take_mean_action(mean, globdat)

        # Get the covariance
        self.take_variance_action(cov, globdat)

        # Compute the prior and posterior standard deviation from the covariance matrix
        std[gn.PRIOR][gn.STATE0]       = np.sqrt(cov[gn.PRIOR][gn.STATE0].diagonal())
        std[gn.PRIOR][gn.EXTFORCE]     = np.sqrt(cov[gn.PRIOR][gn.EXTFORCE].diagonal())
        std[gn.POSTERIOR][gn.STATE0]   = np.sqrt(cov[gn.POSTERIOR][gn.STATE0].diagonal())
        std[gn.POSTERIOR][gn.EXTFORCE] = np.sqrt(cov[gn.POSTERIOR][gn.EXTFORCE].diagonal())

        # Get the log likelihood and store it in globdat
        model.take_action(gpact.GETLOGLIKELIHOOD, params, globdat)
        globdat[gn.GP][gn.LOGLIKELIHOOD] = params[gppn.LOGLIKELIHOOD]

        # Get the samples
        self.take_sample_action(samples, globdat)

        # Get the fields belonging to the samples
        self.take_sample_field_action(samples, globdat)

        return output

    def shutdown(self, globdat):
        pass


def declare(factory):
    factory.declare_module('GPSolver', GPExactModule)
