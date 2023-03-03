import numpy as np

from jive.fem.names import GlobNames as gn

from gpmodule import GPModule

GETUNITMASSMATRIX = 'getUnitMassMatrix'
POSTPROJECT = 'postproject'
PRIORMEAN = 'priorMean'
NSAMPLE = 'nsample'
SEED = 'seed'

class GPSamplerModule(GPModule):

    def run(self, globdat):

        # Run gpmodule (and solvermodule) first
        output = super().run(globdat)

        # Get the sub-dictionaries of the mean, std and samples
        mean = globdat[gn.GP][gn.MEAN]
        std = globdat[gn.GP][gn.STD]
        samples = globdat[gn.GP][gn.SAMPLES]

        # Get the samples
        self.take_sample_action(samples, globdat)

        # Compute the prior and posterior mean from the samples
        mean[gn.PRIOR][gn.STATE0]       = np.mean(samples[gn.PRIOR][gn.STATE0], axis=1)
        mean[gn.PRIOR][gn.EXTFORCE]     = np.mean(samples[gn.PRIOR][gn.EXTFORCE], axis=1)
        mean[gn.POSTERIOR][gn.STATE0]   = np.mean(samples[gn.POSTERIOR][gn.STATE0], axis=1)
        mean[gn.POSTERIOR][gn.EXTFORCE] = np.mean(samples[gn.POSTERIOR][gn.EXTFORCE], axis=1)

        # Compute the prior and posterior standard deviation from the samples
        std[gn.PRIOR][gn.STATE0]       = np.std(samples[gn.PRIOR][gn.STATE0], axis=1)
        std[gn.PRIOR][gn.EXTFORCE]     = np.std(samples[gn.PRIOR][gn.EXTFORCE], axis=1)
        std[gn.POSTERIOR][gn.STATE0]   = np.std(samples[gn.POSTERIOR][gn.STATE0], axis=1)
        std[gn.POSTERIOR][gn.EXTFORCE] = np.std(samples[gn.POSTERIOR][gn.EXTFORCE], axis=1)

        return output

    def shutdown(self, globdat):
        pass


def declare(factory):
    factory.declare_module('GPSampler', GPSamplerModule)
