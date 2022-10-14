import numpy as np

from jive.fem.names import GlobNames as gn
from jive.fem.names import GPParamNames as gppn
from jive.fem.names import GPActions as gpact

from jive.app.module import Module

NSAMPLE = 'nsample'
SEED = 'seed'

class GPSamplerModule(Module):

    def init(self, props, globdat):
        myprops = props[self._name]
        self._nsample = int(myprops.get(NSAMPLE,1))
        self._seed = eval(myprops.get(SEED,'None'))
        if not self._seed is None:
            self._seed = int(self._seed)

    def run(self, globdat):
        model = globdat[gn.MODEL]

        # Define a dictionary for the settings of u
        u_params = {}
        u_params[gppn.FIELD] = 'u'
        u_params[gppn.NSAMPLE] = self._nsample
        u_params[gppn.RNG] = np.random.default_rng(self._seed)

        # Define a dictionary for the settings of f
        f_params = {}
        f_params[gppn.FIELD] = 'f'
        f_params[gppn.NSAMPLE] = self._nsample
        f_params[gppn.RNG] = np.random.default_rng(self._seed)

        # Take the appropriate actions for u
        model.take_action(gpact.GETPRIORSAMPLES, u_params, globdat)
        model.take_action(gpact.GETPOSTERIORSAMPLES, u_params, globdat)

        # Take the appropriate actions for f
        model.take_action(gpact.GETPRIORSAMPLES, f_params, globdat)
        model.take_action(gpact.GETPOSTERIORSAMPLES, f_params, globdat)

        globdat['samples_f_prior'] = f_params[gppn.PRIORSAMPLES]
        globdat['samples_u_prior'] = u_params[gppn.PRIORSAMPLES]
        globdat['samples_f_post'] = f_params[gppn.POSTERIORSAMPLES]
        globdat['samples_u_post'] = u_params[gppn.POSTERIORSAMPLES]

        return 'exit'

    def shutdown(self, globdat):
        pass

    def __solve(self, globdat):
        pass


def declare(factory):
    factory.declare_module('GPSampler', GPSamplerModule)
