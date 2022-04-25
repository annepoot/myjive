import numpy as np

from names import GlobNames as gn
from names import GPParamNames as gppn
from names import PropNames as prn
from names import GPActions as gpact

from module import Module

NSAMPLE = 'nsample'

class SamplerModule(Module):

    def init(self, props, globdat):
        myprops = props[self._name]
        self._nsample = int(myprops.get(NSAMPLE,1))

        # Get the appropriate model for this module
        self._modelname = myprops.get(gn.MODEL, gn.MODEL)
        modelprops = props[self._modelname]
        modelfac = globdat[gn.MODELFACTORY]

    def run(self, globdat):
        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[self._modelname]

        # Configure the model again, to make sure K, M and f are stored there as well
        model.configure_fem(globdat)

        # Define a dictionary for the settings of u
        u_params = {}
        u_params[gppn.FIELD] = 'u'
        u_params[gppn.NSAMPLE] = self._nsample

        # Define a dictionary for the settings of f
        f_params = {}
        f_params[gppn.FIELD] = 'f'
        f_params[gppn.NSAMPLE] = self._nsample

        # Take the appropriate actions for u
        model.take_action(gpact.GETPRIORSAMPLES, u_params, globdat)
        model.take_action(gpact.GETPOSTERIORSAMPLES, u_params, globdat)

        # Take the appropriate actions for f
        model.take_action(gpact.GETPRIORSAMPLES, f_params, globdat)
        model.take_action(gpact.GETPOSTERIORSAMPLES, f_params, globdat)
        model.take_action(gpact.GETPRIORCOVARIANCE, f_params, globdat)
        model.take_action(gpact.GETPOSTERIORCOVARIANCE, f_params, globdat)

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
    factory.declare_module('Sampler', SamplerModule)
