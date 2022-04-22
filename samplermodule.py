import numpy as np

from names import GlobNames as gn
from names import GPParamNames as gppn
from names import PropNames as prn
from names import GPActions as gpact

from module import Module

NOBS = 'nobs'
OBSNOISE = 'obsNoise'
ALPHA = 'alpha'
STOREMATRIX = 'storeMatrix'
RANDOMOBS = 'randomObs'

class SamplerModule(Module):

    def init(self, props, globdat):
        myprops = props[self._name]
        self._store_matrix = bool(eval(myprops.get(STOREMATRIX,'False')))
        self._dc = globdat[gn.DOFSPACE].dof_count()

        # Get the appropriate model for this module
        self._modelname = myprops.get(gn.MODEL, gn.MODEL)
        modelprops = props[self._modelname]
        modelfac = globdat[gn.MODELFACTORY]

    def run(self, globdat):
        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[self._modelname]

        # Configure the model again, to make sure K, M and f are stored there as well
        model.configure_fem(globdat)

        # take_action comments go here....

    def shutdown(self, globdat):
        pass

    def __solve(self, globdat):
        pass


def declare(factory):
    factory.declare_module('Sampler', SamplerModule)
