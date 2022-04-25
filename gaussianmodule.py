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

class GaussianModule(Module):

    def init(self, props, globdat):
        myprops = props[self._name]
        self._store_matrix = bool(eval(myprops.get(STOREMATRIX,'False')))

        # Get the appropriate model for this module
        self._modelname = myprops.get(gn.MODEL, gn.MODEL)
        modelprops = props[self._modelname]
        modelfac = globdat[gn.MODELFACTORY]

        # Initialize model
        print('GaussianModule: Creating model...')
        m = modelfac.get_model(modelprops[prn.TYPE], self._modelname)
        m.configure(modelprops, globdat)
        globdat[self._modelname] = m

    def run(self, globdat):
        model = globdat[self._modelname]

        # Configure the model again, to make sure K, M and f are stored there as well
        model.configure_fem(globdat)

        # Define a dictionary for the settings of u
        u_params = {}
        u_params[gppn.FIELD] = 'u'
        u_params[gppn.FULLCOVARIANCE] = False

        # Define a dictionary for the settings of f
        f_params = {}
        f_params[gppn.FIELD] = 'f'
        f_params[gppn.FULLCOVARIANCE] = False

        # Take the appropriate actions for u
        model.take_action(gpact.GETPRIORMEAN, u_params, globdat)
        model.take_action(gpact.GETPOSTERIORMEAN, u_params, globdat)
        model.take_action(gpact.GETPRIORCOVARIANCE, u_params, globdat)
        model.take_action(gpact.GETPOSTERIORCOVARIANCE, u_params, globdat)

        # Take the appropriate actions for f
        model.take_action(gpact.GETPRIORMEAN, f_params, globdat)
        model.take_action(gpact.GETPOSTERIORMEAN, f_params, globdat)
        model.take_action(gpact.GETPRIORCOVARIANCE, f_params, globdat)
        model.take_action(gpact.GETPOSTERIORCOVARIANCE, f_params, globdat)

        # Optionally store stiffness matrix in Globdat
        if ( self._store_matrix ):
            globdat['f_prior'] = f_params[gppn.PRIORMEAN]
            globdat['u_prior'] = u_params[gppn.PRIORMEAN]
            globdat['f_post'] = f_params[gppn.POSTERIORMEAN]
            globdat['u_post'] = u_params[gppn.POSTERIORMEAN]
            globdat['sigma_f_prior'] = f_params[gppn.PRIORCOVARIANCE]
            globdat['sigma_f_post'] = f_params[gppn.POSTERIORCOVARIANCE]
            globdat['sigma_u_prior'] = u_params[gppn.PRIORCOVARIANCE]
            globdat['sigma_u_post'] = u_params[gppn.POSTERIORCOVARIANCE]

        return 'exit'

    def shutdown(self, globdat):
        pass

    def __solve(self, globdat):
        pass


def declare(factory):
    factory.declare_module('Gaussian', GaussianModule)
