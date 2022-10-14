import numpy as np
import scipy.sparse as spsp

from jive.fem.names import GlobNames as gn
from jive.fem.names import ParamNames as pn
from jive.fem.names import Actions as act
from jive.fem.names import GPActions as gpact
from jive.fem.names import GPParamNames as gppn

from jive.solver.solvermodule import SolverModule
from jive.solver.constrainer import Constrainer

GETUNITMASSMATRIX = 'getUnitMassMatrix'
GETFORCERESULTS = 'getForceResults'
GETFULLCOVARIANCE = 'getFullCovariance'

class GPSolverModule(SolverModule):

    def init(self, props, globdat):

        # Initialize solvermodule first
        super().init(props, globdat)

        myprops = props[self._name]
        self._get_unit_mass_matrix = bool(eval(myprops.get(GETUNITMASSMATRIX, 'False')))
        self._get_force_results = bool(eval(myprops.get(GETFORCERESULTS, 'False')))
        self._get_full_covariance = bool(eval(myprops.get(GETFULLCOVARIANCE, 'False')))

    def run(self, globdat):

        # Run solvermodule first
        output = super().run(globdat)

        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[gn.MODEL]

        # Optionally get the mass matrix
        if self._get_unit_mass_matrix:
            M = spsp.csr_array((dc, dc))
            params = {
                pn.MATRIX2: M,
                pn.UNITMATRIX: True
            }
            model.take_action(act.GETMATRIX2, params, globdat)

            # Optionally store mass matrix in Globdat
            if self._store_matrix:
                globdat[gn.MATRIX2] = M

        # Configure the GP based on the fine FEM results
        model.take_action(gpact.CONFIGUREFEM, params, globdat)
        model.take_action(gpact.CONFIGUREPRIOR, params, globdat)

        if self._get_force_results:
            fields = ['u', 'f']
        else:
            fields = ['u']

        for field in fields:

            # Define a dictionary for the settings of u
            params = {}
            params[gppn.FIELD] = field
            params[gppn.FULLCOVARIANCE] = self._get_full_covariance

            # Take the appropriate actions for u
            model.take_action(gpact.GETPRIORMEAN, params, globdat)
            model.take_action(gpact.GETPOSTERIORMEAN, params, globdat)
            model.take_action(gpact.GETPRIORCOVARIANCE, params, globdat)
            model.take_action(gpact.GETPOSTERIORCOVARIANCE, params, globdat)

            # Get the log likelihood
            model.take_action(gpact.GETLOGLIKELIHOOD, params, globdat)

            # Optionally store stiffness matrix in Globdat
            if self._store_matrix:
                globdat[field+'_prior'] = params[gppn.PRIORMEAN]
                globdat[field+'_post'] = params[gppn.POSTERIORMEAN]
                globdat['var_'+field+'_prior'] = params[gppn.PRIORCOVARIANCE]
                globdat['var_'+field+'_post'] = params[gppn.POSTERIORCOVARIANCE]
                globdat['logLikelihood'] = params[gppn.LOGLIKELIHOOD]

        return output

    def shutdown(self, globdat):
        pass

    def __solve(self, globdat):
        pass


def declare(factory):
    factory.declare_module('GPSolver', GPSolverModule)
