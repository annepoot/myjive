import numpy as np

from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act
from names import GPActions as gpact
from names import GPParamNames as gppn

from solvermodule import SolverModule
from constrainer import Constrainer

GETUNITMASSMATRIX = 'getUnitMassMatrix'

class GPSolverModule(SolverModule):

    def init(self, props, globdat):

        # Initialize solvermodule first
        super().init(props, globdat)

        myprops = props[self._name]
        self._get_unit_mass_matrix = bool(eval(myprops.get(GETUNITMASSMATRIX,'False')))

    def run(self, globdat):

        # Run solvermodule first
        output = super().run(globdat)

        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[gn.MODEL]

        # Optionally get the mass matrix
        if self._get_unit_mass_matrix:
            M = np.zeros((dc, dc))
            params = {pn.MATRIX2: M}
            model.take_action(act.GETUNITMATRIX2, params, globdat)

            # Optionally store mass matrix in Globdat
            if self._store_matrix:
                globdat[gn.MATRIX2] = M

        # Configure the GP based on the fine FEM results
        model.take_action(gpact.CONFIGUREFEM, params, globdat)

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
            globdat['var_f_prior'] = f_params[gppn.PRIORCOVARIANCE]
            globdat['var_f_post'] = f_params[gppn.POSTERIORCOVARIANCE]
            globdat['var_u_prior'] = u_params[gppn.PRIORCOVARIANCE]
            globdat['var_u_post'] = u_params[gppn.POSTERIORCOVARIANCE]

        return output

    def shutdown(self, globdat):
        pass

    def __solve(self, globdat):
        pass


def declare(factory):
    factory.declare_module('GPSolver', GPSolverModule)
