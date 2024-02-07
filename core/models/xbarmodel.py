import numpy as np

from jive.fem.names import Actions as act
from jive.fem.names import ParamNames as pn
from jive.fem.names import GlobNames as gn
from core.models.barmodel import BarModel
import jive.util.proputils as pu

ELEMENTS = 'elements'
EA = 'EA'
k = 'k'
RHOA = 'rhoA'
SHAPE = 'shape'
TYPE = 'type'
INTSCHEME = 'intScheme'
DOFTYPES = ['dx']
LOAD = 'q'

class XBarModel(BarModel):
    def take_action(self, action, params, globdat):
        showmsg = True

        # Refer the core actions to the parent class
        super().take_action(action, params, globdat)

        # Add extended actions below
        if action == act.GETUNITMATRIX2:
            self._get_unit_mass_matrix(params, globdat)
        else:
            showmsg = False

        verbose = params.get(pn.VERBOSE, True)

        if showmsg and verbose:
            print('XBarModel taking action', action)

    def configure(self, props, globdat):
        # This function gets only the core values from props

        # Get basic parameter values
        self._EA = pu.soft_cast(props[EA], float)
        self._k = pu.soft_cast(props.get(k, 0), float)
        self._rhoA = pu.soft_cast(props.get(RHOA, 0), float)
        self._q = pu.soft_cast(props.get(LOAD, 0), float)

        # Get shape and element info
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][TYPE], props[SHAPE][INTSCHEME])
        egroup = globdat[gn.EGROUPS][props[ELEMENTS]]
        self._elems = egroup.get_elements()
        self._ielems = egroup.get_indices()
        self._nodes = self._elems.get_nodes()

        # The rest of the configuration happens in configure_noprops
        self._configure_noprops(globdat)

    def _get_unit_mass_matrix(self, params, globdat):

        # Swap out rhoA for 1
        rhoA_ = self._rhoA
        self._rhoA = 1.

        # Get the mass matrix
        self._get_mass_matrix(params, globdat)

        # Restore the original rhoA value
        self._rhoA = rhoA_


def declare(factory):
    factory.declare_model('XBar', XBarModel)
