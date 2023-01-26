import numpy as np

from jive.fem.names import Actions    as act
from jive.fem.names import ParamNames as pn
from jive.fem.names import GlobNames  as gn
from poissonmodel import PoissonModel
import jive.util.proputils as pu

ELEMENTS = 'elements'
KAPPA = 'kappa'
RHO = 'rho'
LOAD = 'q'
SHAPE = 'shape'
TYPE = 'type'
INTSCHEME = 'intScheme'
DOFTYPES = ['u']


class XPoissonModel(PoissonModel):
    def take_action(self, action, params, globdat):
        showmsg = True

        # Refer the core actions to the parent class
        super().take_action(action, params, globdat)

        # Add extended actions below
        if action == act.GETUNITMATRIX2:
            self._get_unit_mass_matrix(params, globdat)
        else:
            showmsg = False

        if showmsg:
            print('XPoissonModel taking action', action)

    def configure(self, props, globdat):
        # This function gets only the core values from props

        # Get basic parameter values
        self._kappa = pu.soft_cast(props[KAPPA], float)
        self._rho = pu.soft_cast(props.get(RHO,0), float)
        self._q = pu.soft_cast(props.get(LOAD,0), float)

        # Get shape and element info
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][TYPE], props[SHAPE][INTSCHEME])
        egroup = globdat[gn.EGROUPS][props[ELEMENTS]]
        self._elems = egroup.get_elements()
        self._ielems = egroup.get_indices()
        self._nodes = self._elems.get_nodes()

        # Make sure the shape rank and mesh rank are identitcal
        if self._shape.global_rank() != globdat[gn.MESHRANK]:
            raise RuntimeError('PoissonModel: Shape rank must agree with mesh rank')

        # The rest of the configuration happens in configure_noprops
        self._configure_noprops(globdat)

    def _get_unit_mass_matrix(self, params, globdat):

        # Swap out rho for 1
        rho_ = self._rho
        self._rho = 1.

        # Get the mass matrix
        self._get_mass_matrix(params, globdat)

        # Restore the original rho value
        self._rho = rho_


def declare(factory):
    factory.declare_model('XPoisson', XPoissonModel)
