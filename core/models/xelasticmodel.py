import numpy as np

from jive.fem.names import Actions as act
from jive.fem.names import ParamNames as pn
from jive.fem.names import GlobNames as gn
from core.models.elasticmodel import ElasticModel
import jive.util.proputils as pu

ELEMENTS = 'elements'
YOUNG = 'young'
RHO = 'rho'
THICKNESS = 'thickness'
POISSON = 'poisson'
SHAPE = 'shape'
INTSCHEME = 'intScheme'
STATE = 'state'
TYPE = 'type'
DOFTYPES = ['dx', 'dy', 'dz']
PE_STATE = 'plane_strain'
PS_STATE = 'plane_stress'


class XElasticModel(ElasticModel):
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
            print('XElasticModel taking action', action)

    def configure(self, props, globdat):
        # This function gets only the core values from props

        # Get basic parameter values
        self._young = pu.soft_cast(props[YOUNG], float)
        self._rho = pu.soft_cast(props.get(RHO,0), float)

        if globdat[gn.MESHRANK] > 1:
            self._poisson = pu.soft_cast(props[POISSON], float)
        if globdat[gn.MESHRANK] == 2:
            self._thickness = pu.soft_cast(props.get(THICKNESS,1), float)
            self._state = props[STATE]
            if not self._state in (PE_STATE, PS_STATE):
                raise RuntimeError('ElasticModel: state in 2d should be plane_strain or plane_stress')

        # Get shape and element info
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][TYPE], props[SHAPE][INTSCHEME])
        egroup = globdat[gn.EGROUPS][props[ELEMENTS]]
        self._elems = egroup.get_elements()
        self._ielems = egroup.get_indices()
        self._nodes = self._elems.get_nodes()

        # Make sure the shape rank and mesh rank are identitcal
        if self._shape.global_rank() != globdat[gn.MESHRANK]:
            raise RuntimeError('ElasticModel: Shape rank must agree with mesh rank')

        # The rest of the configuration happens in configure_noprops
        self._configure_noprops(globdat)

    def _get_unit_mass_matrix(self, params, globdat):

        # Swap out rho for 1
        rho_ = self._rho
        self._rho = 1.

        # Swap out thickness for 1 if in 2d
        if self._rank == 2:
            thickness_ = self._thickness
            self._thickness = 1.

        # Get the mass matrix
        self._get_mass_matrix(params, globdat)

        # Restore the original rho value
        self._rho = rho_

        # Restore the original thickness value if in 2d
        if self._rank == 2:
            self._thickness = thickness_


def declare(factory):
    factory.declare_model('XElastic', XElasticModel)
