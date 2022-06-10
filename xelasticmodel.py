import numpy as np

from names import Actions as act
from names import GlobNames as gn
from elasticmodel import ElasticModel
import proputils as pu

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
        print('XElasticModel taking action', action)

        # Refer the core actions to the parent class
        super().take_action(action, params, globdat)

        # Add extended actions below
        if action == act.GETUNITMATRIX2:
            self._get_unit_mass_matrix(params, globdat)

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
        self._elems = [globdat[gn.ESET][e] for e in egroup]

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

    def _get_D_matrix(self, ipcoords):
        D = np.zeros((self._strcount, self._strcount))
        E = pu.evaluate(self._young, {'x':ipcoords[0], 'y':ipcoords[1]})

        if self._rank == 1:
            D[[0]] = E
            return D

        nu = pu.evaluate(self._poisson, {'x':ipcoords[0], 'y':ipcoords[1]})
        g = 0.5 * E / (1. + nu)

        if self._rank == 3:
            a = E * (1. - nu) / ((1. + nu) * (1. - 2. * nu))
            b = E * nu / ((1. + nu) * (1. - 2. * nu))
            D[:,:] = [
                [a, b, b, .0, .0, .0],
                [b, a, b, .0, .0, .0],
                [b, b, a, .0, .0, .0],
                [.0, .0, .0, g, .0, .0],
                [.0, .0, .0, .0, g, .0],
                [.0, .0, .0, .0, .0, g]
            ]

        elif self._rank == 2:
            if self._state == PE_STATE:
                a = E * (1. - nu) / ((1. + nu) * (1. - 2. * nu))
                b = E * nu / ((1. + nu) * (1. - 2. * nu))
                D[:,:] = [
                    [a, b, .0],
                    [b, a, .0],
                    [.0, .0, g]
                ]
            else:
                assert (self._state == PS_STATE)
                a = E / (1. - nu * nu)
                b = a * nu
                D[:,:] = [
                    [a, b, .0],
                    [b, a, .0],
                    [.0, .0, g]
                ]

            D *= self._thickness

        return D

    def _get_M_matrix(self, ipcoords):
        rho_ = pu.evaluate(self._rho, {'x':ipcoords[0], 'y':ipcoords[1]})
        M = rho_ * np.identity(self._rank)

        if self._rank == 2:
            M *= self._thickness

        return M

    def _get_b_vector(self, ipcoords):
        rho_ = pu.evaluate(self._rho, {'x':ipcoords[0], 'y':ipcoords[1]})

        if self._rank == 3:
            gravity = np.array([0, -1, 0])
            b = rho_ * gravity

        elif self._rank == 2:
            gravity = np.array([0, -1])
            b = rho_ * self._thickness * gravity

        else:
            raise RuntimeError('ElasticModel: self weight can only be computed in 2d or 3d')

        return b


def declare(factory):
    factory.declare_model('XElastic', XElasticModel)
