import numpy as np

from names import Actions    as act
from names import ParamNames as pn
from names import GlobNames  as gn
from names import PropNames  as prn
from poissonmodel import PoissonModel
import proputils as pu

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
        print('XElasticModel taking action', action)

        # Refer the core actions to the parent class
        super().take_action(action, params, globdat)

        # Add extended actions below
        if action == act.GETUNITMATRIX2:
            self._get_unit_mass_matrix(params, globdat)
        if action == act.GETEXTFORCE:
            self._get_body_force(params, globdat)

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

    def _get_body_force(self, params, globdat):

        for elem in self._elems:
            # Get the nodal coordinates of each element
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:self._rank, :]

            # Get the shape functions, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            weights = self._shape.get_integration_weights(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            # Reset the element force vector
            elfor = np.zeros(self._dofcount)

            for ip in range(self._ipcount):
                # Get the N matrix and Q vector for each integration point
                N = self._get_N_matrix(sfuncs[:,ip])
                Q = self._get_Q_vector(ipcoords[:,ip])

                # Compute the element force vector
                elfor += weights[ip] * np.matmul(np.transpose(N), Q)

            # Add the element force vector to the global force vector
            params[pn.EXTFORCE][idofs] += elfor


    def _get_D_matrix(self, ipcoords):
        kappa_ = pu.evaluate(self._kappa, {'x':ipcoords[0], 'y':ipcoords[1]})
        D = kappa_ * np.identity(self._rank)
        return D

    def _get_M_matrix(self, ipcoords):
        rho_ = pu.evaluate(self._rho, {'x':ipcoords[0], 'y':ipcoords[1]})
        M = np.array([[rho_]])
        return M

    def _get_Q_vector(self, ipcoords):
        q_ = pu.evaluate(self._q, {'x':ipcoords[0], 'y':ipcoords[1]})
        Q = np.array([q_])
        return Q

def declare(factory):
    factory.declare_model('XPoisson', XPoissonModel)
