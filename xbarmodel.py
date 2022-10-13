import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn
from barmodel import BarModel
import proputils as pu

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
        print('XBarModel taking action', action)

        # Refer the core actions to the parent class
        super().take_action(action, params, globdat)

        # Add extended actions below
        if action == act.GETEXTFORCE:
            self._get_body_force(params, globdat)
        elif action == act.GETUNITMATRIX2:
            self._get_unit_mass_matrix(params, globdat)

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
        EA_ = pu.evaluate(self._EA, ipcoords, 1)
        D = np.array([[EA_]])
        return D

    def _get_K_matrix(self, ipcoords):
        k_ = pu.evaluate(self._k, ipcoords, 1)
        K = np.array([[k_]])
        return K

    def _get_M_matrix(self, ipcoords):
        rhoA_ = pu.evaluate(self._rhoA, ipcoords, 1)
        M = np.array([[rhoA_]])
        return M

    def _get_Q_vector(self, ipcoords):
        q_ = pu.evaluate(self._q, ipcoords, 1)
        Q = np.array([q_])
        return Q


def declare(factory):
    factory.declare_model('XBar', XBarModel)
