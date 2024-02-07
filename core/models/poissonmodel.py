import numpy as np

from jive.fem.names import Actions as act
from jive.fem.names import ParamNames as pn
from jive.fem.names import GlobNames as gn
from jive.fem.names import PropNames as prn
from jive.model.model import Model

ELEMENTS = "elements"
KAPPA = "kappa"
RHO = "rho"
SHAPE = "shape"
TYPE = "type"
INTSCHEME = "intScheme"
DOFTYPES = ["u"]


class PoissonModel(Model):
    def take_action(self, action, params, globdat):
        showmsg = True

        if action == act.GETMATRIX0:
            self._get_matrix(params, globdat)
        elif action == act.GETMATRIX2:
            self._get_mass_matrix(params, globdat)
        else:
            showmsg = False

        verbose = params.get(pn.VERBOSE, True)

        if showmsg and verbose:
            print("PoissonModel taking action", action)

    def configure(self, props, globdat):
        # This function gets only the core values from props

        # Get basic parameter values
        self._kappa = float(props[KAPPA])
        self._rho = float(props.get(RHO, 0))

        # Get shape and element info
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(
            props[SHAPE][TYPE], props[SHAPE][INTSCHEME]
        )
        egroup = globdat[gn.EGROUPS][props[ELEMENTS]]
        self._elems = egroup.get_elements()
        self._ielems = egroup.get_indices()
        self._nodes = self._elems.get_nodes()

        # Make sure the shape rank and mesh rank are identitcal
        if self._shape.global_rank() != globdat[gn.MESHRANK]:
            raise RuntimeError("PoissonModel: Shape rank must agree with mesh rank")

        # The rest of the configuration happens in configure_noprops
        self._configure_noprops(globdat)

    def _configure_noprops(self, globdat):
        # This function gets additional info from self and globdat
        # It has been split off from configure() to allow it to be used in inherited classes as well.

        # Get basic dimensionality info
        self._rank = self._shape.global_rank()
        self._ipcount = self._shape.ipoint_count()
        self._dofcount = len(DOFTYPES) * self._shape.node_count()

        # Create a new dof for every node and dof type
        for doftype in DOFTYPES[0 : self._rank]:
            globdat[gn.DOFSPACE].add_type(doftype)
            for inode in self._elems.get_unique_nodes_of(self._ielems):
                globdat[gn.DOFSPACE].add_dof(inode, doftype)

    def _get_matrix(self, params, globdat):
        for ielem in self._ielems:
            # Get the nodal coordinates of each element
            inodes = self._elems.get_elem_nodes(ielem)
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES)
            coords = self._nodes.get_some_coords(inodes)

            # Get the gradients, weights and coordinates of each integration point
            grads, weights = self._shape.get_shape_gradients(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            # Reset the element stiffness matrix
            elmat = np.zeros((self._dofcount, self._dofcount))

            for ip in range(self._ipcount):
                # Get the B and D matrices for each integration point
                B = self._get_B_matrix(grads[:, :, ip])
                D = self._get_D_matrix(ipcoords[:, ip])

                # Compute the element stiffness matrix
                elmat += weights[ip] * np.matmul(np.transpose(B), np.matmul(D, B))

            # Add the element stiffness matrix to the global stiffness matrix
            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elmat

    def _get_mass_matrix(self, params, globdat):
        unit_matrix = params.get(pn.UNITMATRIX, False)

        if unit_matrix:
            M = np.identity(1)

        for ielem in self._ielems:
            # Get the nodal coordinates of each element
            inodes = self._elems.get_elem_nodes(ielem)
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0 : self._rank])
            coords = self._nodes.get_some_coords(inodes)

            # Get the shape functions, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            weights = self._shape.get_integration_weights(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            # Reset the element mass matrix
            elmat = np.zeros((self._dofcount, self._dofcount))

            for ip in range(self._ipcount):
                # Get the N and M matrices for each integration point
                N = self._get_N_matrix(sfuncs[:, ip])

                if not unit_matrix:
                    M = self._get_M_matrix(ipcoords[:, ip])

                # Compute the element mass matrix
                elmat += weights[ip] * np.matmul(np.transpose(N), np.matmul(M, N))

            # Add the element mass matrix to the global mass matrix
            params[pn.MATRIX2][np.ix_(idofs, idofs)] += elmat

    def _get_N_matrix(self, sfuncs):
        N = np.zeros((1, self._dofcount))
        N[0, :] = sfuncs.transpose()
        return N

    def _get_B_matrix(self, grads):
        B = np.zeros((self._rank, self._dofcount))
        B[:, :] = grads.transpose()
        return B

    def _get_D_matrix(self, ipcoords):
        D = self._kappa * np.identity(self._rank)
        return D

    def _get_M_matrix(self, ipcoords):
        M = np.array([[self._rho]])
        return M


def declare(factory):
    factory.declare_model("Poisson", PoissonModel)
