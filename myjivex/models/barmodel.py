import numpy as np

from myjive.names import GlobNames as gn
from myjive.model.model import Model
from myjive.util.proputils import mandatory_argument, mandatory_dict, optional_argument

TYPE = "type"
INTSCHEME = "intScheme"
DOFTYPES = ["dx"]

__all__ = ["BarModel"]


class BarModel(Model):
    def GETMATRIX0(self, K, globdat, **kwargs):
        K = self._get_matrix(K, globdat, **kwargs)
        return K

    def GETMATRIX2(self, M, globdat, **kwargs):
        M = self._get_mass_matrix(M, globdat, **kwargs)
        return M

    def configure(self, globdat, **props):

        # Get props
        shapeprops = mandatory_dict(
            self, props, "shape", mandatory_keys=[TYPE, INTSCHEME]
        )
        elements = mandatory_argument(self, props, "elements")
        self._EA = mandatory_argument(self, props, "EA")
        self._k = optional_argument(self, props, "k", default=0.0)
        self._rhoA = optional_argument(self, props, "rhoA", default=0.0)

        # Get shape and element info
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(
            shapeprops[TYPE], shapeprops[INTSCHEME]
        )
        egroup = globdat[gn.EGROUPS][elements]
        self._elems = egroup.get_elements()
        self._ielems = egroup.get_indices()
        self._nodes = self._elems.get_nodes()

        # The rest of the configuration happens in configure_noprops
        self._configure_noprops(globdat)

    def _configure_noprops(self, globdat):
        # This function gets additional info from self and globdat
        # It has been split off from configure() to allow it to be used in inherited classes as well.

        # Get basic dimensionality info
        self._rank = self._shape.global_rank()
        self._ipcount = self._shape.ipoint_count()
        self._nodecount = self._shape.node_count()
        self._dofcount = self._rank * self._nodecount
        self._strcount = self._rank * (self._rank + 1) // 2  # 1-->1, 2-->3, 3-->6

        # Create a new dof for every node and dof type
        nodes = np.unique([node for elem in self._elems for node in elem.get_nodes()])
        for doftype in DOFTYPES[0 : self._rank]:
            globdat[gn.DOFSPACE].add_type(doftype)
            for node in nodes:
                globdat[gn.DOFSPACE].add_dof(node, doftype)

    def _get_matrix(self, K, globdat):
        if K is None:
            dc = globdat[gn.DOFSPACE].dof_count()
            K = np.zeros((dc, dc))

        for ielem in self._ielems:
            # Get the nodal coordinates of each element
            inodes = self._elems.get_elem_nodes(ielem)
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0 : self._rank])
            coords = self._nodes.get_some_coords(inodes)

            # Get the shape functions, gradients, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            # Reset the element stiffness matrix
            elmat = np.zeros((self._dofcount, self._dofcount))

            for ip in range(self._ipcount):
                # Get the B, N, D and K matrices for each integration point
                B_elem = self._get_B_matrix(grads[:, :, ip])
                N_elem = self._get_N_matrix(sfuncs[:, ip])
                D_elem = self._get_D_matrix(ipcoords[:, ip])
                K_elem = self._get_K_matrix(ipcoords[:, ip])

                # Compute the element stiffness matrix
                elmat += weights[ip] * (
                    np.matmul(np.transpose(B_elem), np.matmul(D_elem, B_elem))
                    + np.matmul(np.transpose(N_elem), np.matmul(K_elem, N_elem))
                )

            # Add the element stiffness matrix to the global stiffness matrix
            K[np.ix_(idofs, idofs)] += elmat

        return K

    def _get_mass_matrix(self, M, globdat, unit_matrix=False):
        if M is None:
            dc = globdat[gn.DOFSPACE].dof_count()
            M = np.zeros((dc, dc))

        if unit_matrix:
            M_elem = np.identity(1)

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
                N_elem = self._get_N_matrix(sfuncs[:, ip])

                if not unit_matrix:
                    M_elem = self._get_M_matrix(ipcoords[:, ip])

                # Compute the element mass matrix
                elmat += weights[ip] * np.matmul(
                    np.transpose(N_elem), np.matmul(M_elem, N_elem)
                )

            # Add the element mass matrix to the global mass matrix
            M[np.ix_(idofs, idofs)] += elmat

        return M

    def _get_N_matrix(self, sfuncs):
        N = np.zeros((1, self._nodecount))
        N[0, :] = sfuncs.transpose()
        return N

    def _get_B_matrix(self, grads):
        B = np.zeros((1, self._nodecount))
        B[0, :] = grads.transpose()
        return B

    def _get_D_matrix(self, ipcoords):
        D = np.array([[self._EA]])
        return D

    def _get_K_matrix(self, ipcoords):
        K = np.array([[self._k]])
        return K

    def _get_M_matrix(self, ipcoords):
        M = np.array([[self._rhoA]])
        return M
