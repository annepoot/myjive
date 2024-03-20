import numpy as np

from myjive.names import GlobNames as gn
from myjive.model.model import Model
import myjive.util.proputils as pu
from myjive.util.proputils import mandatory_argument, mandatory_dict, mandatory_list

TYPE = "type"
INTSCHEME = "intScheme"

__all__ = ["LoadModel"]


class LoadModel(Model):
    def GETEXTFORCE(self, f_ext, globdat, **kwargs):
        f_ext = self._get_body_force(f_ext, globdat, **kwargs)
        return f_ext

    def configure(self, globdat, **props):

        # Get props
        shapeprops = mandatory_dict(
            self, props, "shape", mandatory_keys=[TYPE, INTSCHEME]
        )
        elements = mandatory_argument(self, props, "elements")
        self._doftypes = mandatory_list(self, props, "dofs")
        self._loads = mandatory_list(self, props, "values")

        # Get shape and element info
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(
            shapeprops[TYPE], shapeprops[INTSCHEME]
        )
        egroup = globdat[gn.EGROUPS][elements]
        self._elems = egroup.get_elements()
        self._ielems = egroup.get_indices()
        self._nodes = self._elems.get_nodes()

        # Make sure the shape rank and mesh rank are identitcal
        if self._shape.global_rank() != globdat[gn.MESHRANK]:
            raise RuntimeError("LoadModel: Shape rank must agree with mesh rank")

        # Get basic dimensionality info
        self._rank = self._shape.global_rank()
        self._ipcount = self._shape.ipoint_count()

        # Get the relevant dofs
        for i, load in enumerate(self._loads):
            self._loads[i] = pu.soft_cast(load, float)

        # Make sure the doftypes and loads have the same size
        if len(self._doftypes) != len(self._loads):
            raise ValueError("LoadModel: dofs and values must have the same size")

        # Get the dofcount (of only the relevant dofs!)
        self._loadcount = len(self._doftypes)
        self._dofcount = self._loadcount * self._shape.node_count()

    def _get_body_force(self, f_ext, globdat):
        if f_ext is None:
            f_ext = np.zeros(globdat[gn.DOFSPACE].dof_count())

        for ielem in self._ielems:
            # Get the nodal coordinates of each element
            inodes = self._elems.get_elem_nodes(ielem)
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, self._doftypes)
            coords = self._nodes.get_some_coords(inodes)

            # Get the shape functions, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            weights = self._shape.get_integration_weights(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            # Reset the element force vector
            elfor = np.zeros(self._dofcount)

            for ip in range(self._ipcount):
                # Get the N matrix and b vector for each integration point
                N = self._get_N_matrix(sfuncs[:, ip])
                b = self._get_b_vector(ipcoords[:, ip])

                # Compute the element force vector
                elfor += weights[ip] * np.matmul(np.transpose(N), b)

            # Add the element force vector to the global force vector
            f_ext[idofs] += elfor

        return f_ext

    def _get_N_matrix(self, sfuncs):
        N = np.zeros((self._loadcount, self._dofcount))
        for i in range(self._loadcount):
            N[i, i :: self._loadcount] = sfuncs.transpose()
        return N

    def _get_b_vector(self, ipcoords):
        b = np.zeros((self._loadcount))
        for i in range(self._loadcount):
            b[i] = pu.evaluate(self._loads[i], ipcoords, self._rank)
        return b
