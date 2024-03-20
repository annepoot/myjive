import numpy as np

from myjive.names import GlobNames as gn
from myjive.model.model import Model
from myjive.util.proputils import mandatory_argument, mandatory_dict

TYPE = "type"
INTSCHEME = "intScheme"
DOFTYPES = ["phi", "dy"]

__all__ = ["TimoshenkoModel"]


class TimoshenkoModel(Model):
    def GETMATRIX0(self, K, globdat, **kwargs):
        K = self._get_matrix(K, globdat, **kwargs)
        return K

    def configure(self, globdat, **props):

        # Get props
        shapeprops = mandatory_dict(
            self, props, "shape", mandatory_keys=[TYPE, INTSCHEME]
        )
        elements = mandatory_argument(self, props, "elements")
        self._EI = mandatory_argument(self, props, "EI")
        self._GAs = mandatory_argument(self, props, "GAs")

        # Get shape and element info
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(
            shapeprops[TYPE], shapeprops[INTSCHEME]
        )

        egroup = globdat[gn.EGROUPS][elements]
        self._elems = egroup.get_elements()
        self._ielems = egroup.get_indices()
        self._nodes = self._elems.get_nodes()

        self._ipcount = self._shape.ipoint_count()
        self._dofcount = 2 * self._shape.node_count()

        nodes = np.unique([node for elem in self._elems for node in elem.get_nodes()])
        for doftype in DOFTYPES:
            globdat[gn.DOFSPACE].add_type(doftype)
            for node in nodes:
                globdat[gn.DOFSPACE].add_dof(node, doftype)

    def _get_matrix(self, K, globdat):
        if K is None:
            dc = globdat[gn.DOFSPACE].dof_count()
            K = np.zeros((dc, dc))

        for ielem in self._ielems:
            inodes = self._elems.get_elem_nodes(ielem)
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES)
            coords = self._nodes.get_some_coords(inodes)

            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords)

            elmat = np.zeros((4, 4))
            for ip in range(self._ipcount):
                B_theta = np.zeros((1, 4))
                N_theta = np.zeros((1, 4))
                B_v = np.zeros((1, 4))
                N_v = np.zeros((1, 4))
                B_theta[:, 0::2] = grads[:, :, ip].transpose()
                B_v[:, 1::2] = grads[:, :, ip].transpose()
                N_theta[:, 0::2] = sfuncs[:, ip].transpose()
                N_v[:, 1::2] = sfuncs[:, ip].transpose()

                elmat += weights[ip] * (
                    np.matmul(B_theta.transpose() * self._EI, B_theta)
                    + np.matmul(N_theta.transpose() * self._GAs, N_theta)
                    - np.matmul(N_theta.transpose() * self._GAs, B_v)
                    - np.matmul(B_v.transpose() * self._GAs, N_theta)
                    + np.matmul(B_v.transpose() * self._GAs, B_v)
                )

            K[np.ix_(idofs, idofs)] += elmat

        return K
