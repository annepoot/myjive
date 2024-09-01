import numpy as np

from myjive.names import GlobNames as gn
from myjive.model.model import Model
from myjive.util.proputils import check_dict

TYPE = "type"
INTSCHEME = "intScheme"
DOFTYPES = ["phi", "dy"]

__all__ = ["TimoshenkoModel"]


class TimoshenkoModel(Model):
    def GETMATRIX0(self, K, globdat, **kwargs):
        K = self._get_matrix(K, globdat, **kwargs)
        return K

    @Model.save_config
    def configure(self, globdat, *, shape, elements, EI, GAs):
        # Validate input arguments
        check_dict(self, shape, [TYPE, INTSCHEME])
        self._EI = EI
        self._GAs = GAs

        # Get shape and element info
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(shape[TYPE], shape[INTSCHEME])
        egroup = globdat[gn.EGROUPS][elements]
        self._elems = egroup.get_elements()
        self._ielems = egroup.get_indices()
        self._nodes = self._elems.get_nodes()

        self._ipcount = self._shape.ipoint_count()
        self._dofcount = 2 * self._shape.node_count()

        for doftype in DOFTYPES:
            globdat[gn.DOFSPACE].add_type(doftype)
            for node in self._elems.get_unique_nodes_of(self._ielems):
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
                B_theta[:, 0::2] = grads[ip]
                B_v[:, 1::2] = grads[ip]
                N_theta[:, 0::2] = sfuncs[ip]
                N_v[:, 1::2] = sfuncs[ip]

                elmat += weights[ip] * (
                    (B_theta.T * self._EI) @ B_theta
                    + (N_theta.T * self._GAs) @ N_theta
                    - (N_theta.T * self._GAs) @ B_v
                    - (B_v.T * self._GAs) @ N_theta
                    + (B_v.T * self._GAs) @ B_v
                )

            K[np.ix_(idofs, idofs)] += elmat

        return K
