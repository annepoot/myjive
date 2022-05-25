import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn
from names import PropNames as prn
from model import Model

ELEMENTS = 'elements'
EA = 'EA'
YOUNG = 'E'
AREA = 'A'
k = 'k'
RHO = 'rho'
SHAPE = 'shape'
INTSCHEME = 'intScheme'
DOFTYPES = ['dx']


class BarModel(Model):
    def take_action(self, action, params, globdat):
        print('BarModel taking action', action)

        if action == act.GETMATRIX0:
            self._get_matrix(params, globdat)
        elif action == act.GETMATRIX2:
            self._get_mass_matrix(params, globdat)

    def configure(self, props, globdat):

        if EA in props:
            self._E = float(props[EA])
            self._A = 1.0
        else:
            self._E = props[YOUNG]
            self._A = props[AREA]

        self._k = float(props.get(k, 0))
        self._rho = float(props.get(RHO, 0))

        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][prn.TYPE], props[SHAPE][INTSCHEME])
        egroup = globdat[gn.EGROUPS][props[ELEMENTS]]
        self._elems = [globdat[gn.ESET][e] for e in egroup]

        self._rank = self._shape.global_rank()
        self._ipcount = self._shape.ipoint_count()
        self._nodecount = self._shape.node_count()
        self._dofcount = self._rank * self._nodecount
        self._strcount = self._rank * (self._rank + 1) // 2

        nodes = np.unique([node for elem in self._elems for node in elem.get_nodes()])
        for doftype in DOFTYPES[0:self._rank]:
            globdat[gn.DOFSPACE].add_type(doftype)
            for node in nodes:
                globdat[gn.DOFSPACE].add_dof(node, doftype)

    def _get_matrix(self, params, globdat):
        D = np.array([[self._E * self._A]])
        K = np.array([[self._k]])

        for elem in self._elems:
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:self._rank, :]
            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords)

            elmat = np.zeros((self._dofcount, self._dofcount))
            for ip in range(self._ipcount):
                B = np.zeros((1, self._nodecount))
                N = np.zeros((1, self._nodecount))
                B = grads[:, :, ip].transpose()
                N[0, :] = sfuncs[:, ip].transpose()

                elmat += weights[ip] * (np.matmul(np.transpose(B), np.matmul(D, B))
                                        + np.matmul(np.transpose(N), np.matmul(K, N)))

            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elmat

    def _get_mass_matrix(self, params, globdat):
        M = np.array([[self._rho * self._A]])

        for elem in self._elems:
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:self._rank, :]
            sfuncs = self._shape.get_shape_functions()
            weights = self._shape.get_integration_weights(coords)

            elmat = np.zeros((self._dofcount, self._dofcount))
            for ip in range(self._ipcount):
                N = np.zeros((1, self._nodecount))
                N[0, :] = sfuncs[:, ip].transpose()

                elmat += weights[ip] * np.matmul(np.transpose(N), np.matmul(M, N))

            params[pn.MATRIX2][np.ix_(idofs, idofs)] += elmat


def declare(factory):
    factory.declare_model('Bar', BarModel)
