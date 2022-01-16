import numpy as np

from names import Actions    as act
from names import ParamNames as pn
from names import GlobNames  as gn
from names import PropNames  as prn
from model import *

ELEMENTS = 'elements'
KAPPA = 'kappa'
SHAPE = 'shape'
INTSCHEME = 'intScheme'
DOFTYPES = ['u']


class PoissonModel(Model):
    def take_action(self, action, params, globdat):
        print('PoissonModel taking action', action)
        lim = 2
        if action == act.GETMATRIX0:
            self._get_matrix(params, globdat)

    def configure(self, props, globdat):
        self._kappa = float(props[KAPPA])
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][prn.TYPE], props[SHAPE][INTSCHEME])
        egroup = globdat[gn.EGROUPS][props[ELEMENTS]]
        self._elems = [globdat[gn.ESET][e] for e in egroup]

        if self._shape.global_rank() != globdat[gn.MESHRANK]:
            raise RuntimeError('PoissonModel: Shape rank must agree with mesh rank')

        self._rank = self._shape.global_rank()
        self._ipcount = self._shape.ipoint_count()
        self._dofcount = len(DOFTYPES) * self._shape.node_count()

        nodes = np.unique([node for elem in self._elems for node in elem.get_nodes()])
        for doftype in DOFTYPES:
            globdat[gn.DOFSPACE].add_type(doftype)
            for node in nodes:
                globdat[gn.DOFSPACE].add_dof(node, doftype)

    def _get_matrix(self, params, globdat):
        kappa = np.eye(self._rank) * self._kappa
        for elem in self._elems:
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES)
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:self._rank, :]
            grads, weights = self._shape.get_shape_gradients(coords)

            elmat = np.zeros((self._dofcount, self._dofcount))
            for ip in range(self._ipcount):
                B = np.transpose(grads[:, :, ip])
                elmat += weights[ip] * np.matmul(np.transpose(B), np.matmul(kappa, B))

            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elmat


def declare(factory):
    factory.declare_model('Poisson', PoissonModel)
