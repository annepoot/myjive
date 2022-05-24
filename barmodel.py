import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn
from names import PropNames as prn
from model import Model

ELEMENTS = 'elements'
EA = 'EA'
E = 'E'
A = 'A'
k = 'k'
q = 'q'
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
        elif action == act.GETEXTFORCE:
            self._get_body_force(params, globdat)

    def configure(self, props, globdat):
        if EA in props:
            self._E = props[EA]
            self._A = 1.0
        else:
            self._E = props[E]
            self._A = props[A]

        try:
            self._A = float(self._A)
        except:
            pass

        try:
            self._E = float(self._E)
        except:
            pass

        self._k = float(props.get(k, 0))
        self._rho = float(props.get(RHO,0))
        self._q = float(props.get(q, 0))

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
        const_E = type(self._E) is float
        const_A = type(self._A) is float
        const_k = type(self._k) is float

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

                if not const_E or not const_A or not const_k:
                    x = self._shape.get_global_point(self._shape.get_integration_points()[:,ip], coords)[0]

                E_ = self._E if const_E else eval(self._E)
                A_ = self._A if const_A else eval(self._A)
                k_ = self._k if const_k else eval(self._k)

                D = np.array([[E_ * A_]])
                K = np.array([[k_]])

                elmat += weights[ip] * (np.matmul(np.transpose(B), np.matmul(D, B))
                                        + np.matmul(np.transpose(N), np.matmul(K, N)))

            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elmat

    def _get_mass_matrix(self, params, globdat):
        const_rho = type(self._rho) is float
        const_A = type(self._A) is float

        if const_rho and const_A:
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

                if not const_rho or not const_A:
                    x = self._shape.get_global_point(self._shape.get_integration_points()[:,ip], coords)[0]

                rho_ = self._rho if const_rho else eval(self._rho)
                A_ = self._A if const_A else eval(self._A)

                M = np.array([[rho_ * A_]])

                elmat += weights[ip] * np.matmul(np.transpose(N), np.matmul(M, N))

            params[pn.MATRIX2][np.ix_(idofs, idofs)] += elmat

    def _get_body_force(self, params, globdat):
        Q = np.array([self._q])
        for elem in self._elems:
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:self._rank, :]
            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords)

            elfor = np.zeros(self._dofcount)
            for ip in range(self._ipcount):
                N = np.zeros((1, self._dofcount * self._rank))
                N[0, :] = sfuncs[:, ip].transpose()
                elfor += weights[ip] * np.matmul(np.transpose(N), Q)

            params[pn.EXTFORCE][idofs] += elfor


def declare(factory):
    factory.declare_model('Bar', BarModel)
