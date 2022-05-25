import numpy as np

import proputils as pu
from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn
from barmodel import BarModel

LOAD = 'q'
DOFTYPES = ['dx']


class XBarModel(BarModel):
    def take_action(self, action, params, globdat):
        print('XBarModel taking action', action)

        # Refer the core actions to the parent class
        super().take_action(action, params, globdat)

        # Add extended actions below
        if action == act.GETEXTFORCE:
            self._get_body_force(params, globdat)

    def configure(self, props, globdat):

        # Perform the configuration of the parent class
        super().configure(props, globdat)

        # Get the body force on the bar
        self._q = float(props.get(LOAD, 0))

        # Try to turn all properties into floats
        self._E = pu.try_float(self._E)
        self._A = pu.try_float(self._A)
        self._k = pu.try_float(self._k)
        self._rho = pu.try_float(self._rho)
        self._q = pu.try_float(self._q)

    def _get_matrix(self, params, globdat):
        E_ = params.get(pn.YOUNG, self._E)
        A_ = params.get(pn.AREA, self._A)
        k_ = params.get(pn.SPRING, self._k)
        eval_params = {}

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

                if type(E_) is str or type(A_) is str or type(k_) is str:
                    x = self._shape.get_global_point(self._shape.get_integration_points()[:,ip], coords)[0]
                    eval_params['x'] = x

                E = pu.evaluate(E_, eval_params)
                A = pu.evaluate(A_, eval_params)
                k = pu.evaluate(k_, eval_params)

                D = np.array([[E * A]])
                K = np.array([[k]])

                elmat += weights[ip] * (np.matmul(np.transpose(B), np.matmul(D, B))
                                        + np.matmul(np.transpose(N), np.matmul(K, N)))

            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elmat

    def _get_mass_matrix(self, params, globdat):
        rho_ = params.get(pn.RHO, self._rho)
        A_ = params.get(pn.AREA, self._A)
        eval_params = {}

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

                if type(rho_) is str or type(A_) is str:
                    x = self._shape.get_global_point(self._shape.get_integration_points()[:,ip], coords)[0]
                    eval_params['x'] = x

                rho = pu.evaluate(rho_, eval_params)
                A = pu.evaluate(A_, eval_params)

                M = np.array([[rho * A]])

                elmat += weights[ip] * np.matmul(np.transpose(N), np.matmul(M, N))

            params[pn.MATRIX2][np.ix_(idofs, idofs)] += elmat

    def _get_body_force(self, params, globdat):
        q_ = params.get(pn.LOAD, self._q)
        eval_params = {}

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

                if type(q_) is str:
                    x = self._shape.get_global_point(self._shape.get_integration_points()[:,ip], coords)[0]
                    eval_params['x'] = x

                q = pu.evaluate(q_, eval_params)

                Q = np.array([q])

                elfor += weights[ip] * np.matmul(np.transpose(N), Q)

            params[pn.EXTFORCE][idofs] += elfor


def declare(factory):
    factory.declare_model('XBar', XBarModel)
