import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn
from names import PropNames as prn
from model import *

ELEMENTS = 'elements'
SUBTYPE = 'subtype'
LINEAR = 'linear'
NONLIN = 'nonlin'
EA = 'EA'
GAs = 'GAs'
EI = 'EI'
SHAPE = 'shape'
INTSCHEME = 'intScheme'
DOFTYPES = ['dx', 'dy', 'phi']


class FrameModel(Model):
    def take_action(self, action, params, globdat):
        print('Model taking action', action)

        if action == act.GETMATRIX0:
            self._get_matrix(params, globdat)

        elif action == act.GETMATRIXLB:
            self._get_matrix_lb(params, globdat)
            
        elif action == act.GETINTFORCE:
            self._get_int_force(params, globdat)

    def configure(self, props, globdat):
        self._subtype = str(props[SUBTYPE])
        self._EA = float(props[EA])
        self._GAs = float(props[GAs])
        self._EI = float(props[EI])
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][prn.TYPE], props[SHAPE][INTSCHEME])
        egroup = globdat[gn.EGROUPS][props[ELEMENTS]]
        self._elems = [globdat[gn.ESET][e] for e in egroup]

        self._ipcount = self._shape.ipoint_count()
        self._dofcount = 3 * self._shape.node_count()
        self._strcount = 3

        nodes = np.unique([node for elem in self._elems for node in elem.get_nodes()])
        for doftype in DOFTYPES:
            globdat[gn.DOFSPACE].add_type(doftype)
            for node in nodes:
                globdat[gn.DOFSPACE].add_dof(node, doftype)

    def _get_matrix(self, params, globdat):
        D = self._get_D_matrix()
        for elem in self._elems:
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES)
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=0)[:, :]

            d0 = coords[1, :] - coords[0, :]
            phi = np.arctan2(d0[1], d0[0])
            l_0 = np.linalg.norm(d0)
            coords1d = np.array([0, l_0])

            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords1d)
            elmat = np.zeros((6, 6))

            if self._subtype == LINEAR:
                for ip in range(self._ipcount):
                    N = sfuncs[:, ip]
                    dN = grads[:, 0, ip]

                    B = self._get_B_matrix(N=N, dN=dN, omega=phi)
                    elmat += weights[ip] * np.matmul(B.transpose(), np.matmul(D, B))

            elif self._subtype == NONLIN:

                if self._shape.node_count() > 2:
                    raise ImplementationError('nonlinear strain only implemented for 2-node element')

                # TODO: use shape functions for evaluating psi and l_ locally
                ue = [globdat[gn.STATE0][i] for i in idofs]
                d = d0 + ue[3:4] - ue[0:1]
                l_ = np.linalg.norm(d)
                psi = np.arccos((d0[0] * d[0] + d0[1] * d[1]) / l_0 / l_)

                  
                for ip in range(self._ipcount):
                    N = sfuncs[:, ip]
                    dN = grads[:, 0, ip]

                    theta = np.matmul(N, ue[2::3])
                    kappa = np.matmul(dN, ue[2::3])
                    omega = phi + theta
                    gamma = (l_ * (np.cos(theta) * np.sin(psi) + np.sin(theta) * np.cos(psi))) / l_0
                    eps = (l_ * (np.cos(theta) * np.cos(psi) + np.sin(theta) * np.sin(psi))) / l_0 - 1
                    evec = [eps, gamma, kappa]
                    svec = np.matmul(D, evec)

                    B = self._get_B_matrix(N=N, dN=dN, omega=omega, gamma=gamma, eps=eps)
                    WN = self._get_WN_matrix(N=N, dN=dN, omega=omega, eps=eps)
                    WV = self._get_WV_matrix(N=N, dN=dN, omega=omega, gamma=gamma)
                    Kmat = weights[ip] * np.matmul(B.transpose(), np.matmul(D, B))
                    Kgeo = svec[0] * WN + svec[1] * WV  # *l0/l0
                    elmat += Kmat + Kgeo

            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elmat
            
    def _get_matrix_lb(self, params, globdat):
        D = self._get_D_matrix()
        for elem in self._elems:
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES)
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=0)[:, :]

            d0 = coords[1, :] - coords[0, :]
            phi = np.arctan2(d0[1], d0[0])
            l_0 = np.linalg.norm(d0)
            coords1d = np.array([0, l_0])

            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords1d)
            elmatM = np.zeros((6, 6))
            elmatG = np.zeros((6, 6))

            for ip in range(self._ipcount):
                N = sfuncs[:, ip]
                dN = grads[:, 0, ip]

                B = self._get_B_matrix(N=N, dN=dN, omega=phi)

                ue = [globdat[gn.STATE0][i] for i in idofs]
                evec = np.matmul(B, ue)
                svec = np.matmul(D, evec)

                WN = self._get_WN_matrix(N=N, dN=dN, omega=phi, eps=0)
                WV = self._get_WV_matrix(N=N, dN=dN, omega=phi, gamma=0)

                elmatM += weights[ip] * np.matmul(B.transpose(), np.matmul(D, B))
                elmatG += svec[0] * WN + svec[1] * WV  # *l0/l0

            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elmatM
            params[pn.MATRIX1][np.ix_(idofs, idofs)] += elmatG

    def _get_int_force(self, params, globdat):
        D = self._get_D_matrix()
        for elem in self._elems:
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES)
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=0)[:, :]

            d0 = coords[1, :] - coords[0, :]
            phi = np.arctan2(d0[1], d0[0])
            l_0 = np.linalg.norm(d0)
            coords1d = np.array([0, l_0])

            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords1d)
            elfor = np.zeros(6)

            if self._subtype == LINEAR:
                ue = [globdat[gn.STATE0][i] for i in idofs]
                for ip in range(self._ipcount):
                    N = sfuncs[:, ip]
                    dN = grads[:, 0, ip]

                    B = self._get_B_matrix(N=N, dN=dN, omega=phi)
                    elmat = weights[ip] * np.matmul(B.transpose(), np.matmul(D, B))
                    elfor += np.matmul(elmat, ue)

            if self._subtype == NONLIN:
                ue = [globdat[gn.STATE0][i] for i in idofs]
                d = d0 + ue[3:4] - ue[0:1]
                l_ = np.linalg.norm(d)
                psi = np.arccos((d0[0] * d[0] + d0[1] * d[1]) / l_0 / l_)
                for ip in range(self._ipcount):
                    N = sfuncs[:, ip]
                    dN = grads[:, 0, ip]

                    theta = np.matmul(N, [ue[2], ue[5]])
                    kappa = np.matmul(dN, [ue[2], ue[5]])
                    omega = phi + theta
                    gamma = (l_ * (np.cos(theta) * np.sin(psi) + np.sin(theta) * np.cos(psi))) / l_0
                    eps = (l_ * (np.cos(theta) * np.cos(psi) + np.sin(theta) * np.sin(psi))) / l_0 - 1
                    evec = [eps, gamma, kappa]
                    svec = np.matmul(D, evec)

                    B = self._get_B_matrix(N=N, dN=dN, omega=omega, gamma=gamma, eps=eps)
                    elfor += weights[ip] * np.matmul(B.transpose(), svec)

            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elfor

    def _get_B_matrix(self, N, dN, omega, gamma=0, eps=0):

        B = np.zeros((self._strcount,self._dofcount))
        for inode in range(self._shape.node_count()):
            i = 3 * inode
            c = np.cos(omega) * dN[inode]
            s = np.sin(omega) * dN[inode]
            B[:,i:(i+3)] = np.array([[ c, s, N[inode]*gamma],
                                     [-s, c, -N[inode]*(1+eps)],
                                     [ 0, 0, dN[inode]]])
        return B

    def _get_WN_matrix(self, N, dN, omega, eps=0):
        WN = np.zeros((6, 6))

        dn = dN[1]
        l0 = 1/dn
        c = np.cos(omega) 
        s = np.sin(omega)
        n1 = N[0]
        n2 = N[1]
        WN[0, 2] = n1 * s
        WN[0, 5] = n2 * s
        WN[1, 2] = -n1 * c
        WN[1, 5] = -n2 * c
        WN[2, 2] = -n1 ** 2 * l0 * (1 + eps)
        WN[2, 3] = -n1 * s
        WN[2, 4] = n1 * c
        WN[2, 5] = -n1 * n2 * l0 * (1 + eps)
        WN[3, 5] = -n2 * s
        WN[4, 5] = n2 * c
        WN[5, 5] = -n2 ** 2 * l0 * (1 + eps)
        WN = WN + WN.transpose() - np.diag(np.diag(WN))  # Add symmetric values in lower triangle
        return WN

    def _get_WV_matrix(self, N, dN, omega, gamma=0):
        WV = np.zeros((6, 6))

        dn = dN[1]
        l0 = 1/dn
        c = np.cos(omega)
        s = np.sin(omega)
        n1 = N[0]
        n2 = N[1]
        WV[0, 2] = n1 * c
        WV[0, 5] = n2 * c
        WV[1, 2] = n1 * s
        WV[1, 5] = n2 * s
        WV[2, 2] = -n1 ** 2 * l0 * gamma
        WV[2, 3] = -n1 * c
        WV[2, 4] = n1 * s
        WV[2, 5] = -n1 * n2 * l0 * gamma
        WV[3, 5] = -n2 * c
        WV[4, 5] = n2 * s
        WV[5, 5] = -n2 ** 2 * l0 * gamma
        WV = WV + WV.transpose() - np.diag(np.diag(WV))  # Add symmetric values in lower triangle
        return WV

    def _get_D_matrix(self):
        D = np.diag([self._EA, self._GAs, self._EI])
        return D


def declare(factory):
    factory.declare_model('Frame', FrameModel)
