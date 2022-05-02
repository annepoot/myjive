import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn
from names import PropNames as prn
from model import *
from node import Node

ELEMENTS = 'elements'
SUBTYPE = 'subtype'
LINEAR = 'linear'
NONLIN = 'nonlin'
EA = 'EA'
GAs = 'GAs'
EI = 'EI'
PLASTIC = 'plastic'
MP = 'Mp'
SHAPE = 'shape'
INTSCHEME = 'intScheme'
DOFTYPES = ['dx', 'dy', 'phi']


class FrameModel(Model):
    def take_action(self, action, params, globdat):

        if action == act.GETMATRIX0:
            self._get_matrix(params, globdat)

        elif action == act.GETMATRIXLB:
            self._get_matrix_lb(params, globdat)

        elif action == act.GETINTFORCE:
            self._get_int_force(params, globdat)

        elif action == act.CHECKCOMMIT:
            self._check_commit(params, globdat)

    def configure(self, props, globdat):
        self._subtype = str(props[SUBTYPE])
        self._EA = float(props[EA])
        self._GAs = float(props[GAs])
        self._EI = float(props[EI])

        self._plastic = bool(eval(props.get(PLASTIC, 'False')))
        if self._plastic:
            self._mp = float(props.get(MP, None))
        self._nhinges = 0
        self._hingedofs = []
        self._hingemoments = []

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
                    raise NotImplementedError('nonlinear strain only implemented for 2-node element')

                # TODO: use shape functions for evaluating psi and l_ locally
                ue = [globdat[gn.STATE0][i] for i in idofs]

                d = d0 + ue[3:5] - ue[0:2]
                lcps = (d0[0] * d[0] + d0[1] * d[1]) / l_0
                lsps = (d0[0] * d[1] - d0[1] * d[0]) / l_0

                for ip in range(self._ipcount):
                    N = sfuncs[:, ip]
                    dN = grads[:, 0, ip]

                    theta = np.matmul(N, ue[2::3])
                    kappa = np.matmul(dN, ue[2::3])
                    omega = phi + theta
                    gamma = (np.cos(theta) * lsps - np.sin(theta) * lcps) / l_0
                    eps = (np.cos(theta) * lcps + np.sin(theta) * lsps) / l_0 - 1
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

            ue = [globdat[gn.STATE0][i] for i in idofs]

            if self._subtype == LINEAR:
                for ip in range(self._ipcount):
                    N = sfuncs[:, ip]
                    dN = grads[:, 0, ip]

                    B = self._get_B_matrix(N=N, dN=dN, omega=phi)
                    elmat = weights[ip] * np.matmul(B.transpose(), np.matmul(D, B))
                    elfor += np.matmul(elmat, ue)

            if self._subtype == NONLIN:
                d = d0 + ue[3:5] - ue[0:2]
                lcps = (d0[0] * d[0] + d0[1] * d[1]) / l_0
                lsps = (d0[0] * d[1] - d0[1] * d[0]) / l_0
                for ip in range(self._ipcount):
                    N = sfuncs[:, ip]
                    dN = grads[:, 0, ip]

                    theta = np.matmul(N, [ue[2], ue[5]])
                    kappa = np.matmul(dN, [ue[2], ue[5]])
                    omega = phi + theta
                    gamma = (np.cos(theta) * lsps - np.sin(theta) * lcps) / l_0
                    eps = (np.cos(theta) * lcps + np.sin(theta) * lsps) / l_0 - 1
                    evec = [eps, gamma, kappa]
                    svec = np.matmul(D, evec)

                    B = self._get_B_matrix(N=N, dN=dN, omega=omega, gamma=gamma, eps=eps)
                    elfor += weights[ip] * np.matmul(B.transpose(), svec)

            params[pn.INTFORCE][np.ix_(idofs)] += elfor

        if self._nhinges > 0:
            params[pn.INTFORCE][np.ix_(self._hingedofs)] = self._hingemoments

    def _check_commit(self, params, globdat):
        globdat[gn.ACCEPTED] = True
        if self._plastic:
            if self._mp is None:
                raise RuntimeError('Plastic moment ' + MP + ' has not been defined')

            moments = self._get_moments(globdat)

            # Get element and node where highest moment occurs
            maxarg = np.argmax(abs(moments))
            maxratio = np.max(abs(moments)) / self._mp
            hingeelem = maxarg // 2
            enodes = globdat[gn.ESET][hingeelem].get_nodes()
            hingenode = enodes[maxarg % 2]
            if abs(np.min(moments)) < abs(np.max(moments)):
                sign = 1
            else:
                sign = -1

            if maxratio > 1.:
                globdat[gn.ACCEPTED] = False
                self._add_plastic_hinge(globdat, params, hingenode, hingeelem, sign)

    def _add_plastic_hinge(self, globdat, params, hingenode, hingeelem, sign):
        print('Adding plastic hinge on node %i (in element %i)' % (hingenode, hingeelem))

        # Add node
        oldnode = hingenode
        coords = globdat[gn.NSET][oldnode].get_coords()
        globdat[gn.NSET].append(Node(coords))
        newnode = len(globdat[gn.NSET]) - 1

        # Duplicate dofs and add new phi dof
        globdat[gn.DOFSPACE].set_dof(oldnode, newnode, 'dx')
        globdat[gn.DOFSPACE].set_dof(oldnode, newnode, 'dy')
        globdat[gn.DOFSPACE].add_dof(newnode, 'phi')

        # Initialize new dof
        oldphidof = globdat[gn.DOFSPACE].get_dof(oldnode, 'phi')
        newphidof = globdat[gn.DOFSPACE].get_dof(newnode, 'phi')
        phi_state0 = globdat[gn.OLDSTATE0][oldphidof]
        phi_oldstate0 = globdat[gn.BACKUPSTATE0][oldphidof]
        globdat[gn.STATE0] = np.append(globdat[gn.OLDSTATE0], phi_state0)
        globdat[gn.OLDSTATE0] = np.append(globdat[gn.BACKUPSTATE0], phi_oldstate0)

        # Update element connectivity
        globdat[gn.ESET][hingeelem].change_node(oldnode, newnode)
        self._elems = globdat[gn.ESET]

        # Modify hinge variables
        self._nhinges += 1
        self._hingemoments.append(-self._mp * sign)
        self._hingemoments.append(self._mp * sign)
        self._hingedofs.append(newphidof)
        self._hingedofs.append(oldphidof)

    def _get_moments(self, globdat):
        moments = np.zeros(2 * len(self._elems))
        D = self._get_D_matrix()
        for i, elem in enumerate(self._elems):
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES)
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=0)[:, :]

            d0 = coords[1, :] - coords[0, :]
            l_0 = np.linalg.norm(d0)
            coords1d = np.array([0, l_0])

            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords1d)
            m1 = 0
            m2 = 0

            ue = [globdat[gn.STATE0][i] for i in idofs]

            d = d0 + ue[3:5] - ue[0:2]
            lcps = (d0[0] * d[0] + d0[1] * d[1]) / l_0
            lsps = (d0[0] * d[1] - d0[1] * d[0]) / l_0

            for ip in range(self._ipcount):
                N = sfuncs[:, ip]
                dN = grads[:, 0, ip]

                theta = np.matmul(N, [ue[2], ue[5]])
                kappa = np.matmul(dN, [ue[2], ue[5]])
                gamma = (np.cos(theta) * lsps - np.sin(theta) * lcps) / l_0
                eps = (np.cos(theta) * lcps + np.sin(theta) * lsps) / l_0 - 1
                evec = [eps, gamma, kappa]
                svec = np.matmul(D, evec)
                g = gamma / 2
                e = (1 + eps) / 2
                t = 1 / l_0

                m1 += weights[ip] * np.matmul([g, -e, -t], svec)
                m2 += weights[ip] * np.matmul([g, -e, t], svec)

            moments[2 * i] = m1
            moments[2 * i + 1] = m2
        return moments

    def _get_B_matrix(self, N, dN, omega, gamma=0, eps=0):
        B = np.zeros((self._strcount, self._dofcount))
        for inode in range(self._shape.node_count()):
            i = 3 * inode
            c = np.cos(omega) * dN[inode]
            s = np.sin(omega) * dN[inode]
            B[:, i:(i + 3)] = np.array([[c, s, N[inode] * gamma],
                                        [-s, c, -N[inode] * (1 + eps)],
                                        [0, 0, dN[inode]]])
        return B

    def _get_WN_matrix(self, N, dN, omega, eps=0):
        WN = np.zeros((6, 6))
        dn = dN[1]
        l0 = 1 / dn
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
        l0 = 1 / dn
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
