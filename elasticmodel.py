import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn
from names import PropNames as prn
from model import Model

ELEMENTS = 'elements'
YOUNG = 'young'
RHO = 'rho'
THICKNESS = 'thickness'
POISSON = 'poisson'
SHAPE = 'shape'
INTSCHEME = 'intScheme'
STATE = 'state'
DOFTYPES = ['dx', 'dy', 'dz']
PE_STATE = 'plane_strain'
PS_STATE = 'plane_stress'


class ElasticModel(Model):
    def take_action(self, action, params, globdat):
        print('ElasticModel taking action', action)
        if action == act.GETMATRIX0:
            self._get_matrix(params, globdat)
        elif action == act.GETMATRIX2:
            self._get_mass_matrix(params, globdat)
        elif action == act.GETEXTFORCE:
            self._get_body_force(params, globdat)
        elif action == act.GETTABLE:
            if 'stress' in params[pn.TABLENAME]:
                self._get_stresses(params, globdat)

    def configure(self, props, globdat):
        self._young = float(props[YOUNG])
        self._rho = float(props.get(RHO,0))
        if globdat[gn.MESHRANK] == 2:
            self._thickness = float(props.get(THICKNESS,1))
        if globdat[gn.MESHRANK] > 1:
            self._poisson = float(props[POISSON])
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][prn.TYPE], props[SHAPE][INTSCHEME])
        egroup = globdat[gn.EGROUPS][props[ELEMENTS]]
        self._elems = [globdat[gn.ESET][e] for e in egroup]

        if self._shape.global_rank() != globdat[gn.MESHRANK]:
            raise RuntimeError('ElasticModel: Shape rank must agree with mesh rank')

        self._rank = self._shape.global_rank()
        self._ipcount = self._shape.ipoint_count()
        self._dofcount = self._rank * self._shape.node_count()
        self._strcount = self._rank * (self._rank + 1) // 2;  # 1-->1, 2-->3, 3-->6

        if self._rank == 2:
            self._state = props[STATE]
            if not self._state in (PE_STATE, PS_STATE):
                raise RuntimeError('ElasticModel: state in 2d should be plane_strain or plane_stress')
        else:
            self._state = ''

        nodes = np.unique([node for elem in self._elems for node in elem.get_nodes()])
        for doftype in DOFTYPES[0:self._rank]:
            globdat[gn.DOFSPACE].add_type(doftype)
            for node in nodes:
                globdat[gn.DOFSPACE].add_dof(node, doftype)

    def _get_matrix(self, params, globdat):
        D = self._get_D_matrix()
        for elem in self._elems:
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:self._rank, :]
            grads, weights = self._shape.get_shape_gradients(coords)

            elmat = np.zeros((self._dofcount, self._dofcount))
            for ip in range(self._ipcount):
                B = self._get_B_matrix(grads[:, :, ip])
                elmat += weights[ip] * np.matmul(np.transpose(B), np.matmul(D, B)) * self._thickness

            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elmat

    def _get_mass_matrix(self, params, globdat):
        M = params.get(pn.RHO, self._rho) * np.identity(self._rank)
        for elem in self._elems:
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:self._rank, :]
            sfuncs = self._shape.get_shape_functions()
            weights = self._shape.get_integration_weights(coords)

            elmat = np.zeros((self._dofcount, self._dofcount))
            for ip in range(self._ipcount):
                N = self._get_N_matrix(sfuncs[:, ip])
                elmat += weights[ip] * np.matmul(np.transpose(N), np.matmul(M, N))

            params[pn.MATRIX2][np.ix_(idofs, idofs)] += elmat

    def _get_body_force(self, params, globdat):
        if self._rank == 2:
            gravity = np.array([0, -1])
            b = self._rho * gravity
            for elem in self._elems:
                inodes = elem.get_nodes()
                idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
                coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:self._rank, :]
                sfuncs = self._shape.get_shape_functions()
                weights = self._shape.get_integration_weights(coords)

                elfor = np.zeros(self._dofcount)
                for ip in range(self._ipcount):
                    N = self._get_N_matrix(sfuncs[:, ip])
                    elfor += weights[ip] * self._thickness * np.matmul(np.transpose(N), b)

                params[pn.EXTFORCE][idofs] += elfor
        else:
            pass

    def _get_stresses(self, params, globdat):
        D = self._get_D_matrix()
        table = params[pn.TABLE]
        tbwts = params[pn.TABLEWEIGHTS]

        if 'xx' not in table:
            table['xx'] = np.zeros(len(globdat[gn.NSET]))
        if self._rank > 1:
            if 'yy' not in table:
                table['yy'] = np.zeros(len(globdat[gn.NSET]))
            if 'xy' not in table:
                table['xy'] = np.zeros(len(globdat[gn.NSET]))
        elif self._rank > 2:
            if 'zz' not in table:
                table['zz'] = np.zeros(len(globdat[gn.NSET]))
            if 'yz' not in table:
                table['yz'] = np.zeros(len(globdat[gn.NSET]))
            if 'zx' not in table:
                table['zx'] = np.zeros(len(globdat[gn.NSET]))

        for elem in self._elems:
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:self._rank, :]
            grads, weights = self._shape.get_shape_gradients(coords)

            eldisp = globdat[gn.STATE0][idofs]
            sfuncs = self._shape.get_shape_functions()

            elsig = np.zeros((self._shape.node_count(), self._strcount))
            elwts = np.zeros(self._shape.node_count())

            for ip in range(self._ipcount):
                B = self._get_B_matrix(grads[:, :, ip])
                ptsig = np.matmul(D, np.matmul(B, eldisp))
                elsig += np.outer(sfuncs[:, ip], ptsig)
                elwts += sfuncs[:, ip].flatten()

            tbwts[inodes] += elwts

            table['xx'][inodes] += elsig[:, 0]

            if self._rank == 2:
                table['yy'][inodes] += elsig[:, 1]
                table['xy'][inodes] += elsig[:, 2]
            elif self._rank == 3:
                table['yy'][inodes] += elsig[:, 1]
                table['zz'][inodes] += elsig[:, 2]
                table['xy'][inodes] += elsig[:, 3]
                table['yz'][inodes] += elsig[:, 4]
                table['zx'][inodes] += elsig[:, 5]

    def _get_N_matrix(self, sfuncs):
        N = np.zeros((self._rank, self._dofcount))
        for i in range(self._rank):
            N[i,i::self._rank] = sfuncs.transpose()
        return N

    def _get_B_matrix(self, grads):
        B = np.zeros((self._strcount, self._dofcount))
        if self._rank == 1:
            B = grads.transpose()
        elif self._rank == 2:
            for inode in range(self._shape.node_count()):
                i = 2 * inode
                gi = grads[inode, :]
                B[0:3, i:(i + 2)] = [
                    [gi[0], 0.],
                    [0., gi[1]],
                    [gi[1], gi[0]]
                ]
        elif self._rank == 3:
            B = np.zeros((6, self._dofcount))
            for inode in range(self._shape.node_count()):
                i = 3 * inode
                gi = grads[inode, :]
                B[0:6, i:(i + 3)] = [
                    [gi[0], 0., 0.],
                    [0., gi[1], 0.],
                    [0., 0., gi[2]],
                    [gi[1], gi[0], 0.],
                    [0., gi[2], gi[1]],
                    [gi[2], 0., gi[0]]
                ]
        return B

    def _get_D_matrix(self):
        D = np.zeros((self._strcount, self._strcount))
        E = self._young
        if self._rank == 1:
            D[[0]] = E
            return D
        nu = self._poisson
        g = 0.5 * E / (1. + nu)
        if self._rank == 3:
            a = E * (1. - nu) / ((1. + nu) * (1. - 2. * nu))
            b = E * nu / ((1. + nu) * (1. - 2. * nu))
            D = [
                [a, b, b, .0, .0, .0],
                [b, a, b, .0, .0, .0],
                [b, b, a, .0, .0, .0],
                [.0, .0, .0, g, .0, .0],
                [.0, .0, .0, .0, g, .0],
                [.0, .0, .0, .0, .0, g]
            ]
        elif self._rank == 2:
            if self._state == PE_STATE:
                a = E * (1. - nu) / ((1. + nu) * (1. - 2. * nu))
                b = E * nu / ((1. + nu) * (1. - 2. * nu))
                D = [
                    [a, b, .0],
                    [b, a, .0],
                    [.0, .0, g]
                ]
            else:
                assert (self._state == PS_STATE)
                a = E / (1. - nu * nu)
                b = a * nu
                D = [
                    [a, b, .0],
                    [b, a, .0],
                    [.0, .0, g]
                ]
        return D


def declare(factory):
    factory.declare_model('Elastic', ElasticModel)
