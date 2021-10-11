import numpy as np

from names import Actions    as act
from names import ParamNames as pn
from names import GlobNames  as gn
from names import PropNames  as prn
from model import *

ELEMENTS  = 'elements'
YOUNG     = 'young'
POISSON   = 'poisson'
SHAPE     = 'shape'
INTSCHEME = 'intScheme'
STATE     = 'state'
DOFTYPES  = ['dx','dy','dz']
PE_STATE  = 'plane_strain'
PS_STATE  = 'plane_stress'

class ElasticModel (Model):
    def take_action (self, action, params, globdat):
        print('ElasticModel taking action',action)
        if action == act.GETMATRIX0:
            self._get_matrix(params, globdat)

    def configure (self, props, globdat):
        self._young   = float(props[YOUNG])
        self._poisson = float(props[POISSON])
        self._shape   = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][prn.TYPE], props[SHAPE][INTSCHEME])
        egroup        = globdat[gn.EGROUPS][props[ELEMENTS]]
        self._elems   = [globdat[gn.ESET][e] for e in egroup]

        if self._shape.global_rank() != globdat[gn.MESHRANK]:
            raise RuntimeError ('ElasticModel: Shape rank must agree with mesh rank')

        self._rank     = self._shape.global_rank()
        self._ipcount  = self._shape.ipoint_count() 
        self._dofcount = self._rank * self._shape.node_count()
        self._strcount = self._rank * ( self._rank + 1 ) // 2;  # 1-->1, 2-->3, 3-->6

        if self._rank == 2:
            self._state = props[STATE]
            if not self._state in (PE_STATE,PS_STATE):
                raise RuntimeError ('ElasticModel: state in 2d should be plane_strain or plane_stress')
        else:
            self._state = ''

        nodes = np.unique ([node for elem in self._elems for node in elem.get_nodes()])
        for doftype in DOFTYPES[0:self._rank]:
            globdat[gn.DOFSPACE].add_type (doftype)
            for node in nodes:
                globdat[gn.DOFSPACE].add_dof(node,doftype)

    def _get_matrix (self, params, globdat):
        D = self._get_D_matrix()
        for elem in self._elems:
            inodes = elem.get_nodes()
            idofs  = globdat[gn.DOFSPACE].get_dofs(inodes,DOFTYPES[0:self._rank])
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes],axis=1)[0:self._rank,:]
            grads, weights = self._shape.get_shape_gradients(coords)

            elmat = np.zeros((self._dofcount,self._dofcount))
            for ip in range(self._ipcount):
                B = self._get_B_matrix ( grads[:,:,ip] )
                elmat += weights[ip] * np.matmul(np.transpose(B),np.matmul(D,B))
            
            params[pn.MATRIX0][np.ix_(idofs,idofs)] += elmat

    def _get_B_matrix ( self, grads ):
        B = np.zeros((self._strcount,self._dofcount))
        if self._rank == 1:
            B = grads
        elif self._rank == 2:
            for inode in range(self._shape.node_count()):
                i = 2*inode
                gi = grads[inode,:]
                B[0:3,i:(i+2)] = [
                    [gi[0],   0.],
                    [   0., gi[1]],
                    [gi[1], gi[0]]
                ]
        elif self._rank == 3:
            B = np.zeros((6,self._dofcount))
            for inode in range(self._shape.node_count()):
                i = 3*inode
                gi = grads[inode,:]
                B[0:6,i:(i+3)] = [
                    [gi[0],    0.,    0.],
                    [   0., gi[1],    0.],
                    [   0.,    0., gi[2]],
                    [gi[1], gi[0],    0.],
                    [   0., gi[2], gi[1]],
                    [gi[2],    0., gi[0]]
                ]
        return B

    def _get_D_matrix ( self ):
        D = np.zeros((self._strcount,self._strcount))
        E = self._young
        nu = self._poisson
        g = 0.5*E / (1.+nu)
        if self._rank == 1:
            D = self._young
        elif self._rank == 3:
            a = E*(1.-nu) / ((1.+nu)*(1.-2.*nu))
            b = E*nu / ((1.+nu)*(1.-2.*nu))
            D = [
                [ a, b, b,.0,.0,.0],
                [ b, a, b,.0,.0,.0],
                [ b, b, a,.0,.0,.0],
                [.0,.0,.0, g,.0,.0],
                [.0,.0,.0,.0, g,.0],
                [.0,.0,.0,.0,.0, g]
            ]
        elif self._rank == 2:
            if self._state == PE_STATE:
                a = E*(1.-nu) / ((1.+nu)*(1.-2.*nu))
                b = E*nu / ((1.+nu)*(1.-2.*nu))
                D = [
                    [ a, b,.0],
                    [ b, a,.0],
                    [.0,.0, g]
                ]
            else: 
                assert(self._state == PS_STATE)
                a = E/(1.-nu*nu)
                b = a*nu
                D = [
                    [ a, b,.0],
                    [ b, a,.0],
                    [.0,.0, g]
                ]
        return D

def declare (factory):
    factory.declare_model ('Elastic', ElasticModel)

