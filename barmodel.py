import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames  as gn
from names import PropNames  as prn
from model import *

ELEMENTS  = 'elements'
YOUNG     = 'young'
SHAPE     = 'shape'
INTSCHEME = 'intScheme'
DOFTYPES  = ['dx']


class BarModel (Model):
    def take_action (self, action, params, globdat):
        print('BarModel taking action', action)

        if action == act.GETMATRIX0:
            self.__stiffness(params, globdat)

    def configure (self, props, globdat):
        self._young = float(props[YOUNG])
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][prn.TYPE], props[SHAPE][INTSCHEME])
        egroup = globdat[gn.EGROUPS][props[ELEMENTS]]
        self._elems = [globdat[gn.ESET][e] for e in egroup]

        self._rank     = self._shape.global_rank()
        self._ipcount  = self._shape.ipoint_count() 
        self._dofcount = self._rank * self._shape.node_count()
        self._strcount = self._rank * ( self._rank + 1 ) // 2

        nodes = np.unique ([node for elem in self._elems for node in elem.get_nodes()])
        for doftype in DOFTYPES[0:self._rank]:
            globdat[gn.DOFSPACE].add_type (doftype)
            for node in nodes:
                globdat[gn.DOFSPACE].add_dof(node,doftype)
        

    def __stiffness (self, params, globdat):
        D = np.array([[self._young]])
        for elem in self._elems:
            inodes = elem.get_nodes()
            idofs  = globdat[gn.DOFSPACE].get_dofs(inodes,DOFTYPES[0:self._rank])
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes],axis=1)[0:self._rank,:]
            grads, weights = self._shape.get_shape_gradients(coords)

            elmat = np.zeros((self._dofcount,self._dofcount))
            for ip in range(self._ipcount):
                B = grads[:,:,ip].transpose()
                elmat += weights[ip] * np.matmul(np.transpose(B),np.matmul(D,B))
            
            params[pn.MATRIX0][np.ix_(idofs,idofs)] += elmat
     


def declare (factory):
    factory.declare_model ('Bar', BarModel)

