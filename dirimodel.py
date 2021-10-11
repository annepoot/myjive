import numpy as np

from names import Actions    as act
from names import ParamNames as pn
from names import GlobNames  as gn
from model import *

GROUPS = 'groups'
DOFS   = 'dofs'
VALS   = 'values'

class DirichletModel (Model):
    def take_action (self, action, params, globdat):
        if action == act.GETCONSTRAINTS:
            self._get_constraints(params, globdat)

    def configure (self, props, globdat):
        self._groups = props[GROUPS].strip('[').strip(']').split(',')
        self._dofs   = props[DOFS].strip('[').strip(']').split(',')
        self._vals   = props[VALS].strip('[').strip(']').split(',')

    def _get_constraints (self, params, globdat):
        ds = globdat[gn.DOFSPACE]
        for group,dof,val in zip(self._groups,self._dofs,self._vals):
            for node in globdat[gn.NGROUPS][group]: 
                idof = ds.get_dof (node,dof)    
                params[pn.CONSTRAINTS].add_constraint(idof,float(val))

def declare (factory):
    factory.declare_model ('Dirichlet', DirichletModel)

