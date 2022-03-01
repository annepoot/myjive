import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn
from model import *

import proputils as pu

GROUPS = 'groups'
DOFS = 'dofs'
VALS = 'values'


class NeumannModel(Model):
    def take_action(self, action, params, globdat):
        if action == act.GETEXTFORCE:
            self._get_ext_force(params, globdat)

    def configure(self, props, globdat):
        self._groups = pu.parse_list(props[GROUPS])
        self._dofs = pu.parse_list(props[DOFS])
        self._vals = pu.parse_list(props[VALS],float)

    def _get_ext_force(self, params, globdat):
        ds = globdat[gn.DOFSPACE]
        for group, dof, val in zip(self._groups, self._dofs, self._vals):
            for node in globdat[gn.NGROUPS][group]:
                idof = ds.get_dof(node, dof)
                params[pn.EXTFORCE][idof] += val


def declare(factory):
    factory.declare_model('Neumann', NeumannModel)
