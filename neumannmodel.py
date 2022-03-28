import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn
from model import *

import proputils as pu

GROUPS = 'groups'
DOFS = 'dofs'
VALS = 'values'
INCR = 'loadIncr'


class NeumannModel(Model):
    def take_action(self, action, params, globdat):
        if action == act.GETEXTFORCE:
            self._get_ext_force(params, globdat)
        if action == act.ADVANCEFORCE:
            self._advance_step_force(params, globdat)

    def configure(self, props, globdat):
        self._groups = pu.parse_list(props[GROUPS])
        self._dofs = pu.parse_list(props[DOFS])
        self._vals = pu.parse_list(props[VALS], float)
        self._initLoad = self._vals
        if INCR in props:
            self._loadIncr = pu.parse_list(props[INCR], float)

    def _get_ext_force(self, params, globdat):
        ds = globdat[gn.DOFSPACE]
        for group, dof, val in zip(self._groups, self._dofs, self._vals):
            for node in globdat[gn.NGROUPS][group]:
                idof = ds.get_dof(node, dof)
                params[pn.EXTFORCE][idof] += val

    def _advance_step_force(self, params, globdat):
        self._vals = np.array(self._initLoad) + np.array(globdat[gn.TIMESTEP] * self._loadIncr)


def declare(factory):
    factory.declare_model('Neumann', NeumannModel)
