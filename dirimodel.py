import numpy as np

from names import Actions as act
from names import ParamNames as pn
from names import GlobNames as gn
from model import *

import proputils as pu

GROUPS = 'groups'
DOFS = 'dofs'
VALS = 'values'
INCR = 'dispIncr'


class DirichletModel(Model):
    def take_action(self, action, params, globdat):
        if action == act.GETCONSTRAINTS:
            self._get_constraints(params, globdat)
        if action == act.ADVANCE:
            self._advance_step_constraints(params, globdat)

    def configure(self, props, globdat):
        self._groups = pu.parse_list(props[GROUPS])
        self._dofs = pu.parse_list(props[DOFS])
        self._vals = pu.parse_list(props[VALS], float)
        self._initDisp = self._vals
        if INCR in props:
            self._dispIncr = pu.parse_list(props[INCR], float)
        else:
            self._dispIncr = np.zeros(len(self._vals))


    def _get_constraints(self, params, globdat):
        ds = globdat[gn.DOFSPACE]
        for group, dof, val in zip(self._groups, self._dofs, self._vals):
            for node in globdat[gn.NGROUPS][group]:
                idof = ds.get_dof(node, dof)
                params[pn.CONSTRAINTS].add_constraint(idof, val)

    def _advance_step_constraints(self, params, globdat):
        self._vals = np.array(self._initDisp) + globdat[gn.TIMESTEP] * np.array(self._dispIncr)


def declare(factory):
    factory.declare_model('Dirichlet', DirichletModel)
