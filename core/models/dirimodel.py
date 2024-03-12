import numpy as np

from jive.names import GlobNames as gn
from jive.model.model import Model

import jive.util.proputils as pu

GROUPS = "groups"
DOFS = "dofs"
VALS = "values"
INCR = "dispIncr"

__all__ = ["DirichletModel"]


class DirichletModel(Model):
    def GETCONSTRAINTS(self, c, globdat, **kwargs):
        c = self._get_constraints(c, globdat, **kwargs)
        return c

    def ADVANCE(self, globdat):
        self._advance_step_constraints(globdat)

    def configure(self, props, globdat):
        self._groups = pu.parse_list(props[GROUPS])
        self._dofs = pu.parse_list(props[DOFS])
        self._vals = pu.parse_list(props[VALS], float)
        self._initDisp = self._vals
        if INCR in props:
            self._dispIncr = pu.parse_list(props[INCR], float)
        else:
            self._dispIncr = np.zeros(len(self._vals))

    def _get_constraints(self, c, globdat):
        ds = globdat[gn.DOFSPACE]
        for group, dof, val in zip(self._groups, self._dofs, self._vals):
            for node in globdat[gn.NGROUPS][group]:
                idof = ds.get_dof(node, dof)
                c.add_dirichlet(idof, val)
        return c

    def _advance_step_constraints(self, globdat):
        self._vals = np.array(self._initDisp) + globdat[gn.TIMESTEP] * np.array(
            self._dispIncr
        )
