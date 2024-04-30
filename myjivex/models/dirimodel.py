import numpy as np

from myjive.names import GlobNames as gn
from myjive.model.model import Model
from myjive.util.proputils import check_list

__all__ = ["DirichletModel"]


class DirichletModel(Model):
    def GETCONSTRAINTS(self, c, globdat, **kwargs):
        c = self._get_constraints(c, globdat, **kwargs)
        return c

    def ADVANCE(self, globdat):
        self._advance_step_constraints(globdat)

    @Model.save_config
    def configure(self, globdat, *, groups, dofs, values, dispIncr=None):
        # Validate input arguments
        check_list(self, groups)
        check_list(self, dofs)
        check_list(self, values)
        self._groups = groups
        self._dofs = dofs
        self._vals = values

        self._initDisp = self._vals
        if dispIncr is None:
            self._dispIncr = np.zeros(len(self._vals))
        else:
            self._dispIncr = dispIncr

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
