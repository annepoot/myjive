import numpy as np

from myjive.names import GlobNames as gn
from myjive.model.model import Model
from myjive.util.proputils import mandatory_list, optional_argument

__all__ = ["NeumannModel"]


class NeumannModel(Model):
    def GETCONSTRAINTS(self, c, globdat, **kwargs):
        c = self._get_constraints(c, globdat, **kwargs)
        return c

    def GETEXTFORCE(self, f_ext, globdat, **kwargs):
        f_ext = self._get_ext_force(f_ext, globdat, **kwargs)
        return f_ext

    def GETNEUMANNFORCE(self, f_neum, globdat, **kwargs):
        f_neum = self._get_neumann_force(f_neum, globdat, **kwargs)
        return f_neum

    def GETUNITFORCE(self, f_unit, globdat, **kwargs):
        f_unit = self._get_unit_force(f_unit, globdat, **kwargs)
        return f_unit

    def ADVANCE(self, globdat):
        self._advance_step(globdat)

    def configure(self, globdat, **props):

        # Get props
        self._groups = mandatory_list(self, props, "groups")
        self._dofs = mandatory_list(self, props, "dofs")
        self._vals = mandatory_list(self, props, "values")
        loadIncr = optional_argument(self, props, "loadIncr")

        self._initLoad = self._vals
        if loadIncr is None:
            self._loadIncr = np.zeros(len(self._vals))
        else:
            self._loadIncr = loadIncr

    def _get_constraints(self, c, globdat):
        ds = globdat[gn.DOFSPACE]
        for group, dof, val in zip(self._groups, self._dofs, self._vals):
            for node in globdat[gn.NGROUPS][group]:
                idof = ds.get_dof(node, dof)
                c.add_neumann(idof, val)
        return c

    def _get_ext_force(self, f_ext, globdat):
        if f_ext is None:
            f_ext = np.zeros(globdat[gn.DOFSPACE].dof_count())

        ds = globdat[gn.DOFSPACE]
        for group, dof, val in zip(self._groups, self._dofs, self._vals):
            for node in globdat[gn.NGROUPS][group]:
                idof = ds.get_dof(node, dof)
                f_ext[idof] += val

        return f_ext

    def _get_neumann_force(self, f_neum, globdat):
        if f_neum is None:
            f_neum = np.zeros(globdat[gn.DOFSPACE].dof_count())

        ds = globdat[gn.DOFSPACE]
        for group, dof, val in zip(self._groups, self._dofs, self._vals):
            for node in globdat[gn.NGROUPS][group]:
                idof = ds.get_dof(node, dof)
                f_neum[idof] += val

        return f_neum

    def _get_unit_force(self, f_unit, globdat):
        if f_unit is None:
            f_unit = np.zeros(globdat[gn.DOFSPACE].dof_count())

        ds = globdat[gn.DOFSPACE]
        for group, dof, incr in zip(self._groups, self._dofs, self._loadIncr):
            for node in globdat[gn.NGROUPS][group]:
                idof = ds.get_dof(node, dof)
                f_unit[idof] += incr

        return f_unit

    def _advance_step(self, globdat):
        self._vals = np.array(self._initLoad) + globdat[gn.TIMESTEP] * np.array(
            self._loadIncr
        )
