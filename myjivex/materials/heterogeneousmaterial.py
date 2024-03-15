from .isotropicmaterial import IsotropicMaterial
from .isotropicmaterial import (
    E_PROP,
    NU_PROP,
    RHO_PROP,
    AREA_PROP,
    ANMODEL_PROP,
)
import myjive.util.proputils as pu

__all__ = ["HeterogeneousMaterial"]


class HeterogeneousMaterial(IsotropicMaterial):
    def configure(self, props, globdat):
        self._anmodel = props.get(ANMODEL_PROP, self._anmodel)
        assert self._is_valid_anmodel(self._anmodel), (
            "Analysis model " + self._anmodel + " not valid for rank " + str(self._rank)
        )

        self._E = props.get(E_PROP, self._E)
        self._nu = props.get(NU_PROP, self._nu)
        self._rho = props.get(RHO_PROP, self._rho)

        self._E = pu.soft_cast(self._E, float)
        self._nu = pu.soft_cast(self._nu, float)
        self._rho = pu.soft_cast(self._rho, float)

        if self._rank == 1:
            self._area = props.get(AREA_PROP, self._area)
            self._area = pu.soft_cast(self._area, float)

    def stiff_at_point(self, ipoint=None):
        return self._compute_stiff_matrix(ipoint)

    def mass_at_point(self, ipoint=None):
        return self._compute_mass_matrix(ipoint)

    def _get_E(self, ipoint=None):
        return pu.evaluate(self._E, ipoint, self._rank)

    def _get_nu(self, ipoint=None):
        return pu.evaluate(self._nu, ipoint, self._rank)

    def _get_rho(self, ipoint=None):
        return pu.evaluate(self._rho, ipoint, self._rank)

    def _get_area(self, ipoint=None):
        return pu.evaluate(self._area, ipoint, self._rank)
