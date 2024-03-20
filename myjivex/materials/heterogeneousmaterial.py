from .isotropicmaterial import IsotropicMaterial
from myjive.util.proputils import mandatory_argument, optional_argument
import myjive.util.proputils as pu
import numpy as np

__all__ = ["HeterogeneousMaterial"]


class HeterogeneousMaterial(IsotropicMaterial):

    def configure(self, globdat, **props):

        # Get props
        self._rank = mandatory_argument(self, props, "rank")
        self._anmodel = mandatory_argument(self, props, "anmodel")
        self._E = optional_argument(self, props, "E", default=1.0)
        self._nu = optional_argument(self, props, "nu", default=0.0)
        self._rho = optional_argument(self, props, "rho", default=0.0)
        self._area = optional_argument(self, props, "area", default=1.0)

        self._strcount = self._rank * (self._rank + 1) // 2
        assert self._is_valid_anmodel(self._anmodel), (
            "Analysis model " + self._anmodel + " not valid for rank " + str(self._rank)
        )

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
