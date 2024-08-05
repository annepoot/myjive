from .isotropicmaterial import IsotropicMaterial
import myjive.util.proputils as pu

__all__ = ["HeterogeneousMaterial"]


class HeterogeneousMaterial(IsotropicMaterial):
    @IsotropicMaterial.save_config
    def configure(
        self, globdat, *, rank, anmodel, E=1.0, nu=0.0, rho=0.0, area=1.0, params={}
    ):
        # Validate input arguments
        self._rank = rank
        self._anmodel = anmodel
        self._E = E
        self._nu = nu
        self._rho = rho
        self._area = area

        self._eval_dict = pu.get_core_eval_dict()
        self._eval_dict.update(params)

        self._strcount = self._rank * (self._rank + 1) // 2
        assert self._is_valid_anmodel(self._anmodel), (
            "Analysis model " + self._anmodel + " not valid for rank " + str(self._rank)
        )

    def stiff_at_point(self, ipoint=None):
        return self._compute_stiff_matrix(ipoint)

    def mass_at_point(self, ipoint=None):
        return self._compute_mass_matrix(ipoint)

    def _get_E(self, ipoint=None):
        return pu.evaluate(self._E, ipoint, self._rank, extra_dict=self._eval_dict)

    def _get_nu(self, ipoint=None):
        return pu.evaluate(self._nu, ipoint, self._rank, extra_dict=self._eval_dict)

    def _get_rho(self, ipoint=None):
        return pu.evaluate(self._rho, ipoint, self._rank, extra_dict=self._eval_dict)

    def _get_area(self, ipoint=None):
        return pu.evaluate(self._area, ipoint, self._rank, extra_dict=self._eval_dict)
