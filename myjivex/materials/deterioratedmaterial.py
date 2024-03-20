from .heterogeneousmaterial import HeterogeneousMaterial
from myjive.names import GlobNames as gn
from scipy.stats import norm
import numpy as np
from myjive.util.proputils import mandatory_argument, optional_argument
import myjive.util.proputils as pu


__all__ = ["DeterioratedMaterial"]


class DeterioratedMaterial(HeterogeneousMaterial):
    def configure(self, globdat, **props):
        super().configure(globdat, **props)

        # Get props
        self._ndet = mandatory_argument(self, props, "deteriorations")
        self._detscale = optional_argument(self, props, "scale", default=1.0)
        self._locx = optional_argument(self, props, "locX", default="x")
        self._locy = optional_argument(self, props, "locY", default="y")
        self._stdx = optional_argument(self, props, "stdX", default=1.0)
        self._stdy = optional_argument(self, props, "stdY", default=1.0)
        self._seed = optional_argument(self, props, "seed", default=None)

        self._detlocs = np.zeros((self._rank, self._ndet))
        self._detrads = np.zeros((self._rank, self._ndet))

        self._generate_deteriorations(globdat)

        globdat["detlocs"] = self._detlocs
        globdat["detrads"] = self._detrads

    def stiff_at_point(self, ipoint=None):
        return self._compute_stiff_matrix(ipoint)

    def mass_at_point(self, ipoint=None):
        return self._compute_mass_matrix(ipoint)

    def _get_E(self, ipoint=None):
        E = super()._get_E(ipoint)
        tiny = E * 1e-10
        scale = self._detscale * E

        # Subtract all deteriorations
        for i in range(self._ndet):
            det = norm.pdf(ipoint, loc=self._detlocs[:, i], scale=self._detrads[:, i])
            E -= np.prod(det) * scale

            if E <= 0:
                E = tiny
                break

        return E

    def _generate_deteriorations(self, globdat):
        elems = globdat[gn.ESET]

        np.random.seed(self._seed)

        for i in range(self._ndet):
            # randomly select an element
            ielem = np.random.randint(0, len(elems) - 1)
            elem = elems[ielem]
            inodes = elem.get_nodes()
            coords = globdat[gn.NSET].get_some_coords(inodes)

            center_coords = np.mean(coords, axis=1)

            # Generate the deterioration using the center coordinates of the element
            self._detlocs[0, i] = pu.evaluate(
                self._locx, center_coords, self._rank, extra_dict={"np": np}
            )
            self._detlocs[1, i] = pu.evaluate(
                self._locy, center_coords, self._rank, extra_dict={"np": np}
            )

            # Generate the standard deviations of the deterioration in two directions
            self._detrads[0, i] = pu.evaluate(
                self._stdx, center_coords, self._rank, extra_dict={"np": np}
            )
            self._detrads[1, i] = pu.evaluate(
                self._stdy, center_coords, self._rank, extra_dict={"np": np}
            )
