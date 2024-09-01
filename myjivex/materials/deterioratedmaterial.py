from .heterogeneousmaterial import HeterogeneousMaterial
from myjive.names import GlobNames as gn
from scipy.stats import norm
import numpy as np
import myjive.util.proputils as pu


__all__ = ["DeterioratedMaterial"]


class DeterioratedMaterial(HeterogeneousMaterial):
    @HeterogeneousMaterial.save_config
    def configure(
        self,
        globdat,
        *,
        deteriorations,
        scale=1.0,
        locX="x",
        locY="y",
        stdX=1.0,
        stdY=1.0,
        seed=None,
        **otherprops
    ):
        super().configure(globdat, **otherprops)

        # Get props
        self._ndet = deteriorations
        self._detscale = scale
        self._locx, self._locy = locX, locY
        self._stdx, self._stdy = stdX, stdY
        self._seed = seed

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
            inodes = elems[ielem]
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
