from .heterogeneousmaterial import HeterogeneousMaterial
from jive.names import GlobNames as gn
from scipy.stats import norm
import numpy as np
import jive.util.proputils as pu

DETER_PROP = "deteriorations"
SCALE_PROP = "scale"
SEED_PROP = "seed"
LOC_X_PROP = "locX"
LOC_Y_PROP = "locY"
STD_X_PROP = "stdX"
STD_Y_PROP = "stdY"
SHAPE = "shape"
TYPE = "type"
INTSCHEME = "intScheme"

__all__ = ["DeterioratedMaterial"]


class DeterioratedMaterial(HeterogeneousMaterial):
    def configure(self, props, globdat):
        super().configure(props, globdat)

        self._ndet = int(props[DETER_PROP])
        self._detscale = float(props.get(SCALE_PROP, 1.0))

        self._seed = None
        if SEED_PROP in props:
            self._seed = int(props[SEED_PROP])

        # Get the location and standard deviation of the deterioration from the props file
        self._locx = props.get(LOC_X_PROP, "x")
        self._locy = props.get(LOC_Y_PROP, "y")
        self._stdx = props.get(STD_X_PROP, 1.0)
        self._stdy = props.get(STD_Y_PROP, 1.0)

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
