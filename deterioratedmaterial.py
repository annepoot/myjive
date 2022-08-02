from heterogeneousmaterial import HeterogeneousMaterial
from names import GlobNames as gn
import numpy as np

DETER_PROP = 'deteriorations'
SCALE_PROP = 'scale'
MAXSTD_PROP = 'stdMax'
SHAPE = 'shape'
TYPE = 'type'
INTSCHEME = 'intScheme'


class DeterioratedMaterial(HeterogeneousMaterial):

    def configure(self, props, globdat):

        super().configure(props)

        self._ndet = props[DETER_PROP]
        self._detscale = props[SCALE_PROP]
        self._maxstd = props[MAXSTD_PROP]

        self._detlocs = np.zeros((self._rank, self._ndet))
        self._detrads = np.zeros((self._rank, self._ndet))

        self._generate_deteriorations(self._ndet, globdat)

    def stiff_at_point(self, ipoint=None):
        return self._compute_stiff_matrix(ipoint)

    def mass_at_point(self, ipoint=None):
        return self._compute_mass_matrix(ipoint)

    def _get_E(self, ipoint=None):
        return pu.evaluate(self._E, ipoint, self._rank)

    def _generate_deteriorations(self, globdat):
        elems = globdat[gn.ESET]
        shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][TYPE], props[SHAPE][INTSCHEME])

        np.random.seed(0)

        for i in range(self._ndet):
            # randomly select an element
            ielem = np.random.randint(0, len(elems)-1)
            elem = elems[ielem]
            inodes = elem.get_nodes()
            coords = np.stack([globdat[gn.NSET][i].get_coords() for i in inodes], axis=1)[0:self._rank, :]

            # Put the center of the deterioration in the center of the element
            self._detlocs[0,i] = np.mean(coords[:,0])
            self._detlocs[1,i] = np.mean(coords[:,1])

            # Generate a random standard deviation in two directions
            self._detrads[0,i] = np.random.uniform(0, self._maxstd)
            self._detrads[1,i] = np.random.uniform(0, self._maxstd)
