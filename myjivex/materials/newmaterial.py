from .isotropicmaterial import IsotropicMaterial
from .heterogeneousmaterial import HeterogeneousMaterial
from .deterioratedmaterial import DeterioratedMaterial

__all__ = ["new_material"]


def new_material(typ):

    if typ == "Isotropic":
        mat = IsotropicMaterial()
    elif typ == "Heterogeneous":
        mat = HeterogeneousMaterial()
    elif typ == "Deteriorated":
        mat = DeterioratedMaterial()
    else:
        raise ValueError(typ + " is not a valid material")

    return mat
