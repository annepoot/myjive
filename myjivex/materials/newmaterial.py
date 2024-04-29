from .isotropicmaterial import IsotropicMaterial
from .heterogeneousmaterial import HeterogeneousMaterial
from .deterioratedmaterial import DeterioratedMaterial

__all__ = ["new_material"]


def new_material(typ, name):
    if typ == "Isotropic":
        mat = IsotropicMaterial(name)
    elif typ == "Heterogeneous":
        mat = HeterogeneousMaterial(name)
    elif typ == "Deteriorated":
        mat = DeterioratedMaterial(name)
    else:
        raise ValueError(typ + " is not a valid material")

    return mat
