from .isotropicmaterial import IsotropicMaterial
from .heterogeneousmaterial import HeterogeneousMaterial
from .deterioratedmaterial import DeterioratedMaterial
from .material import TYPE, RANK

__all__ = ["new_material"]


def new_material(props):
    typ = props[TYPE]
    rank = int(props[RANK])

    if typ == "Isotropic":
        mat = IsotropicMaterial(rank)
    elif typ == "Heterogeneous":
        mat = HeterogeneousMaterial(rank)
    elif typ == "Deteriorated":
        mat = DeterioratedMaterial(rank)
    else:
        raise ValueError(typ + " is not a valid material")

    return mat
