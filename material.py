from isotropicmaterial import IsotropicMaterial
from heterogeneousmaterial import HeterogeneousMaterial

MATERIAL = 'material'
TYPE = 'type'
RANK = 'rank'


def new_material(name, props):
    myprops = props.get(MATERIAL)
    typ = myprops[TYPE]
    rank = myprops[RANK]

    if typ == 'Isotropic':
        mat = IsotropicMaterial(rank)
    elif typ == 'Heterogeneous':
        mat = HeterogeneousMaterial(rank)
    else:
        raise ValueError(typ + ' is not a valid material')

    mat.configure(props)

    return mat


class Material:

    def __init__(self, rank):
        pass

    def configure(self, props):
        pass

    def get_config(self):
        pass
