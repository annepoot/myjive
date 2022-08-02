
MATERIAL = 'material'
TYPE = 'type'
RANK = 'rank'


def new_material(props):
    typ = props[TYPE]
    rank = int(props[RANK])

    if typ == 'Isotropic':
        from isotropicmaterial import IsotropicMaterial
        mat = IsotropicMaterial(rank)
    elif typ == 'Heterogeneous':
        from heterogeneousmaterial import HeterogeneousMaterial
        mat = HeterogeneousMaterial(rank)
    else:
        raise ValueError(typ + ' is not a valid material')

    return mat


class Material:

    def __init__(self, rank):
        pass

    def configure(self, props):
        pass

    def get_config(self):
        pass
