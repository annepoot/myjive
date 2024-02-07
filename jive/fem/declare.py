from jive.fem.names import GlobNames as gn

from jive.fem import shape
from jive.fem import paramshapes


def declare_shapes(globdat):
    factory = globdat.get(gn.SHAPEFACTORY, shape.ShapeFactory())

    paramshapes.declare(factory)

    globdat[gn.SHAPEFACTORY] = factory
