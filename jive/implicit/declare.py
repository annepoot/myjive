from jive.fem.names import GlobNames as gn

from jive.app import module
from jive.implicit import arclenmodule
from jive.implicit import solvermodule
from jive.implicit import linsolvemodule


def declare_modules(globdat):
    factory = globdat.get(gn.MODULEFACTORY, module.ModuleFactory())

    solvermodule.declare(factory)
    linsolvemodule.declare(factory)
    arclenmodule.declare(factory)

    globdat[gn.MODULEFACTORY] = factory
