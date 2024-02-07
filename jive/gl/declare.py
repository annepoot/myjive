from jive.fem.names import GlobNames as gn

from jive.app import module
from jive.gl import frameviewmodule
from jive.gl import graphmodule
from jive.gl import viewmodule


def declare_modules(globdat):
    factory = globdat.get(gn.MODULEFACTORY, module.ModuleFactory())

    viewmodule.declare(factory)
    frameviewmodule.declare(factory)
    graphmodule.declare(factory)

    globdat[gn.MODULEFACTORY] = factory
