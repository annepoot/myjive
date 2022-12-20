from jive.fem.names import GlobNames as gn

from jive.model import model
from jive.model import multimodel

def declare_models(globdat):
    factory = globdat.get(gn.MODELFACTORY, model.ModelFactory())

    multimodel.declare(factory)

    globdat[gn.MODELFACTORY] = factory
