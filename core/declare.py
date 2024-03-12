from jive.names import GlobNames as gn
from jive.app import ModuleFactory
from jive.model import ModelFactory

from core.models import (
    BarModel,
    DirichletModel,
    ElasticModel,
    LoadModel,
    NeumannModel,
    PoissonModel,
    SolidModel,
    TimoshenkoModel
)

from core.modules import VTKOutModule


def declare_models(globdat):
    factory = globdat.get(gn.MODELFACTORY, ModelFactory())

    BarModel.declare(factory)
    DirichletModel.declare(factory)
    NeumannModel.declare(factory)
    PoissonModel.declare(factory)
    ElasticModel.declare(factory)
    SolidModel.declare(factory)
    TimoshenkoModel.declare(factory)
    LoadModel.declare(factory)

    globdat[gn.MODELFACTORY] = factory


def declare_modules(globdat):
    factory = globdat.get(gn.MODULEFACTORY, ModuleFactory())

    VTKOutModule.declare(factory)

    globdat[gn.MODULEFACTORY] = factory
