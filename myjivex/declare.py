from myjive.names import GlobNames as gn
from myjive.app import ModuleFactory
from myjive.model import ModelFactory

from .models import (
    BarModel,
    DirichletModel,
    ElasticModel,
    LoadModel,
    NeumannModel,
    PoissonModel,
    SolidModel,
    TimoshenkoModel,
)

from .modules import VTKOutModule

__all__ = ["declare_all", "declare_models", "declare_modules"]


def declare_all(globdat):
    # Declare all core models and modules in one go
    declare_models(globdat)
    declare_modules(globdat)


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
