from jive.fem.names import GlobNames as gn

from jive.app import module
from jive.model import model

from core.models import (
    BarModel,
    DirichletModel,
    ElasticModel,
    LoadModel,
    NeumannModel,
    PoissonModel,
    SolidModel,
    TimoshenkoModel,
    XBarModel,
    XElasticModel,
    XPoissonModel,
)

from core.modules import (
    NonlinModule,
    VTKOutModule,
    LinBuckModule,
    LoadDispModule,
)


def declare_models(globdat):
    factory = globdat.get(gn.MODELFACTORY, model.ModelFactory())

    BarModel.declare(factory)
    DirichletModel.declare(factory)
    NeumannModel.declare(factory)
    PoissonModel.declare(factory)
    ElasticModel.declare(factory)
    SolidModel.declare(factory)
    TimoshenkoModel.declare(factory)
    LoadModel.declare(factory)
    XBarModel.declare(factory)
    XElasticModel.declare(factory)
    XPoissonModel.declare(factory)

    globdat[gn.MODELFACTORY] = factory


def declare_modules(globdat):
    factory = globdat.get(gn.MODULEFACTORY, module.ModuleFactory())

    NonlinModule.declare(factory)
    VTKOutModule.declare(factory)
    LinBuckModule.declare(factory)
    LoadDispModule.declare(factory)

    globdat[gn.MODULEFACTORY] = factory
