from jive.fem.names import GlobNames as gn

from jive.app import module
from jive.model import model

from core.models import barmodel
from core.models import dirimodel
from core.models import neumannmodel
from core.models import poissonmodel
from core.models import elasticmodel
from core.models import solidmodel
from core.models import timoshenkomodel
from core.models import loadmodel
from core.models import xbarmodel
from core.models import xelasticmodel
from core.models import xpoissonmodel

from core.modules import nonlinmodule
from core.modules import vtkoutmodule
from core.modules import linbuckmodule
from core.modules import loaddispmodule

def declare_models(globdat):
    factory = globdat.get(gn.MODELFACTORY, model.ModelFactory())

    barmodel.declare(factory)
    dirimodel.declare(factory)
    neumannmodel.declare(factory)
    poissonmodel.declare(factory)
    elasticmodel.declare(factory)
    solidmodel.declare(factory)
    timoshenkomodel.declare(factory)
    loadmodel.declare(factory)
    xbarmodel.declare(factory)
    xelasticmodel.declare(factory)
    xpoissonmodel.declare(factory)

    globdat[gn.MODELFACTORY] = factory


def declare_modules(globdat):
    factory = globdat.get(gn.MODULEFACTORY, module.ModuleFactory())

    nonlinmodule.declare(factory)
    vtkoutmodule.declare(factory)
    linbuckmodule.declare(factory)
    loaddispmodule.declare(factory)

    globdat[gn.MODULEFACTORY] = factory
