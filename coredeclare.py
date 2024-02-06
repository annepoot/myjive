from jive.fem.names import GlobNames as gn

from jive.app import module
from jive.model import model

import barmodel
import dirimodel
import neumannmodel
import poissonmodel
import elasticmodel
import solidmodel
import timoshenkomodel
import loadmodel
import xbarmodel
import xelasticmodel
import xpoissonmodel

import nonlinmodule
import vtkoutmodule
import linbuckmodule
import loaddispmodule

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
