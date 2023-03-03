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
import gpmodel
import gpfmodel
import gpenkfmodel
import xbarmodel
import xelasticmodel
import xpoissonmodel

import gpinitmodule
import gpmodule
import gpexactmodule
import gpsamplermodule
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
    gpmodel.declare(factory)
    gpfmodel.declare(factory)
    gpenkfmodel.declare(factory)
    xbarmodel.declare(factory)
    xelasticmodel.declare(factory)
    xpoissonmodel.declare(factory)

    globdat[gn.MODELFACTORY] = factory


def declare_modules(globdat):
    factory = globdat.get(gn.MODULEFACTORY, module.ModuleFactory())

    gpinitmodule.declare(factory)
    gpmodule.declare(factory)
    gpexactmodule.declare(factory)
    gpsamplermodule.declare(factory)
    nonlinmodule.declare(factory)
    vtkoutmodule.declare(factory)
    linbuckmodule.declare(factory)
    loaddispmodule.declare(factory)

    globdat[gn.MODULEFACTORY] = factory
