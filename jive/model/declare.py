from jive.fem.names import GlobNames as gn

from jive.app import module
from jive.app import initmodule
from jive.app import outputmodule
from jive.gl import frameviewmodule
from jive.gl import graphmodule
from jive.gl import viewmodule
from jive.model import model
from jive.model import multimodel
from jive.fem import shape
from jive.fem import paramshapes
from jive.implicit import arclenmodule
from jive.implicit import solvermodule
from jive.implicit import linsolvemodule
from jive.solver import solver
from jive.solver import directsolver
from jive.solver import sparsecholesky
from jive.solver import cholmod
from jive.solver import iterativesolver
from jive.solver import cg
from jive.solver import preconditioner
from jive.solver import diagprecon
from jive.solver import icholprecon

import barmodel
import dirimodel
import neumannmodel
import poissonmodel
import elasticmodel
import solidmodel
import timoshenkomodel
import gpmodel
import gpfmodel
import gpenkfmodel
import xbarmodel
import xelasticmodel
import xpoissonmodel

import gpinitmodule
import gpexactmodule
import gpsamplermodule
import nonlinmodule
import vtkoutmodule
import linbuckmodule
import loaddispmodule

def declare_models(globdat):
    factory = model.ModelFactory()

    multimodel.declare(factory)
    barmodel.declare(factory)
    dirimodel.declare(factory)
    neumannmodel.declare(factory)
    poissonmodel.declare(factory)
    elasticmodel.declare(factory)
    solidmodel.declare(factory)
    timoshenkomodel.declare(factory)
    gpmodel.declare(factory)
    gpfmodel.declare(factory)
    gpenkfmodel.declare(factory)
    xbarmodel.declare(factory)
    xelasticmodel.declare(factory)
    xpoissonmodel.declare(factory)

    globdat[gn.MODELFACTORY] = factory


def declare_modules(globdat):
    factory = module.ModuleFactory()

    initmodule.declare(factory)
    solvermodule.declare(factory)
    linsolvemodule.declare(factory)
    gpinitmodule.declare(factory)
    gpexactmodule.declare(factory)
    gpsamplermodule.declare(factory)
    nonlinmodule.declare(factory)
    arclenmodule.declare(factory)
    outputmodule.declare(factory)
    vtkoutmodule.declare(factory)
    linbuckmodule.declare(factory)
    viewmodule.declare(factory)
    frameviewmodule.declare(factory)
    loaddispmodule.declare(factory)
    graphmodule.declare(factory)

    globdat[gn.MODULEFACTORY] = factory


def declare_shapes(globdat):
    factory = shape.ShapeFactory()

    paramshapes.declare(factory)

    globdat[gn.SHAPEFACTORY] = factory


def declare_solvers(globdat):
    factory = solver.SolverFactory()

    iterativesolver.declare(factory)
    cg.declare(factory)
    directsolver.declare(factory)
    sparsecholesky.declare(factory)
    cholmod.declare(factory)

    globdat[gn.SOLVERFACTORY] = factory


def declare_precons(globdat):
    factory = preconditioner.PreconFactory()

    diagprecon.declare(factory)
    icholprecon.declare(factory)

    globdat[gn.PRECONFACTORY] = factory
