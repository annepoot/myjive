from names import GlobNames as gn

import model
import module
import shape

import multimodel
import barmodel
import dirimodel
import neumannmodel
import poissonmodel
import elasticmodel
import timoshenkomodel
import samplemodel

import initmodule
import solvermodule
import gaussianmodule
import outputmodule
import vtkoutmodule
import viewmodule

import paramshapes


def declare_models(globdat):
    factory = model.ModelFactory()

    multimodel.declare(factory)
    barmodel.declare(factory)
    dirimodel.declare(factory)
    neumannmodel.declare(factory)
    poissonmodel.declare(factory)
    elasticmodel.declare(factory)
    timoshenkomodel.declare(factory)
    samplemodel.declare(factory)

    globdat[gn.MODELFACTORY] = factory


def declare_modules(globdat):
    factory = module.ModuleFactory()

    initmodule.declare(factory)
    solvermodule.declare(factory)
    gaussianmodule.declare(factory)
    outputmodule.declare(factory)
    vtkoutmodule.declare(factory)
    viewmodule.declare(factory)

    globdat[gn.MODULEFACTORY] = factory


def declare_shapes(globdat):
    factory = shape.ShapeFactory()

    paramshapes.declare(factory)

    globdat[gn.SHAPEFACTORY] = factory
