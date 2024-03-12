from jive.names import GlobNames as gn

from jive.app import ModuleFactory, InitModule, OutputModule
from jive.implicit import LinsolveModule, SolverModule
from jive.model import ModelFactory
from jive.solver import (
    SolverFactory,
    CGSolver,
    CholmodSolver,
    DirectSolver,
    IterativeSolver,
    SparseCholeskySolver,
    PreconFactory,
    DiagPrecon,
    ICholPrecon,
    IdPrecon,
)
from jive.fem import (
    ShapeFactory,
    Tri3Shape,
    Tri6Shape,
    Quad4Shape,
    Quad9Shape,
    Line2Shape,
    Line3Shape,
)

import core.declare as core_declare


def declare_all(globdat):
    # Declare all standard jive models and modules in one go
    declare_models(globdat)
    declare_modules(globdat)
    declare_solvers(globdat)
    declare_precons(globdat)
    declare_shapes(globdat)

    # Declare all custom models and modules as well
    core_declare.declare_models(globdat)
    core_declare.declare_modules(globdat)


def declare_models(globdat):
    pass


def declare_modules(globdat):
    factory = globdat.get(gn.MODULEFACTORY, ModuleFactory())

    InitModule.declare(factory)

    OutputModule.declare(factory)
    LinsolveModule.declare(factory)
    SolverModule.declare(factory)

    globdat[gn.MODULEFACTORY] = factory


def declare_solvers(globdat):
    factory = globdat.get(gn.SOLVERFACTORY, SolverFactory())

    CGSolver.declare(factory)
    CholmodSolver.declare(factory)
    DirectSolver.declare(factory)
    IterativeSolver.declare(factory)
    SparseCholeskySolver.declare(factory)

    globdat[gn.SOLVERFACTORY] = factory


def declare_precons(globdat):
    factory = globdat.get(gn.PRECONFACTORY, PreconFactory())

    DiagPrecon.declare(factory)
    ICholPrecon.declare(factory)
    IdPrecon.declare(factory)

    globdat[gn.PRECONFACTORY] = factory


def declare_shapes(globdat):
    factory = globdat.get(gn.SHAPEFACTORY, ShapeFactory())

    Tri3Shape.declare(factory)
    Tri6Shape.declare(factory)
    Quad4Shape.declare(factory)
    Quad9Shape.declare(factory)
    Line2Shape.declare(factory)
    Line3Shape.declare(factory)

    globdat[gn.SHAPEFACTORY] = factory
