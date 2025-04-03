from .names import GlobNames as gn

__all__ = [
    "declare_all",
    "declare_models",
    "declare_modules",
    "declare_solvers",
    "declare_precons",
    "declare_shapes",
]


def declare_all(globdat, extra_declares=[]):
    # Declare all standard jive models and modules in one go
    declare_models(globdat)
    declare_modules(globdat)
    declare_solvers(globdat)
    declare_precons(globdat)
    declare_shapes(globdat)

    # Declare all custom models and modules as well
    for extra_declare in extra_declares:
        extra_declare(globdat)


def declare_models(globdat):
    from .model import ModelFactory

    factory = globdat.get(gn.MODELFACTORY, ModelFactory())

    globdat[gn.MODELFACTORY] = factory


def declare_modules(globdat):
    from .app import ModuleFactory, InitModule, OutputModule
    from .implicit import LinsolveModule, SolverModule

    factory = globdat.get(gn.MODULEFACTORY, ModuleFactory())

    InitModule.declare(factory)

    OutputModule.declare(factory)
    LinsolveModule.declare(factory)
    SolverModule.declare(factory)

    globdat[gn.MODULEFACTORY] = factory


def declare_solvers(globdat):
    from .solver import (
        SolverFactory,
        CGSolver,
        CholmodSolver,
        DirectSolver,
        IterativeSolver,
        SparseCholeskySolver,
    )

    factory = globdat.get(gn.SOLVERFACTORY, SolverFactory())

    CGSolver.declare(factory)
    CholmodSolver.declare(factory)
    DirectSolver.declare(factory)
    IterativeSolver.declare(factory)
    SparseCholeskySolver.declare(factory)

    globdat[gn.SOLVERFACTORY] = factory


def declare_precons(globdat):
    from .solver import PreconFactory, DiagPrecon, ICholPrecon, IdPrecon

    factory = globdat.get(gn.PRECONFACTORY, PreconFactory())

    DiagPrecon.declare(factory)
    ICholPrecon.declare(factory)
    IdPrecon.declare(factory)

    globdat[gn.PRECONFACTORY] = factory


def declare_shapes(globdat):
    from .fem import (
        ShapeFactory,
        Tri3Shape,
        Tri6Shape,
        Quad4Shape,
        Quad9Shape,
        Line2Shape,
        Line3Shape,
    )

    factory = globdat.get(gn.SHAPEFACTORY, ShapeFactory())

    Tri3Shape.declare(factory)
    Tri6Shape.declare(factory)
    Quad4Shape.declare(factory)
    Quad9Shape.declare(factory)
    Line2Shape.declare(factory)
    Line3Shape.declare(factory)

    globdat[gn.SHAPEFACTORY] = factory
