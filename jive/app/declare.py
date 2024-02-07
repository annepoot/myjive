from jive.fem.names import GlobNames as gn

from jive.app import module
from jive.app import initmodule
from jive.app import outputmodule

import jive.app.declare as app_declare
import jive.fem.declare as fem_declare
import jive.gl.declare as gl_declare
import jive.implicit.declare as implicit_declare
import jive.model.declare as model_declare
import jive.solver.declare as solver_declare
import core.declare as core_declare


def declare_all(globdat):
    # Declare all standard jive models and modules in one go
    app_declare.declare_modules(globdat)
    gl_declare.declare_modules(globdat)
    implicit_declare.declare_modules(globdat)
    model_declare.declare_models(globdat)
    fem_declare.declare_shapes(globdat)
    solver_declare.declare_solvers(globdat)
    solver_declare.declare_precons(globdat)

    # Declare all custom models and modules as well
    core_declare.declare_models(globdat)
    core_declare.declare_modules(globdat)


def declare_modules(globdat):
    factory = globdat.get(gn.MODULEFACTORY, module.ModuleFactory())

    initmodule.declare(factory)
    outputmodule.declare(factory)

    globdat[gn.MODULEFACTORY] = factory
