from jive.fem.names import GlobNames as gn

from jive.solver import solver
from jive.solver import directsolver
from jive.solver import sparsecholesky
from jive.solver import cholmod
from jive.solver import iterativesolver
from jive.solver import cg

from jive.solver import preconditioner
from jive.solver import idprecon
from jive.solver import diagprecon
from jive.solver import icholprecon

def declare_solvers(globdat):
    factory = globdat.get(gn.SOLVERFACTORY, solver.SolverFactory())

    iterativesolver.declare(factory)
    cg.declare(factory)
    directsolver.declare(factory)
    sparsecholesky.declare(factory)
    cholmod.declare(factory)

    globdat[gn.SOLVERFACTORY] = factory


def declare_precons(globdat):
    factory = globdat.get(gn.PRECONFACTORY, preconditioner.PreconFactory())

    idprecon.declare(factory)
    diagprecon.declare(factory)
    icholprecon.declare(factory)

    globdat[gn.PRECONFACTORY] = factory
