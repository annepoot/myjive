# BFEM
A Bayesian Finite Element Method.
The fundamental ideas setup of BFEM can be found [here](https://hal.archives-ouvertes.fr/hal-03719275/).

## Contributors
Anne Poot (a.poot-1@tudelft.nl), Pierre Kerfriden (pierre.kerfriden@minesparis.psl.eu), Iuri Rocha (i.rocha@tudeft.nl), Frans van der Meer (f.p.vandermeer@tudelft.nl)

## Examples
All BFEM-related examples are stored in a directory that starts with `gp`.
- `gpbeam`: BFEM applied to a three-point bending problem (similar to `beam`)
- `gpbar`: BFEM applied to an axially loaded bar problem (similar to `bar`)
- `gppoisson`: BFEM applied to a poisson problem
- `gptapered`: similar to `gpbar`, but the bar is tapered, and inhomogeneous Dirichlet boundary conditions have been applied
- `gpensemble`: similar to `gptapered`, but unlike the other examples an Ensemble Kalman Filter is used instead of an exact Kalman update

## GP Modules
- `gpinitmodule`: initializes coarse mesh (child class of `initmodule`)
- `gpsolvermodule`: assembles the prior and performs the Kalman update to obtain the posterior (child class of `solvermodule`)
- `gpsamplermodule`: draws samples from the prior and posterior

## GP Models
- `gpmodel`: fundamental BFEM model, performing all calculations exactly and making no assumptions on the prior
- `gpfmodel`: assumes that the prior can be defined on $`K u`$ (i.e. $`f`$) instead of $`u`$ (child class of `gpmodel`)
- `gpenkfmodel`: implements the posterior update using an Ensemble Kalman Filter instead of the exact Kalman update (child class of `gpmodel`)

# pyJive
A python Finite Element library inspired by jive

## Contributors
Iuri Rocha (i.rocha@tudelft.nl), Frans van der Meer (f.p.vandermeer@tudelft.nl), Andres Martinez Colan

## Philosophy
The library works with *modules* and *models*. Modules define the flow of the program, e.g. the `solvermodule` defines that a matrix and RHS vector need to be assembled after which a linear system of equation is solved. Models define how the main task are implemented, e.g. the `elasticmodel` that assembles the stiffness matrix for an elasticity problem.

Modules and models interact through the `take_action` function. 

## Input
Problem-specific input is defined in a properties file (with .pro extension). The repository contains sample input files. 

## Examples
- `beam`: two-dimensional elastic analysis of a three-point bending problem
- `bar`: one-dimensional elastic analysis of an axially loaded bar problem
- `timoshenko`: one-dimensional elastic analysis of a timoshenko beam problem

## Important classes
- `dofspace`: degree of freedom management to provide mapping between dof index and node number and dof type 
- `shape`: abstract base class that defines the element type 
- `constrainer`: class that handles the application of constrains to the as-assembled stiffness matrix and force vector, and returns constrained versions of them ready for solving 

## Implemented shapes
Generates shape functions, gradients and weights. The following classes are implemented in `paramshapes.py`:
- `Tri3Shape`: 3-node linear triangular shape for 1-point Gauss integration
- `Line2Shape`: 2-node linear line shape for 1-point and 2-point Gauss integration
- `Line3Shape`: 3-node quadratic line shape for 2-point Gauss integration

## Utility functions
- `declare.py`: this is where the available models and modules are defined (this is needed to be able to construct a problem dependent module-model scheme)
- `miniJive.py`: simple master script that parses an input .pro file through `proputils.py` and calls `main.py`
- `main.py`: defines chain of modules

## Modules
- `initmodule`: initializes `globdat` object with global data including mesh and nodegroups. Accepts mesh files from `gmsh` and `meshio` as well as manually generated mesh files (see syntaxis on `initmodule.py`)
- `solvermodule`: assembles and solves system of equations on a given number of steps
- `vtkoutmodule`: writes output to vtk file (compatible with paraview)
- `viewmodule`: visualization of full field data
- `loaddispmodule`: stores load displacement data for specified nodegroups in globdat
- `graphmodule`: plots data for instance data stored by loaddispmodule

## Models
- `multimodel`: provides a fork from a single model into a list of models
- `elasticmodel`: assembles matrix for elasticity problem
- `barmodel`: assembles matrix for a 1D bar problem
- `timoshenkomodel`: assembles matrix for a Timoshenko beam problem
- `poissonmodel`: assembles matrix for Poisson equation
- `dirimodel`: implements Dirichlet boundary conditions
- `neumannmodel`: implements Neumann boundary conditions
