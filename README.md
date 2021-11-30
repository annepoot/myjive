# pyJive
A python Finite Element library inspired by jive

## Contributors
Iuri Rocha (i.rocha@tudelft.nl)
Frans van der Meer (f.p.vandermeer@tudelft.nl)

## Philosophy
The library works with *modules* and *models*. Modules define the flow of the program, e.g. the `solvermodule` defines that a matrix and RHS vector need to be assembled after which a linear system of equation is solved. Models define how the main task are implemented, e.g. the `elasticmodel` that assembles the stiffness matrix for an elasticity problem.

Modules and models interact through the `take_action` function. 

## Input
Problem-specific input is defined in a properties file (with .pro extension). The repository contains sample input files. 

## Examples
- `beam`: two-dimensional elastic analysis of a three-point bending problem 

## Important classes
- `dofspace`: degree of freedom management to provide mapping between dof index and node number and dof type 
- `shape`: abstract base clase that defines the element type 
- 

## Utility functions
- `declare.py`: this is where the available models and modules are defined (this is needed to be able to construct a problem dependent module-model scheme)
- `miniJive.py`: simple master script that parses an input file and calls `main.py`
- `main.py`: defines chain of modules

## Modules
- `initmodule`: initializes `globdat` object with global data including mesh and nodegroups
- `solvermodule`: assembles and solves system of equations
- `vtkoutmodule`: writes output to vtk file (compatible with paraview)

## Models
- `multimodel`: provides a fork from a single model into a list of models
- `elasticmodel`: assembles matrix for elasticity problem
- `poissonmodel`: assembles matrix for Poisson equation
- `dirimodel`: implements Dirichlet boundary conditions
- `neumannmodel`: implements Neumann boundary conditions
