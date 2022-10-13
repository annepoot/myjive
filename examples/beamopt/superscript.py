import os, sys

import numpy as np
import pygmsh
from scipy.optimize import differential_evolution, Bounds

sys.path.append('../')

import proputils as pu
import main
from   names import GlobNames as gn

# Parameters defining the optimization problem
length    = 10.0
maxstress = 50.0
maxheight = 2.0
meshsize  = 0.2
nvars     = 5
rseed     = 10
bounds = Bounds ([0.1] * nvars, [maxheight] * nvars)

# Utility class for hiding FEM output
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

# Function defining the penalized objective function
def objective_function (vars):
    coords  = np.array ([[0.0, 0.0]])
    elsizes = []

    # Bottom
    for i in range(nvars-1):
        xp = float(i+1) * length/float(nvars-1)
        yp = 0.0
        coords = np.vstack([coords,[xp,yp]])

    # Top
    for i in range(nvars):
        xp = length - float(i) * length/float(nvars-1)
        yp = vars[i]
        coords = np.vstack ([coords, [xp,yp]])
        elsizes.append (yp*meshsize)

    #x = coords[:,0]
    #y = coords[:,1]

    # Gauss's formula for area of a general polygon (volume minimization)
    obj = 0.5*np.abs(np.dot(coords[:,0],np.roll(coords[:,1],1))-np.dot(coords[:,1],np.roll(coords[:,0],1)))

    with pygmsh.geo.Geometry() as geom:
        geom.add_polygon(
            coords.tolist(),
            mesh_size = list(reversed(elsizes)) + elsizes
        )
        mesh = geom.generate_mesh()

    # Read input file and set mesh
    props = pu.parse_file ('beam.pro')
    props['init']['mesh'] = mesh

    # Run Jive
    with HiddenPrints():
        globdat = main.jive (props)

    # Get the maximum horizontal stress seen by the beam
    stresses = globdat[gn.TABLES]['stress']['stress_xx']

    cons = np.max(np.abs(stresses)) - 50.0

    return obj + max(cons,0.0)

# Function for running the final design, including VTK output for Paraview
def run_final_design (vars):
    coords  = np.array ([[0.0, 0.0]])
    elsizes = []

    # Bottom
    for i in range(nvars-1):
        xp = float(i+1) * length/float(nvars-1)
        yp = 0.0
        coords = np.vstack([coords,[xp,yp]])

    # Top
    for i in range(nvars):
        xp = length - float(i) * length/float(nvars-1)
        yp = vars[i]
        coords = np.vstack ([coords, [xp,yp]])
        elsizes.append (yp*meshsize)

    with pygmsh.geo.Geometry() as geom:
        geom.add_polygon(
            coords.tolist(),
            mesh_size = list(reversed(elsizes)) + elsizes
        )
        mesh = geom.generate_mesh()

    # Run Jive
    props = pu.parse_file ('beam.pro')

    props['init']['mesh'] = mesh
    props['vtkout']['file'] = 'results'

    with HiddenPrints():
        globdat = main.jive (props)

# Run a meta-heuristic algorithm
result = differential_evolution (
                                 objective_function, 
                                 bounds=bounds,
                                 workers=4,
                                 seed=rseed, 
                                 popsize=10,
                                 maxiter=100,
                                 polish=False,
                                 disp=True
                                )

print('End of optimization. Final design', result.x, 'with objective function', result.fun)

# Run the optimum to obtain VTK results
run_final_design(result.x)
