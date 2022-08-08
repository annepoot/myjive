import pandas as pd
import pyvista as pv

df = pd.read_csv('regression-data.csv')

reader = pv.get_reader('vtk-output/stiffness500.vtu')
mesh = reader.read()
mesh.plot(scalars='stiffness', cpos='xy')
