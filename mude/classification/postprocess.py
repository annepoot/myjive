import pandas as pd
import pyvista as pv

df = pd.read_csv("classification-data.csv")

for i in range(100, 1000, 100):
    reader = pv.get_reader("vtk-output/stiffness{}.vtu".format(i + 1))
    mesh = reader.read()
    mesh.plot(scalars="stiffness", cpos="xy")
