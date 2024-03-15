import numpy as np

from myjive.app import Module
from myjive.names import GlobNames as gn
import myjive.util.proputils as pu
from myjive.util import to_xtable

FILENAME = "file"
TABLES = "tables"

__all__ = ["VTKOutModule"]


class VTKOutModule(Module):
    def init(self, props, globdat):
        self._fname = ""
        self._tnames = []

        if self._name in props:
            myprops = props.get(self._name)

            if FILENAME in myprops:
                self._fname = myprops[FILENAME]
            if TABLES in myprops:
                self._tnames = pu.parse_list(myprops[TABLES])

    def run(self, globdat):
        nodes = globdat[gn.NSET]
        elems = globdat[gn.ESET]
        disp = globdat[gn.STATE0]
        dofs = globdat[gn.DOFSPACE]
        types = dofs.get_types()

        self._write_tables(self._tnames, globdat)

        if self._fname:
            print("VTKOutModule: Writing output to file...")

            fname = self._fname + str(globdat[gn.TIMESTEP]) + ".vtu"

            with open(fname, "w") as out:
                out.write('<VTKFile type="UnstructuredGrid"  version="0.1">\n')
                out.write("<UnstructuredGrid>\n")
                out.write(
                    '<Piece NumberOfPoints="'
                    + str(len(nodes))
                    + '" NumberOfCells="'
                    + str(len(elems))
                    + '">\n'
                )
                out.write("<Points>\n")
                out.write(
                    '<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n'
                )
                for node in nodes:
                    coords3d = np.zeros(3)
                    coords3d[: node.rank()] = node.get_coords()
                    out.write(" ".join(map(str, coords3d)) + "\n")
                out.write("</DataArray>\n")
                out.write("</Points>\n")
                out.write("<Cells>\n")
                out.write(
                    '<DataArray type="Int32" Name="connectivity" format="ascii">\n'
                )
                for elem in elems:
                    out.write(" ".join(map(str, elem.get_nodes())) + "\n")
                out.write("</DataArray>\n")
                out.write('<DataArray type="Int32" Name="offsets" format="ascii">\n')
                i = 0
                for elem in elems:
                    i += len(elem.get_nodes())
                    out.write(str(i) + "\n")
                out.write("</DataArray>\n")
                out.write('<DataArray type="UInt8" Name="types" format="ascii">\n')
                for elem in elems:
                    assert (
                        len(elem.get_nodes()) == 3
                    )  # only writing type=5 (triangle) for now
                    out.write("5\n")
                out.write("</DataArray>\n")
                out.write("</Cells>\n")
                out.write('<PointData Vectors="fields">\n')
                out.write(
                    '<DataArray type="Float64" Name="U" NumberOfComponents="3" format="ascii">\n'
                )
                for inode in range(len(nodes)):
                    idofs = dofs.get_dofs([inode], types)
                    out.write(" ".join(map(str, disp[idofs])))
                    out.write((3 - len(idofs)) * " 0.0" + "\n")
                out.write("</DataArray>\n")
                for name, table in globdat[gn.TABLES].items():
                    for comp in table:
                        if comp == "":
                            out.write(
                                '<DataArray type="Float64" Name="'
                                + name
                                + '" NumberOfComponents="1" format="ascii">\n'
                            )
                        else:
                            out.write(
                                '<DataArray type="Float64" Name="'
                                + name
                                + "_"
                                + comp
                                + '" NumberOfComponents="1" format="ascii">\n'
                            )
                        for inode in range(len(nodes)):
                            out.write(str(table[comp][inode]) + "\n")
                        out.write("</DataArray>\n")
                out.write("</PointData>\n")
                out.write("</Piece>\n")
                out.write("</UnstructuredGrid>\n")
                out.write("</VTKFile>\n")
        return "ok"

    def shutdown(self, globdat):
        pass

    def _write_tables(self, table_names, globdat):
        models = globdat[gn.MODELS]

        if gn.TABLES not in globdat:
            globdat[gn.TABLES] = {}

        for name in table_names:
            table = None
            tbwts = None

            for model in self.get_relevant_models("GETTABLE", models):
                table, tbwts = model.GETTABLE(name, table, tbwts, globdat)

            to_xtable(table)

            for jcol in range(table.column_count()):
                values = table.get_col_values(None, jcol)
                table.set_col_values(None, jcol, values / tbwts)

            table.to_table()
            globdat[gn.TABLES][name] = table
