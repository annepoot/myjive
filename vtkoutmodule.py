import numpy as np

from module import *
from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

FILENAME = 'file'
TABLES = 'tables'


class VTKOutModule(Module):

    def init(self, props, globdat):
        self._fname = ''
        self._tnames = []

        if self._name in props:
            myprops = props.get(self._name)

            if FILENAME in myprops:
                self._fname = myprops[FILENAME]
            if TABLES in myprops:
                self._tnames = myprops[TABLES].strip('[').strip(']').split(',')

    def run(self, globdat):

        nodes = globdat[gn.NSET]
        elems = globdat[gn.ESET]
        disp = globdat[gn.STATE0]
        dofs = globdat[gn.DOFSPACE]
        types = dofs.get_types()

        self._write_tables(self._tnames, globdat)

        if self._fname:
            print('VTKOutModule: Writing output to file...')

            fname = self._fname + str(globdat[gn.TIMESTEP]) + '.vtu'

            with open(fname, 'w') as out:
                out.write('<VTKFile type="UnstructuredGrid"  version="0.1">\n')
                out.write('<UnstructuredGrid>\n')
                out.write('<Piece NumberOfPoints="' + str(len(nodes)) + '" NumberOfCells="' + str(len(elems)) + '">\n')
                out.write('<Points>\n')
                out.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
                for node in nodes:
                    out.write(' '.join(map(str, node.get_coords())) + '\n')
                out.write('</DataArray>\n')
                out.write('</Points>\n')
                out.write('<Cells>\n')
                out.write('<DataArray type="Int32" Name="connectivity" format="ascii">\n')
                for elem in elems:
                    out.write(' '.join(map(str, elem.get_nodes())) + '\n')
                out.write('</DataArray>\n')
                out.write('<DataArray type="Int32" Name="offsets" format="ascii">\n')
                i = 0
                for elem in elems:
                    i += len(elem.get_nodes())
                    out.write(str(i) + '\n')
                out.write('</DataArray>\n')
                out.write('<DataArray type="UInt8" Name="types" format="ascii">\n')
                for elem in elems:
                    assert (len(elem.get_nodes()) == 3)  # only writing type=5 (triangle) for now
                    out.write('5\n')
                out.write('</DataArray>\n')
                out.write('</Cells>\n')
                out.write('<PointData Vectors="fields">\n')
                out.write('<DataArray type="Float64" Name="U" NumberOfComponents="3" format="ascii">\n')
                for inode in range(len(nodes)):
                    idofs = dofs.get_dofs([inode], types)
                    out.write(' '.join(map(str, disp[idofs])))
                    out.write((3 - len(idofs)) * ' 0.0' + '\n')
                out.write('</DataArray>\n')
                for table in globdat[gn.TABLES].values():
                    for comp in table:
                        out.write(
                            '<DataArray type="Float64" Name="' + comp + '" NumberOfComponents="1" format="ascii">\n')
                        for inode in range(len(nodes)):
                            out.write(str(table[comp][inode]) + '\n')
                        out.write('</DataArray>\n')
                out.write('</PointData>\n')
                out.write('</Piece>\n')
                out.write('</UnstructuredGrid>\n')
                out.write('</VTKFile>\n')
        return 'ok'

    def shutdown(self, globdat):
        pass

    def _write_tables(self, table_names, globdat):
        nodes = globdat[gn.NSET]
        model = globdat[gn.MODEL]

        globdat[gn.TABLES] = {}

        for name in table_names:
            params = {}
            params[pn.TABLE] = {}
            params[pn.TABLENAME] = name
            params[pn.TABLEWEIGHTS] = np.zeros(len(nodes))

            model.take_action(act.GETTABLE, params, globdat)

            table = params[pn.TABLE]
            tbwts = params[pn.TABLEWEIGHTS]

            for col in table.values():
                col /= params[pn.TABLEWEIGHTS]

            globdat[gn.TABLES][name] = params[pn.TABLE]


def declare(factory):
    factory.declare_module('VTKOut', VTKOutModule)
