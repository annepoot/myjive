import numpy as np

from module import *
from names import GlobNames as gn

class VTKOutModule (Module):

    def init ( self, props, globdat ):
        pass

    def run ( self, globdat ):
        print('VTKOutModule: Writing output to file...')

        fname = 'results' + str(globdat[gn.TIMESTEP]) + '.vtu'

        nodes = globdat[gn.NSET]
        elems = globdat[gn.ESET]
        disp = globdat[gn.STATE0]
        dofs = globdat[gn.DOFSPACE]
        types = dofs.get_types()

        with open (fname, 'w') as out:
            out.write('<VTKFile type="UnstructuredGrid"  version="0.1">\n')
            out.write('<UnstructuredGrid>\n')
            out.write('<Piece NumberOfPoints="'+str(len(nodes))+'" NumberOfCells="'+str(len(elems))+'">\n')
            out.write('<Points>\n')
            out.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
            for node in nodes:
                out.write(' '.join(map(str, node.get_coords()))+'\n')
            out.write('</DataArray>\n')
            out.write('</Points>\n')
            out.write('<Cells>\n')
            out.write('<DataArray type="Int32" Name="connectivity" format="ascii">\n')
            for elem in elems:
                out.write(' '.join(map(str, elem.get_nodes()))+'\n')
            out.write('</DataArray>\n')
            out.write('<DataArray type="Int32" Name="offsets" format="ascii">\n')
            i = 0
            for elem in elems:
                i += len(elem.get_nodes())
                out.write(str(i)+'\n')
            out.write('</DataArray>\n')
            out.write('<DataArray type="UInt8" Name="types" format="ascii">\n')
            for elem in elems:
                assert(len(elem.get_nodes())==3)   # only writing type=5 (triangle) for now
                out.write('5\n')
            out.write('</DataArray>\n')
            out.write('</Cells>\n')
            out.write('<PointData Vectors="U">\n')
            out.write('<DataArray type="Float64" Name="U" NumberOfComponents="3" format="ascii">\n')
            for inode in range(len(nodes)):
                idofs = dofs.get_dofs([inode],types)
                out.write(' '.join(map(str, disp[idofs])))
                out.write((3-len(idofs))*' 0.0' + '\n')
            out.write('</DataArray>\n')
            out.write('</PointData>\n')
            out.write('</Piece>\n')
            out.write('</UnstructuredGrid>\n')
            out.write('</VTKFile>\n')
        return 'ok'

    def shutdown ( self, globdat ):
        pass

def declare (factory):
    factory.declare_module ('VTKOut', VTKOutModule)
