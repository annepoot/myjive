import numpy as np

from model import Model, ModelFactory
from module import Module, ModuleFactory
from node import Node
from element import Element
from dofspace import DofSpace

from names import GlobNames as gn
from names import PropNames as prn

import proputils as pu

MESH = 'mesh'
TYPE = 'type'
FILE = 'file'


class InitModule(Module):
    def init(self, props, globdat):
        myprops = props.get(self._name)

        if not myprops:
            raise KeyError('Properties for InitModule not found')

        # Initialize some parameters

        self._ctol = 1.e-5
        globdat[gn.NGROUPS] = {}
        globdat[gn.EGROUPS] = {}
        modelfac = globdat[gn.MODELFACTORY]

        # Get the appropriate model for this module
        self._modelname = myprops.get(gn.MODEL, gn.MODEL)
        modelprops = props[self._modelname]

        # Initialize DofSpace
        print('InitModule: Creating DofSpace...')
        globdat[gn.DOFSPACE] = DofSpace()

        # Read mesh
        meshprops = myprops[MESH]

        if 'gmsh' in meshprops[TYPE]:
            self._read_gmsh(meshprops[FILE], globdat)
        elif 'manual' in meshprops[TYPE]:
            self._read_mesh(meshprops[FILE], globdat)
        elif 'meshio' in meshprops[TYPE]:
            self._read_meshio(meshprops[FILE], globdat)
        elif 'geo' in meshprops[TYPE]:
            self._read_geo(meshprops[FILE], globdat)
        else:
            raise KeyError('InitModule: Mesh input type unknown')

        # Create node groups
        if gn.NGROUPS in myprops:
            print('InitModule: Creating node groups...')
            groups = pu.parse_list(myprops[gn.NGROUPS])
            self._create_ngroups(groups, myprops, globdat)

        # Create element groups
        if gn.EGROUPS in myprops:
            print('InitModule: Creating element groups...')
            groups = pu.parse_list(myprops[gn.EGROUPS])
            self._create_egroups(groups, globdat)

        # Initialize model
        print('InitModule: Creating model...')
        m = modelfac.get_model(modelprops[prn.TYPE], self._modelname)
        m.configure(modelprops, globdat)
        globdat[self._modelname] = m

    def run(self, globdat):
        return ('ok')

    def shutdown(self, globdat):
        pass

    def _read_gmsh(self, fname, globdat):
        print('InitModule: Reading mesh file', fname, '...')

        if not fname.endswith('.msh'):
            raise RuntimeError('Unexpected mesh file extension')

        nodes = []
        elems = []

        parse_nodes = False
        parse_elems = False
        eltype = 0
        nnodes = 0

        with open(fname) as msh:
            for line in msh:
                sp = line.split()

                if '$Nodes' in line:
                    parse_nodes = True

                elif '$Elements' in line:
                    parse_nodes = False
                    parse_elems = True

                elif parse_nodes and len(sp) > 1:
                    if len(sp[1:]) != 3:
                        raise SyntaxError('InitModule: Three coordinates per node are expected')
                    coords = np.array(sp[1:], dtype=np.float64)
                    nodes.append(Node(coords))

                elif parse_elems and len(sp) > 1:
                    if eltype == 0:
                        eltype = int(sp[1])
                        if eltype == 1:
                            globdat[gn.MESHSHAPE] = 'Line2'
                            globdat[gn.MESHRANK] = 1
                            nnodes = 2
                        elif eltype == 2:
                            globdat[gn.MESHSHAPE] = 'Triangle3'
                            globdat[gn.MESHRANK] = 2
                            nnodes = 3
                        elif eltype == 3:
                            globdat[gn.MESHSHAPE] = 'Quad4'
                            globdat[gn.MESHRANK] = 2
                            nnodes = 4
                        elif eltype == 4:
                            globdat[gn.MESHSHAPE] = 'Tet4'
                            globdat[gn.MESHRANK] = 3
                            nnodes = 4
                        elif eltype == 5:
                            globdat[gn.MESHSHAPE] = 'Brick8'
                            globdat[gn.MESHRANK] = 3
                            nnodes = 8
                        elif eltype == 8:
                            globdat[gn.MESHSHAPE] = 'Line3'
                            globdat[gn.MESHRANK] = 1
                            nnodes = 3
                        else:
                            raise SyntaxError('InitModule: Unsupported element type')
                    elif eltype != int(sp[1]):
                        raise SyntaxError('InitModule: Only one element type per mesh is supported')
                    inodes = np.array(sp[3 + int(sp[2]):], dtype=np.int32) - 1
                    if len(inodes) != nnodes:
                        raise SyntaxError('InitModule: Could not read element with incorrect number of nodes')
                    elems.append(Element(inodes))

        globdat[gn.NSET] = nodes
        globdat[gn.ESET] = elems

        globdat[gn.NGROUPS]['all'] = [*range(len(nodes))]
        globdat[gn.EGROUPS]['all'] = [*range(len(elems))]

    def _read_mesh(self, fname, globdat):
        print('InitModule: Reading manual mesh file', fname, '...')

        nodes = []
        elems = []
        parse_nodes = False
        parse_elems = False

        with open(fname) as msh:
            for line in msh:
                sp = line.split()

                if 'nodes' in line:
                    parse_nodes = True
                    parse_elems = False

                elif 'elements' in line or 'elems' in line:
                    parse_nodes = False
                    parse_elems = True

                elif parse_nodes and len(sp) > 1:
                    coords = np.array(sp[1:], dtype=np.float64)
                    nodes.append(Node(coords))

                elif parse_elems and len(sp) > 0:
                    connectivity = np.array(sp, dtype=np.int16)
                    elems.append(Element(connectivity))

        globdat[gn.NSET] = nodes
        globdat[gn.ESET] = elems

        globdat[gn.NGROUPS]['all'] = [*range(len(nodes))]
        globdat[gn.EGROUPS]['all'] = [*range(len(elems))]

    def _read_geo(self, fname, globdat):
        print('InitModule: Reading geo mesh file', fname, '...')

        nodes = []
        members = []
        nelem = []
        elems = []
        parse_nodes = False
        parse_elems = False

        with open(fname) as msh:
            for line in msh:
                sp = line.split()
                if 'node' in line or 'nodes' in line:
                    parse_nodes = True
                    parse_elems = False

                elif 'member' in line or 'members' in line:
                    parse_nodes = False
                    parse_elems = True

                elif parse_nodes and len(sp) > 1:
                    coords = np.array(sp[1:], dtype=np.float64)
                    nodes.append(Node(coords))

                elif parse_elems and len(sp) > 0:
                    members.append([int(sp[0]), int(sp[1])])
                    nelem.append(int(sp[2]))

        nN = len(nodes)
        inode = nN

        for mem, nel in zip(members, nelem):

            if nel == 1:
                pass

            elif nel > 1:
                x0 = nodes[mem[0]].get_coords()
                x1 = nodes[mem[1]].get_coords()
                dx = (x1 - x0) / nel
                nodes.append(Node(x0 + dx))
                connectivity = np.array([mem[0], inode])
                elems.append(Element(connectivity))  # first element on member

                if nel > 2:
                    for i in range(nel - 2):
                        coords = np.array(x0 + (i + 2) * dx, dtype=np.float64)
                        nodes.append(Node(coords))
                        connectivity = np.array([inode, inode + 1])
                        elems.append(Element(connectivity))  # intermediate elements on member
                        inode += 1

                connectivity = np.array([inode, mem[1]])
                elems.append(Element(connectivity))  # last element on member
                inode += 1

        globdat[gn.NSET] = nodes
        globdat[gn.ESET] = elems
        globdat[gn.NGROUPS]['all'] = [*range(len(nodes))]
        globdat[gn.EGROUPS]['all'] = [*range(len(elems))]

    def _read_meshio(self, mesh, globdat):
        print('Reading mesh from a Meshio object...')

        nodes = []
        elems = []

        for point in mesh.points:
            nodes.append(Node(point))

        if 'triangle' in mesh.cells_dict:
            globdat[gn.MESHSHAPE] = 'Triangle3'
            globdat[gn.MESHRANK] = 2
            for elem in mesh.cells_dict['triangle']:
                elems.append(Element(elem))
        else:
            raise SyntaxError('InitModule: Unsupported Meshio element type')

        globdat[gn.NSET] = nodes
        globdat[gn.ESET] = elems

        globdat[gn.NGROUPS]['all'] = [*range(len(nodes))]
        globdat[gn.EGROUPS]['all'] = [*range(len(elems))]

    def _create_ngroups(self, groups, props, globdat):
        coords = np.stack([node.get_coords() for node in globdat[gn.NSET]], axis=1)
        for g in groups:
            group = np.array(globdat[gn.NGROUPS]['all'])
            gprops = props[g]
            if isinstance(gprops,dict):
                for i, axis in enumerate(['xtype', 'ytype', 'ztype']):
                    if axis in gprops:
                        if gprops[axis] == 'min':
                            ubnd = np.min(coords[i, :]) + self._ctol
                            group = group[coords[i, group] < ubnd]
                        elif gprops[axis] == 'max':
                            lbnd = np.max(coords[i, :]) - self._ctol
                            group = group[coords[i, group] > lbnd]
                        elif gprops[axis] == 'mid':
                            mid = 0.5 * (np.max(coords[i, :]) - np.min(coords[i, :]))
                            lbnd = mid - self._ctol
                            ubnd = mid + self._ctol
                            group = group[coords[i, group] > lbnd]
                            group = group[coords[i, group] < ubnd]
                        else:
                            pass
            else:
                group = pu.parse_list(gprops,int)
                print(gprops, group)

            globdat[gn.NGROUPS][g] = group
            print('InitModule: Created group', g, 'with nodes', group)

    def _create_egroups(self, groups, globdat):
        pass


def declare(factory):
    factory.declare_module('Init', InitModule)
