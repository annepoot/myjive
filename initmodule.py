import numpy as np

from module import Module
from node import Node
from nodeset import XNodeSet
from nodegroup import NodeGroup
from element import Element
from elementset import XElementSet
from elementgroup import ElementGroup
from dofspace import DofSpace

from names import GlobNames as gn
from names import PropNames as prn

import proputils as pu

MESH = 'mesh'
TYPE = 'type'
FILE = 'file'

class InitModule(Module):

    # Predefine some parameters
    _ctol = 1.e-5

    def init(self, props, globdat):
        myprops = props.get(self._name)

        if not myprops:
            raise KeyError('Properties for InitModule not found')

        # Initialize the node/elemenet group dictionaries
        globdat[gn.NGROUPS] = {}
        globdat[gn.EGROUPS] = {}

        modelfac = globdat[gn.MODELFACTORY]
        modelprops = props[gn.MODEL]

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
        m = modelfac.get_model(modelprops[prn.TYPE], gn.MODEL)
        m.configure(modelprops, globdat)
        globdat[gn.MODEL] = m

    def run(self, globdat):
        return 'ok'

    def shutdown(self, globdat):
        pass

    def _read_gmsh(self, fname, globdat):
        print('InitModule: Reading mesh file', fname, '...')

        if not fname.endswith('.msh'):
            raise RuntimeError('Unexpected mesh file extension')

        nodes = XNodeSet()
        elems = XElementSet(nodes)

        parse_nodes = False
        parse_elems = False
        eltype = 0
        nnodes = 0

        with open(fname) as msh:
            lines = msh.readlines()

            # Extract the Nodes and Elements blocks from the gmsh file
            nlines = lines[lines.index('$Nodes\n')+2:lines.index('$EndNodes\n')]
            elines = lines[lines.index('$Elements\n')+2:lines.index('$EndElements\n')]

            # Split the node and element info
            node_ids = np.genfromtxt(nlines, dtype=int)[:,0]
            coords = np.genfromtxt(nlines, dtype=float)[:,1:]

            elem_ids = np.genfromtxt(elines, dtype=int)[:,0]
            elem_info = np.genfromtxt(elines, dtype=int)[:,1:5]
            inodes = np.genfromtxt(elines, dtype=int)[:,5:]

            # Get the element type, and make sure it is the only one
            eltype = elem_info[0,0]
            if not np.all(elem_info[:,0] == eltype):
                raise SyntaxError('InitModule: Only one element type per mesh is supported')

            # Get the info belonging to the element type
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
            elif eltype == 9:
                globdat[gn.MESHSHAPE] = 'Triangle6'
                globdat[gn.MESHRANK] = 2
                nnodes = 6
            else:
                raise SyntaxError('InitModule: Unsupported element type')

            # Make sure that the correct number of nodes and coordinates is passed
            if coords.shape[1] != 3:
                raise SyntaxError('InitModule: Three coordinates per node are expected')
            if inodes.shape[1] != nnodes:
                raise SyntaxError('InitModule: Could not read element with incorrect number of nodes')

            # Add all nodes to the node set
            for i in range(coords.shape[0]):
                nodes.add_node(coords[i,:globdat[gn.MESHRANK]], node_id=node_ids[i])

            # Add all elements to the element set
            for i in range(inodes.shape[0]):
                elems.add_element(nodes.find_nodes(inodes[i,:]), elem_id=elem_ids[i])

        # Convert the XNodeSet and XElementSet to a normal NodeSet and ElementSet
        globdat[gn.NSET] = nodes.to_nodeset()
        globdat[gn.ESET] = elems.to_elementset()

        # Create node and element groups containing all items
        globdat[gn.NGROUPS]['all'] = NodeGroup(nodes, [*range(nodes.size())])
        globdat[gn.EGROUPS]['all'] = ElementGroup(elems, [*range(elems.size())])

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
                    rank = len(sp) - 1

                elif parse_elems and len(sp) > 0:
                    connectivity = np.array(sp, dtype=np.int16)
                    elems.append(Element(connectivity))

        globdat[gn.NSET] = nodes
        globdat[gn.ESET] = elems
        globdat[gn.MESHRANK] = rank

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
        imember = 0

        for mem, nel in zip(members, nelem):

            ie0 = len(elems)

            if nel <= 0:
                pass

            elif nel == 1:
                elems.append(Element(mem))

            else:
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

            ie1 = len(elems)
            globdat[gn.EGROUPS]['member'+str(imember)] = [*range(ie0,ie1)]
            imember += 1

        print('done reading geo ' + str(len(elems)) + ' elements')
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
        coords = globdat[gn.NSET].get_coords()
        for g in groups:
            group = globdat[gn.NGROUPS]['all'].get_indices()
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

            globdat[gn.NGROUPS][g] = NodeGroup(globdat[gn.NSET], group)

            print('InitModule: Created group', g, 'with nodes', group)

    def _create_egroups(self, groups, globdat):
        pass


def declare(factory):
    factory.declare_module('Init', InitModule)
