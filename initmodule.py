import numpy as np

from model    import *
from module   import *
from node     import Node
from element  import Element
from dofspace import DofSpace

from names    import GlobNames as gn
from names    import PropNames as pn

MESH    = 'mesh'

class InitModule (Module):
    def init (self, props, globdat):
        myprops = props.get (self._name)

        if not myprops:
            raise KeyError ('Properties for InitModule not found')

        # Initialize some parameters

        self._ctol = 1.e-5
        globdat[gn.NGROUPS] = {}
        globdat[gn.EGROUPS] = {}
        modelfac   = globdat[gn.MODELFACTORY]
        modelprops = props[gn.MODEL]

        # Initialize DofSpace
        print('InitModule: Creating DofSpace...')
        globdat[gn.DOFSPACE] = DofSpace() 

        # Read mesh
        print('InitModule: Reading mesh file',myprops[MESH],'...')
        self._read_gmsh (myprops[MESH], globdat)

        # Create node groups
        if gn.NGROUPS in myprops:
            print('InitModule: Creating node groups...')
            groups = myprops[gn.NGROUPS].strip('[').strip(']').split(',')
            self._create_ngroups (groups,myprops,globdat)

        # Create element groups
        if gn.EGROUPS in myprops:
            print('InitModule: Creating element groups...')
            groups = myprops[gn.EGROUPS].strip('[').strip(']').split(',')
            self._create_egroups (groups,globdat)

        # Initialize model
        print('InitModule: Creating model...')
        m = modelfac.get_model (modelprops[pn.TYPE],gn.MODEL)
        m.configure (modelprops,globdat)
        globdat[gn.MODEL] = m

    def run (self, globdat):
        return('ok')

    def shutdown (self, globdat):
        pass

    def _read_gmsh (self, fname, globdat):
        if not fname.endswith('.msh'):
            raise RuntimeError ('Unexpected mesh file extension')

        nodes = []
        elems = []

        parse_nodes = False
        parse_elems = False
        eltype      = 0
        nnodes      = 0

        with open (fname) as msh:
            for line in msh:
                sp = line.split()

                if '$Nodes' in line:
                    parse_nodes = True

                elif '$Elements' in line:
                    parse_nodes = False
                    parse_elems = True

                elif parse_nodes and len(sp) > 1:
                    if len(sp[1:]) != 3:
                        raise SyntaxError ('InitModule: Three coordinates per node are expected')
                    coords = np.array (sp[1:],dtype=np.float64)
                    nodes.append (Node (coords))

                elif parse_elems and len(sp) > 1:
                    if eltype == 0:
                        eltype = int(sp[1])
                        if eltype == 1:
                            globdat[gn.MESHSHAPE] = 'Line2'
                            globdat[gn.MESHRANK]  = 1
                            nnodes = 2
                        elif eltype == 2:
                            globdat[gn.MESHSHAPE] = 'Triangle3'
                            globdat[gn.MESHRANK]  = 2
                            nnodes = 3
                        elif eltype == 3:
                            globdat[gn.MESHSHAPE] = 'Quad4'
                            globdat[gn.MESHRANK]  = 2
                            nnodes = 4
                        elif eltype == 4:
                            globdat[gn.MESHSHAPE] = 'Tet4'
                            globdat[gn.MESHRANK]  = 3
                            nnodes = 4
                        elif eltype == 5:
                            globdat[gn.MESHSHAPE] = 'Brick8'
                            globdat[gn.MESHRANK]  = 3
                            nnodes = 8
                        else:
                            raise SyntaxError ('InitModule: Unsupported element type')
                    elif eltype != int(sp[1]):
                        raise SyntaxError ('InitModule: Only one element type per mesh is supported')
                    inodes = np.array ( sp[3+int(sp[2]):], dtype=np.int32 ) - 1
                    if len(inodes) != nnodes:
                        raise SyntaxError ('InitModule: Could not read element with incorrect number of nodes')
                    elems.append (Element (inodes))

        globdat[gn.NSET] = nodes
        globdat[gn.ESET] = elems

        globdat[gn.NGROUPS]['all'] = [*range(len(nodes))]
        globdat[gn.EGROUPS]['all'] = [*range(len(elems))]

    def _create_ngroups (self,groups,props,globdat):
        coords = np.stack([node.get_coords() for node in globdat[gn.NSET]],axis=1)
        for g in groups:
            group = np.array(globdat[gn.NGROUPS]['all'])
            gprops = props[g]

            for i, axis in enumerate(['xtype', 'ytype', 'ztype']):
                if axis in gprops:
                    if gprops[axis] == 'min':
                        ubnd = np.min(coords[i,:]) + self._ctol
                        group = group[coords[i,group] < ubnd]
                    elif gprops[axis] == 'max':
                        lbnd = np.max(coords[i,:]) - self._ctol
                        group = group[coords[i,group] > lbnd]
                    elif gprops[axis] == 'mid':
                        mid = 0.5 * (np.max(coords[i,:]) - np.min(coords[i,:]))
                        lbnd = mid - self._ctol
                        ubnd = mid + self._ctol
                        group = group[coords[i,group] > lbnd]
                        group = group[coords[i,group] < ubnd]
                    else:
                        pass

            globdat[gn.NGROUPS][g] = group
            print('InitModule: Created group',g,'with nodes',group)


    def _create_egroups (self,groups,globdat):
        pass

def declare (factory):
    factory.declare_module ('Init', InitModule)
