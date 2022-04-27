from names import GlobNames as gn
from names import PropNames as prn
from initmodule import InitModule
from dofspace import DofSpace
import proputils as pu

MESH = 'mesh'
TYPE = 'type'
FILE = 'file'

class GPInitModule(InitModule):

    def init(self, props, globdat):
        myprops = props.get(self._name)

        if not myprops:
            raise KeyError('Properties for GPInitModule not found')

        # Initialize some parameters
        self._ctol = 1.e-5
        modelfac = globdat[gn.MODELFACTORY]

        # Get the appropriate model for this module
        self._modelname = myprops.get(gn.MODEL, gn.MODEL)
        modelprops = props[self._modelname]

        # Create a separate dictionary to store the coarse mesh properties
        globdat[gn.COARSEMESH] = {}

        # Initialize the node/elemenet group dictionaries
        globdat[gn.COARSEMESH][gn.NGROUPS] = {}
        globdat[gn.COARSEMESH][gn.EGROUPS] = {}

        # Initialize DofSpace
        print('GPInitModule: Creating DofSpace...')
        globdat[gn.COARSEMESH][gn.DOFSPACE] = DofSpace()

        # Read mesh
        meshprops = myprops[MESH]

        if 'gmsh' in meshprops[TYPE]:
            self._read_gmsh(meshprops[FILE], globdat[gn.COARSEMESH])
        elif 'manual' in meshprops[TYPE]:
            self._read_mesh(meshprops[FILE], globdat[gn.COARSEMESH])
        elif 'meshio' in meshprops[TYPE]:
            self._read_meshio(meshprops[FILE], globdat[gn.COARSEMESH])
        elif 'geo' in meshprops[TYPE]:
            self._read_geo(meshprops[FILE], globdat[gn.COARSEMESH])
        else:
            raise KeyError('GPInitModule: Mesh input type unknown')

        # Create node groups
        if gn.NGROUPS in myprops:
            print('GPInitModule: Creating node groups...')
            groups = pu.parse_list(myprops[gn.NGROUPS])
            self._create_ngroups(groups, myprops, globdat[gn.COARSEMESH])

        # Create element groups
        if gn.EGROUPS in myprops:
            print('GPInitModule: Creating element groups...')
            groups = pu.parse_list(myprops[gn.EGROUPS])
            self._create_egroups(groups, globdat[gn.COARSEMESH])

        # Initialize model
        print('GPInitModule: Creating model...')
        m = modelfac.get_model(modelprops[prn.TYPE], self._modelname)
        m.configure(modelprops, globdat)
        globdat[self._modelname] = m

    def run(self, globdat):
        model = globdat[self._modelname]

        # Configure the model again, to make sure K, M and f are stored there as well
        model.configure_fem(globdat)

        return 'ok'

    def shutdown(self, globdat):
        pass


def declare(factory):
    factory.declare_module('GPInit', GPInitModule)
