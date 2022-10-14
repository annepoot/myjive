from jive.fem.names import GlobNames as gn
from jive.app.initmodule import InitModule
from jive.fem.dofspace import DofSpace
import jive.util.proputils as pu

COARSEMESH = 'coarseMesh'
TYPE = 'type'
FILE = 'file'

class GPInitModule(InitModule):

    def init(self, props, globdat):
        myprops = props.get(self._name)

        if not myprops:
            raise KeyError('Properties for GPInitModule not found')

        # Create a separate dictionary to store the coarse mesh properties
        globdat[gn.COARSEMESH] = {}

        # Initialize the node/elemenet group dictionaries
        globdat[gn.COARSEMESH][gn.NGROUPS] = {}
        globdat[gn.COARSEMESH][gn.EGROUPS] = {}

        # Initialize DofSpace
        print('GPInitModule: Creating DofSpace...')
        globdat[gn.COARSEMESH][gn.DOFSPACE] = DofSpace()

        # Read mesh
        meshprops = myprops[COARSEMESH]

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

        # Create node groups in the coarse mesh
        if gn.NGROUPS in myprops:
            print('GPInitModule: Creating node groups...')
            groups = pu.parse_list(myprops[gn.NGROUPS])
            self._create_ngroups(groups, myprops, globdat[gn.COARSEMESH])

        # Create element groups in the coarse mesh
        if gn.EGROUPS in myprops:
            print('GPInitModule: Creating element groups...')
            groups = pu.parse_list(myprops[gn.EGROUPS])
            self._create_egroups(groups, globdat[gn.COARSEMESH])

        # Initialize initmodule
        # Note that this needs to happen last, otherwise gpmodel.configure() will not have the coarse mesh available
        super().init(props, globdat)

    def run(self, globdat):
        return 'ok'

    def shutdown(self, globdat):
        pass


def declare(factory):
    factory.declare_module('GPInit', GPInitModule)
