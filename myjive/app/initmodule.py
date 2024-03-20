import numpy as np

from .module import Module
from ..fem import XNodeSet
from ..fem import NodeGroup
from ..fem import XElementSet
from ..fem import ElementGroup
from ..fem import DofSpace

from ..names import GlobNames as gn

from ..util.proputils import optional_argument, mandatory_dict

TYPE = "type"
FILE = "file"
MODELS = "models"

__all__ = ["InitModule"]


class InitModule(Module):
    # Predefine some parameters
    _ctol = 1.0e-5

    def init(self, globdat, **props):

        # Get props
        modelprops = mandatory_dict(self, props, "modelprops", mandatory_keys=[MODELS])
        meshprops = mandatory_dict(self, props, "mesh", mandatory_keys=[TYPE, FILE])
        nodeGroups = optional_argument(self, props, "nodeGroups", default=[])
        elemGroups = optional_argument(self, props, "elemGroups", default=[])

        # Initialize the node/elemenet group dictionaries
        globdat[gn.NGROUPS] = {}
        globdat[gn.EGROUPS] = {}

        modelfac = globdat[gn.MODELFACTORY]

        # Initialize DofSpace
        print("InitModule: Creating DofSpace...")
        globdat[gn.DOFSPACE] = DofSpace()

        if "gmsh" in meshprops[TYPE]:
            self._read_gmsh(meshprops[FILE], globdat)
        elif "manual" in meshprops[TYPE]:
            self._read_mesh(meshprops[FILE], globdat)
        elif "meshio" in meshprops[TYPE]:
            self._read_meshio(meshprops[FILE], globdat)
        elif "geo" in meshprops[TYPE]:
            self._read_geo(meshprops[FILE], globdat)
        else:
            raise KeyError("InitModule: Mesh input type unknown")

        # Create node groups
        if nodeGroups is not None:
            print("InitModule: Creating node groups...")
            self._create_ngroups(nodeGroups, globdat, **props)

        # Create element groups
        if elemGroups is not None:
            print("InitModule: Creating element groups...")
            self._create_egroups(elemGroups, globdat)

        # Initialize model
        print("InitModule: Creating models...")
        name_list = modelprops[MODELS]
        model_list = []
        for name in name_list:
            m = modelfac.get_model(modelprops[name][TYPE], name)
            m.configure(globdat, **(modelprops[name]))
            model_list.append(m)
        globdat[gn.MODELS] = model_list

    def run(self, globdat):
        return "ok"

    def shutdown(self, globdat):
        pass

    def _read_gmsh(self, fname, globdat):
        print("InitModule: Reading mesh file", fname, "...")

        if not fname.endswith(".msh"):
            raise RuntimeError("Unexpected mesh file extension")

        # Create a clean nodeset and elementset
        nodes = XNodeSet()
        elems = XElementSet(nodes)

        with open(fname) as msh:
            lines = msh.readlines()

            # Extract the Nodes and Elements blocks from the gmsh file
            nlines = lines[lines.index("$Nodes\n") + 2 : lines.index("$EndNodes\n")]
            elines = lines[
                lines.index("$Elements\n") + 2 : lines.index("$EndElements\n")
            ]

            # If possible, extract the element group block from the gmsh file
            if "$PhysicalNames\n" in lines:
                grouplines = lines[
                    lines.index("$PhysicalNames\n")
                    + 2 : lines.index("$EndPhysicalNames\n")
                ]

                # Split the element group info
                group_ids = np.genfromtxt(grouplines, dtype=int, ndmin=2)[:, 1]
                group_names = np.genfromtxt(grouplines, dtype=str, ndmin=2)[:, 2]
                group_names = np.char.strip(group_names, '"')

            else:
                group_ids = []
                group_names = []

            # Split the node info
            node_ids = np.genfromtxt(nlines, dtype=int, ndmin=2)[:, 0]
            coords = np.genfromtxt(nlines, dtype=float, ndmin=2)[:, 1:]

            # Split the element info
            elem_ids = np.genfromtxt(elines, dtype=int, ndmin=2)[:, 0]
            elem_info = np.genfromtxt(elines, dtype=int, ndmin=2)[:, 1:5]
            inodes = np.genfromtxt(elines, dtype=int, ndmin=2)[:, 5:]

            # Get the element type, and make sure it is the only one
            eltype = elem_info[0, 0]
            nnodes = 0
            if not np.all(elem_info[:, 0] == eltype):
                raise SyntaxError(
                    "InitModule: Only one element type per mesh is supported"
                )

            # Get the info belonging to the element type
            if eltype == 1:
                globdat[gn.MESHSHAPE] = "Line2"
                globdat[gn.MESHRANK] = 1
                nnodes = 2
            elif eltype == 2:
                globdat[gn.MESHSHAPE] = "Triangle3"
                globdat[gn.MESHRANK] = 2
                nnodes = 3
            elif eltype == 3:
                globdat[gn.MESHSHAPE] = "Quad4"
                globdat[gn.MESHRANK] = 2
                nnodes = 4
            elif eltype == 4:
                globdat[gn.MESHSHAPE] = "Tet4"
                globdat[gn.MESHRANK] = 3
                nnodes = 4
            elif eltype == 5:
                globdat[gn.MESHSHAPE] = "Brick8"
                globdat[gn.MESHRANK] = 3
                nnodes = 8
            elif eltype == 8:
                globdat[gn.MESHSHAPE] = "Line3"
                globdat[gn.MESHRANK] = 1
                nnodes = 3
            elif eltype == 9:
                globdat[gn.MESHSHAPE] = "Triangle6"
                globdat[gn.MESHRANK] = 2
                nnodes = 6
            elif eltype == 10:
                globdat[gn.MESHSHAPE] = "Quad9"
                globdat[gn.MESHRANK] = 2
                nnodes = 9
            else:
                raise SyntaxError("InitModule: Unsupported element type")

            # Make sure that the correct number of nodes and coordinates is passed
            if coords.shape[1] != 3:
                raise SyntaxError("InitModule: Three coordinates per node are expected")
            if inodes.shape[1] != nnodes:
                raise SyntaxError(
                    "InitModule: Could not read element with incorrect number of nodes"
                )

            # Add all nodes to the node set
            for i in range(coords.shape[0]):
                nodes.add_node(coords[i, : globdat[gn.MESHRANK]], node_id=node_ids[i])

            # Add all elements to the element set
            for i in range(inodes.shape[0]):
                elems.add_element(nodes.find_nodes(inodes[i, :]), elem_id=elem_ids[i])

            for j, (idx, name) in enumerate(zip(group_ids, group_names)):
                egroup = []
                for i in range(elem_info.shape[0]):
                    if elem_info[i, 2] == idx:
                        egroup.append(i)
                globdat[gn.EGROUPS][name] = ElementGroup(elems, egroup)

        # Convert the XNodeSet and XElementSet to a normal NodeSet and ElementSet
        globdat[gn.NSET] = nodes.to_nodeset()
        globdat[gn.ESET] = elems.to_elementset()

        # Create node and element groups containing all items
        globdat[gn.NGROUPS]["all"] = NodeGroup(nodes, [*range(nodes.size())])
        globdat[gn.EGROUPS]["all"] = ElementGroup(elems, [*range(elems.size())])

    def _read_mesh(self, fname, globdat):
        print("InitModule: Reading manual mesh file", fname, "...")

        # Create a clean nodeset and elementset
        nodes = XNodeSet()
        elems = XElementSet(nodes)

        parse_nodes = False
        parse_elems = False

        with open(fname) as msh:
            for line in msh:
                sp = line.split()

                if "nodes" in line:
                    parse_nodes = True
                    parse_elems = False

                elif "elements" in line or "elems" in line:
                    parse_nodes = False
                    parse_elems = True

                elif parse_nodes and len(sp) > 1:
                    nodes.add_node(sp[1:], sp[0])
                    globdat[gn.MESHRANK] = len(sp) - 1

                elif parse_elems and len(sp) > 0:
                    elems.add_element(sp)

        # Convert the XNodeSet and XElementSet to a normal NodeSet and ElementSet
        globdat[gn.NSET] = nodes.to_nodeset()
        globdat[gn.ESET] = elems.to_elementset()

        # Create node and element groups containing all items
        globdat[gn.NGROUPS]["all"] = NodeGroup(nodes, [*range(nodes.size())])
        globdat[gn.EGROUPS]["all"] = ElementGroup(elems, [*range(elems.size())])

    def _read_geo(self, fname, globdat):
        print("InitModule: Reading geo mesh file", fname, "...")

        # Create a clean nodeset and elementset
        nodes = XNodeSet()
        elems = XElementSet(nodes)

        parse_nodes = False
        parse_elems = False

        members = []
        nelem = []

        with open(fname) as msh:
            for line in msh:
                sp = line.split()
                if "node" in line or "nodes" in line:
                    parse_nodes = True
                    parse_elems = False

                elif "member" in line or "members" in line:
                    parse_nodes = False
                    parse_elems = True

                elif parse_nodes and len(sp) > 1:
                    nodes.add_node(sp[1:], sp[0])

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
                elems.add_element(mem)

            else:
                x0 = nodes[mem[0]].get_coords()
                x1 = nodes[mem[1]].get_coords()
                dx = (x1 - x0) / nel
                nodes.add_node(x0 + dx)
                connectivity = np.array([mem[0], inode])
                elems.add_element(connectivity)  # first element on member

                if nel > 2:
                    for i in range(nel - 2):
                        coords = np.array(x0 + (i + 2) * dx, dtype=np.float64)
                        nodes.add_node(coords)
                        connectivity = np.array([inode, inode + 1])
                        elems.add_element(
                            connectivity
                        )  # intermediate elements on member
                        inode += 1

                connectivity = np.array([inode, mem[1]])
                elems.add_element(connectivity)  # last element on member
                inode += 1

            ie1 = len(elems)
            globdat[gn.EGROUPS]["member" + str(imember)] = ElementGroup(
                elems, [*range(ie0, ie1)]
            )
            imember += 1

        print("done reading geo " + str(len(elems)) + " elements")

        # Convert the XNodeSet and XElementSet to a normal NodeSet and ElementSet
        globdat[gn.NSET] = nodes.to_nodeset()
        globdat[gn.ESET] = elems.to_elementset()

        # Create node and element groups containing all items
        globdat[gn.NGROUPS]["all"] = NodeGroup(nodes, [*range(nodes.size())])
        globdat[gn.EGROUPS]["all"] = ElementGroup(elems, [*range(elems.size())])

    def _read_meshio(self, mesh, globdat):
        print("Reading mesh from a Meshio object...")

        # Create a clean nodeset and elementset
        nodes = XNodeSet()
        elems = XElementSet(nodes)

        for point in mesh.points:
            nodes.add_node(point)

        if "triangle" in mesh.cells_dict:
            globdat[gn.MESHSHAPE] = "Triangle3"
            globdat[gn.MESHRANK] = 2
            for elem in mesh.cells_dict["triangle"]:
                elems.add_element(elem)
        else:
            raise SyntaxError("InitModule: Unsupported Meshio element type")

        # Convert the XNodeSet and XElementSet to a normal NodeSet and ElementSet
        globdat[gn.NSET] = nodes.to_nodeset()
        globdat[gn.ESET] = elems.to_elementset()

        # Create node and element groups containing all items
        globdat[gn.NGROUPS]["all"] = NodeGroup(nodes, [*range(nodes.size())])
        globdat[gn.EGROUPS]["all"] = ElementGroup(elems, [*range(elems.size())])

    def _create_ngroups(self, groups, globdat, **props):
        coords = globdat[gn.NSET].get_coords()
        cmax = np.max(coords, axis=1)
        cmin = np.min(coords, axis=1)
        cmid = 0.5 * (cmax + cmin)
        for g in groups:
            group = globdat[gn.NGROUPS]["all"].get_indices()
            gprops = props[g]
            if isinstance(gprops, dict):
                for i, axis in enumerate(["xtype", "ytype", "ztype"]):
                    if axis in gprops:
                        loc = str(gprops[axis])
                        if loc.replace(".", "").isnumeric():
                            lbnd = float(loc) - self._ctol
                            ubnd = float(loc) + self._ctol
                            group = group[coords[i, group] > lbnd]
                            group = group[coords[i, group] < ubnd]
                        elif gprops[axis] == "min":
                            ubnd = cmin[i] + self._ctol
                            group = group[coords[i, group] < ubnd]
                        elif gprops[axis] == "max":
                            lbnd = cmax[i] - self._ctol
                            group = group[coords[i, group] > lbnd]
                        elif gprops[axis] == "mid":
                            lbnd = cmid[i] - self._ctol
                            ubnd = cmid[i] + self._ctol
                            group = group[coords[i, group] > lbnd]
                            group = group[coords[i, group] < ubnd]
                        else:
                            pass
            else:
                group = gprops

            globdat[gn.NGROUPS][g] = NodeGroup(globdat[gn.NSET], group)

            print("InitModule: Created group", g, "with nodes", group)

    def _create_egroups(self, groups, globdat):
        pass
