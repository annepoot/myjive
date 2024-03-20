import numpy as np

from myjive.names import GlobNames as gn
from myjive.model.model import Model
from ..materials import new_material
import myjive.util.proputils as pu
from myjive.util import Table, to_xtable
from myjive.util.proputils import mandatory_argument, mandatory_dict, optional_argument
from .jit.solidmodel import get_N_matrix_jit, get_B_matrix_jit

TYPE = "type"
INTSCHEME = "intScheme"
DOFTYPES = ["dx", "dy", "dz"]

__all__ = ["SolidModel"]


class SolidModel(Model):
    def GETMATRIX0(self, K, globdat, **kwargs):
        K = self._get_matrix(K, globdat, **kwargs)
        return K

    def GETMATRIX2(self, M, globdat, **kwargs):
        M = self._get_mass_matrix(M, globdat, **kwargs)
        return M

    def GETMATRIXB(self, B, wts, globdat, **kwargs):
        B, wts = self._get_strain_matrix(B, wts, globdat, **kwargs)
        return B, wts

    def GETTABLE(self, name, table, tbwts, globdat, **kwargs):
        if "stress" in name:
            table, tbwts = self._get_stresses(table, tbwts, globdat, **kwargs)
        elif "strain" in name:
            table, tbwts = self._get_strains(table, tbwts, globdat, **kwargs)
        elif "stiffness" in name:
            table, tbwts = self._get_stiffness(table, tbwts, globdat, **kwargs)
        elif "size" in name:
            table, tbwts = self._get_elem_size(table, tbwts, globdat, **kwargs)
        return table, tbwts

    def configure(self, globdat, **props):

        # Get props
        shapeprops = mandatory_dict(
            self, props, "shape", mandatory_keys=[TYPE, INTSCHEME]
        )
        elements = mandatory_argument(self, props, "elements")
        matprops = mandatory_dict(self, props, "material", mandatory_keys=[TYPE])
        self._thickness = optional_argument(self, props, "thickness", default=1.0)

        # Configure the material
        mattype = matprops[TYPE]
        self._mat = new_material(mattype)
        self._mat.configure(globdat, **matprops)

        # Get shape and element info
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(
            shapeprops[TYPE], shapeprops[INTSCHEME]
        )

        egroup = globdat[gn.EGROUPS][elements]
        self._elems = egroup.get_elements()
        self._ielems = egroup.get_indices()
        self._nodes = self._elems.get_nodes()

        # Make sure the shape rank and mesh rank are identitcal
        if self._shape.global_rank() != globdat[gn.MESHRANK]:
            raise RuntimeError("ElasticModel: Shape rank must agree with mesh rank")

        # Get basic dimensionality info
        self._rank = self._shape.global_rank()
        self._ipcount = self._shape.ipoint_count()
        self._dofcount = self._rank * self._shape.node_count()
        self._strcount = self._rank * (self._rank + 1) // 2  # 1-->1, 2-->3, 3-->6

        if self._rank == 2:
            self._thickness = pu.soft_cast(self._thickness, float)

        # Create a new dof for every node and dof type
        for doftype in DOFTYPES[0 : self._rank]:
            globdat[gn.DOFSPACE].add_type(doftype)
            for inode in self._elems.get_unique_nodes_of(self._ielems):
                globdat[gn.DOFSPACE].add_dof(inode, doftype)

    def _get_matrix(self, K, globdat, unit_matrix=False):
        if K is None:
            dc = globdat[gn.DOFSPACE].dof_count()
            K = np.zeros((dc, dc))

        if unit_matrix:
            D_elem = np.identity(self._rank)

        for ielem in self._ielems:
            # Get the nodal coordinates of each element
            inodes = self._elems.get_elem_nodes(ielem)
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0 : self._rank])
            coords = self._nodes.get_some_coords(inodes)

            # Get the gradients, weights and coordinates of each integration point
            grads, weights = self._shape.get_shape_gradients(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            if self._rank == 2 and not unit_matrix:
                weights *= self._thickness

            # Reset the element stiffness matrix
            elmat = np.zeros((self._dofcount, self._dofcount))

            for ip in range(self._ipcount):
                # Get the B and D matrices for each integration point
                B_elem = self._get_B_matrix(grads[:, :, ip])

                if not unit_matrix:
                    D_elem = self._mat.stiff_at_point(ipcoords[:, ip])

                # Compute the element stiffness matrix
                elmat += weights[ip] * np.matmul(
                    np.transpose(B_elem), np.matmul(D_elem, B_elem)
                )

            # Add the element stiffness matrix to the global stiffness matrix
            K[np.ix_(idofs, idofs)] += elmat

        return K

    def _get_mass_matrix(self, M, globdat, unit_matrix=False):
        if M is None:
            dc = globdat[gn.DOFSPACE].dof_count()
            M = np.zeros((dc, dc))

        if unit_matrix:
            M_elem = np.identity(self._rank)

        for ielem in self._ielems:
            # Get the nodal coordinates of each element
            inodes = self._elems.get_elem_nodes(ielem)
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0 : self._rank])
            coords = self._nodes.get_some_coords(inodes)

            # Get the shape functions, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            weights = self._shape.get_integration_weights(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            if self._rank == 2 and not unit_matrix:
                weights *= self._thickness

            # Reset the element mass matrix
            elmat = np.zeros((self._dofcount, self._dofcount))

            for ip in range(self._ipcount):
                # Get the N and M matrices for each integration point
                N_elem = self._get_N_matrix(sfuncs[:, ip])

                if not unit_matrix:
                    M_elem = self._mat.mass_at_point(ipcoords[:, ip])

                # Compute the element mass matrix
                elmat += weights[ip] * np.matmul(
                    np.transpose(N_elem), np.matmul(M_elem, N_elem)
                )

            # Add the element mass matrix to the global mass matrix
            M[np.ix_(idofs, idofs)] += elmat

        return M

    def _get_strain_matrix(self, B, wts, globdat):
        # Add the element weights to the global weights
        nc = self._nodes.size()
        node_count = self._shape.node_count()

        for ielem in self._ielems:
            # Get the nodal coordinates of each element
            inodes = self._elems.get_elem_nodes(ielem)
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0 : self._rank])
            coords = self._nodes.get_some_coords(inodes)

            # Get the gradients, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords)

            # Reset the element strain matrix
            elbmat = np.zeros((node_count * self._strcount, self._dofcount))
            elwts = np.zeros(node_count * self._strcount)

            for ip in range(self._ipcount):
                # Get the B and D matrices for each integration point
                B_elem = self._get_B_matrix(grads[:, :, ip])

                # Compute the element strain and weights
                for i in range(self._strcount):
                    elbmat[i * node_count : (i + 1) * node_count, :] += np.outer(
                        sfuncs[:, ip], B_elem[i, :]
                    )
                    elwts[i * node_count : (i + 1) * node_count] += sfuncs[
                        :, ip
                    ].flatten()

            # Get the node index vector
            node_idx = np.zeros(node_count * self._strcount, dtype=int)
            for i in range(self._strcount):
                node_idx[i * node_count : (i + 1) * node_count] = inodes + nc * i

            # Add the element strain matrix to the global strain matrix
            B[np.ix_(node_idx, idofs)] += elbmat

            # Add the element weights to the global weights
            wts[node_idx] += elwts

        return B, wts

    def _get_strains(self, table, tbwts, globdat, solution=None):
        if table is None:
            nodecount = len(globdat[gn.NSET])
            table = Table(size=nodecount)

        if tbwts is None:
            nodecount = len(globdat[gn.NSET])
            tbwts = np.zeros(nodecount)

        # Convert the table to an XTable and store the original class
        xtable = to_xtable(table)

        # Get the STATE0 vector if no custom displacement field is provided
        if solution is None:
            disp = globdat[gn.STATE0]
        else:
            disp = solution

        # Add the columns of all strain components to the table
        if self._rank == 1:
            jcols = xtable.add_columns(["xx"])
        elif self._rank == 2:
            jcols = xtable.add_columns(["xx", "yy", "xy"])
        elif self._rank == 3:
            jcols = xtable.add_columns(["xx", "yy", "zz", "xy", "yz", "zx"])

        for ielem in self._ielems:
            # Get the nodal coordinates of each element
            inodes = self._elems.get_elem_nodes(ielem)
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0 : self._rank])
            coords = self._nodes.get_some_coords(inodes)

            # Get the shape functions, gradients, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords)

            # Get the nodal displacements
            eldisp = disp[idofs]

            # Reset the element stress matrix and weights
            eleps = np.zeros((self._shape.node_count(), self._strcount))
            elwts = np.zeros(self._shape.node_count())

            for ip in range(self._ipcount):
                # Get the B matrix for each integration point
                B = self._get_B_matrix(grads[:, :, ip])

                # Get the strain of the element in the integration point
                strain = np.matmul(B, eldisp)

                # Compute the element strain and weights
                eleps += np.outer(sfuncs[:, ip], strain)
                elwts += sfuncs[:, ip].flatten()

            # Add the element weights to the global weights
            tbwts[inodes] += elwts

            # Add the element stresses to the global stresses
            xtable.add_block(inodes, jcols, eleps)

        # Convert the table back to the original class
        table = xtable.to_table()

        return table, tbwts

    def _get_stresses(self, table, tbwts, globdat, solution=None):
        if table is None:
            nodecount = len(globdat[gn.NSET])
            table = Table(size=nodecount)

        if tbwts is None:
            nodecount = len(globdat[gn.NSET])
            tbwts = np.zeros(nodecount)

        # Convert the table to an XTable and store the original class
        xtable = to_xtable(table)

        # Get the STATE0 vector if no custom displacement field is provided
        if solution is None:
            disp = globdat[gn.STATE0]
        else:
            disp = solution

        # Add the columns of all stress components to the table
        if self._rank == 1:
            jcols = xtable.add_columns(["xx"])
        elif self._rank == 2:
            jcols = xtable.add_columns(["xx", "yy", "xy"])
        elif self._rank == 3:
            jcols = xtable.add_columns(["xx", "yy", "zz", "xy", "yz", "zx"])

        for ielem in self._ielems:
            # Get the nodal coordinates of each element
            inodes = self._elems.get_elem_nodes(ielem)
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0 : self._rank])
            coords = self._nodes.get_some_coords(inodes)

            # Get the shape functions, gradients, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            if self._rank == 2:
                weights *= self._thickness

            # Get the nodal displacements
            eldisp = disp[idofs]

            # Reset the element stress matrix and weights
            elsig = np.zeros((self._shape.node_count(), self._strcount))
            elwts = np.zeros(self._shape.node_count())

            for ip in range(self._ipcount):
                # Get the B matrix for each integration point
                B = self._get_B_matrix(grads[:, :, ip])

                # Get the strain and stress of the element in the integration point
                strain = np.matmul(B, eldisp)
                stress = self._mat.stress_at_point(strain, ipcoords[:, ip])

                # Compute the element stress and weights
                elsig += np.outer(sfuncs[:, ip], stress)
                elwts += sfuncs[:, ip].flatten()

            # Add the element weights to the global weights
            tbwts[inodes] += elwts

            # Add the element stresses to the global stresses
            xtable.add_block(inodes, jcols, elsig)

        # Convert the table back to the original class
        table = xtable.to_table()

        return table, tbwts

    def _get_stiffness(self, table, tbwts, globdat):
        if table is None:
            nodecount = len(globdat[gn.NSET])
            table = Table(size=nodecount)

        if tbwts is None:
            nodecount = len(globdat[gn.NSET])
            tbwts = np.zeros(nodecount)

        # Convert the table to an XTable and store the original class
        xtable = to_xtable(table)

        # Add the column of the Young's modulus to the table
        jcol = xtable.add_column("")

        for ielem in self._ielems:
            # Get the nodal coordinates of each element
            inodes = self._elems.get_elem_nodes(ielem)
            coords = self._nodes.get_some_coords(inodes)

            # Get the shape functions, gradients, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            ipcoords = self._shape.get_global_integration_points(coords)

            # Reset the element stress matrix and weights
            elyoung = np.zeros((self._shape.node_count()))
            elwts = np.zeros(self._shape.node_count())

            for ip in range(self._ipcount):
                # Get the stiffness in the integration point
                E = self._mat._get_E(ipcoords[:, ip])

                # Compute the element stiffness and weights
                elyoung += E * sfuncs[:, ip]
                elwts += sfuncs[:, ip].flatten()

            # Add the element weights to the global weights
            tbwts[inodes] += elwts

            # Add the element stiffness to the global stiffness
            xtable.add_col_values(inodes, jcol, elyoung)

        # Convert the table back to the original class
        table = xtable.to_table()

        return table, tbwts

    def _get_elem_size(self, table, tbwts, globdat):
        if table is None:
            nodecount = len(globdat[gn.NSET])
            table = Table(size=nodecount)

        if tbwts is None:
            nodecount = len(globdat[gn.NSET])
            tbwts = np.zeros(nodecount)

        # Convert the table to an XTable and store the original class
        xtable = to_xtable(table)

        # Add the columns of all stress components to the table
        jcol = xtable.add_column("")

        for ielem in self._ielems:
            # Get the nodal coordinates of each element
            inodes = self._elems.get_elem_nodes(ielem)
            coords = self._nodes.get_some_coords(inodes)

            # Compute the maximum edge length\
            max_edge = 0.0
            for i, inode in enumerate(inodes):
                for j, jnode in enumerate(inodes):
                    icoords = coords[:, i]
                    jcoords = coords[:, j]
                    edge = np.sqrt(np.sum((icoords - jcoords) ** 2))
                    if edge > max_edge:
                        max_edge = edge

            # Reset the element size and weights
            elsize = max_edge * np.ones(self._shape.node_count())
            elwts = np.ones(self._shape.node_count())

            # Add the element weights to the global weights
            tbwts[inodes] += elwts

            # Add the element stresses to the global stresses
            xtable.add_col_values(inodes, jcol, elsize)

        # Convert the table back to the original class
        table = xtable.to_table()

        return table, tbwts

    def _get_N_matrix(self, sfuncs):
        return get_N_matrix_jit(sfuncs, self._dofcount, self._rank)

    def _get_B_matrix(self, grads):
        return get_B_matrix_jit(
            grads, self._strcount, self._dofcount, self._shape.node_count(), self._rank
        )
