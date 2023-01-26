import numpy as np

from jive.fem.names import Actions as act
from jive.fem.names import ParamNames as pn
from jive.fem.names import GlobNames as gn
from jive.model.model import Model
from material import new_material
from jive.util.xtable import XTable
import jive.util.proputils as pu

ELEMENTS = 'elements'
RHO = 'rho'
SHAPE = 'shape'
TYPE = 'type'
INTSCHEME = 'intScheme'
DOFTYPES = ['dx', 'dy', 'dz']
MATERIAL = 'material'
THICKNESS_PROP = 'thickness'


class SolidModel(Model):
    def take_action(self, action, params, globdat):
        showmsg = True

        if action == act.GETMATRIX0:
            self._get_matrix(params, globdat)
        elif action == act.GETMATRIX2:
            self._get_mass_matrix(params, globdat)
        elif action == act.GETTABLE:
            if 'stress' in params[pn.TABLENAME]:
                self._get_stresses(params, globdat)
            elif 'strain' in params[pn.TABLENAME]:
                self._get_strains(params, globdat)
            elif 'stiffness' in params[pn.TABLENAME]:
                self._get_stiffness(params, globdat)
            else:
                showmsg = False
        else:
            showmsg = False

        if showmsg:
            print('SolidModel taking action', action)

    def configure(self, props, globdat):

        # Configure the material
        matprops = props[MATERIAL]
        self._mat = new_material(matprops)
        self._mat.configure(matprops, globdat)

        # Get shape and element info
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][TYPE], props[SHAPE][INTSCHEME])
        egroup = globdat[gn.EGROUPS][props[ELEMENTS]]
        self._elems = egroup.get_elements()
        self._ielems = egroup.get_indices()
        self._nodes = self._elems.get_nodes()

        # Make sure the shape rank and mesh rank are identitcal
        if self._shape.global_rank() != globdat[gn.MESHRANK]:
            raise RuntimeError('ElasticModel: Shape rank must agree with mesh rank')

        # Get basic dimensionality info
        self._rank = self._shape.global_rank()
        self._ipcount = self._shape.ipoint_count()
        self._dofcount = self._rank * self._shape.node_count()
        self._strcount = self._rank * (self._rank + 1) // 2  # 1-->1, 2-->3, 3-->6

        self._thickness = 1.0
        if self._rank == 2:
            self._thickness = props.get(THICKNESS_PROP, self._thickness)
            self._thickness = pu.soft_cast(self._thickness, float)

        # Create a new dof for every node and dof type
        for doftype in DOFTYPES[0:self._rank]:
            globdat[gn.DOFSPACE].add_type(doftype)
            for inode in self._elems.get_unique_nodes_of(self._ielems):
                globdat[gn.DOFSPACE].add_dof(inode, doftype)

    def _get_matrix(self, params, globdat):

        for ielem in self._ielems:
            # Get the nodal coordinates of each element
            inodes = self._elems.get_elem_nodes(ielem)
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = self._nodes.get_some_coords(inodes)

            # Get the gradients, weights and coordinates of each integration point
            grads, weights = self._shape.get_shape_gradients(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            if self._rank == 2:
                weights *= self._thickness

            # Reset the element stiffness matrix
            elmat = np.zeros((self._dofcount, self._dofcount))

            for ip in range(self._ipcount):
                # Get the B and D matrices for each integration point
                B = self._get_B_matrix(grads[:, :, ip])
                D = self._mat.stiff_at_point(ipcoords[:, ip])

                # Compute the element stiffness matrix
                elmat += weights[ip] * np.matmul(np.transpose(B), np.matmul(D, B))

            # Add the element stiffness matrix to the global stiffness matrix
            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elmat

    def _get_mass_matrix(self, params, globdat):
        unit_matrix = params.get(pn.UNITMATRIX, False)

        if unit_matrix:
            M = np.identity(self._rank)

        for ielem in self._ielems:
            # Get the nodal coordinates of each element
            inodes = self._elems.get_elem_nodes(ielem)
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = self._nodes.get_some_coords(inodes)

            # Get the shape functions, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            weights = self._shape.get_integration_weights(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            if self._rank == 2:
                weights *= self._thickness

            # Reset the element mass matrix
            elmat = np.zeros((self._dofcount, self._dofcount))

            for ip in range(self._ipcount):
                # Get the N and M matrices for each integration point
                N = self._get_N_matrix(sfuncs[:, ip])

                if not unit_matrix:
                    M = self._mat.mass_at_point(ipcoords[:, ip])

                # Compute the element mass matrix
                elmat += weights[ip] * np.matmul(np.transpose(N), np.matmul(M, N))

            # Add the element mass matrix to the global mass matrix
            params[pn.MATRIX2][np.ix_(idofs, idofs)] += elmat

    def _get_strains(self, params, globdat):
        xtable = params[pn.TABLE]
        tbwts = params[pn.TABLEWEIGHTS]

        # Convert the table to an XTable and store the original class
        cls_ = xtable.__class__
        xtable.__class__ = XTable

        # Get the STATE0 vector if no custom displacement field is provided
        if pn.SOLUTION in params:
            disp = params[pn.SOLUTION]
        else:
            disp = globdat[gn.STATE0]

        # Add the columns of all stress components to the table
        if self._rank == 1:
            jcols = xtable.add_columns(['xx'])
        elif self._rank == 2:
            jcols = xtable.add_columns(['xx', 'yy', 'xy'])
        elif self._rank == 3:
            jcols = xtable.add_columns(['xx', 'yy', 'zz', 'xy', 'yz', 'zx'])

        for ielem in self._ielems:
            # Get the nodal coordinates of each element
            inodes = self._elems.get_elem_nodes(ielem)
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = self._nodes.get_some_coords(inodes)

            # Get the shape functions, gradients, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

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

        # Divide the strains by the shape function weights
        for jcol in jcols:
            values = xtable.get_col_values(None, jcol)
            xtable.set_col_values(None, jcol, values / tbwts)

        # Convert the table back to the original class
        xtable.__class__ = cls_

    def _get_stresses(self, params, globdat):
        xtable = params[pn.TABLE]
        tbwts = params[pn.TABLEWEIGHTS]

        # Convert the table to an XTable and store the original class
        cls_ = xtable.__class__
        xtable.__class__ = XTable

        # Get the STATE0 vector if no custom displacement field is provided
        if pn.SOLUTION in params:
            disp = params[pn.SOLUTION]
        else:
            disp = globdat[gn.STATE0]

        # Add the columns of all stress components to the table
        if self._rank == 1:
            jcols = xtable.add_columns(['xx'])
        elif self._rank == 2:
            jcols = xtable.add_columns(['xx', 'yy', 'xy'])
        elif self._rank == 3:
            jcols = xtable.add_columns(['xx', 'yy', 'zz', 'xy', 'yz', 'zx'])

        for ielem in self._ielems:
            # Get the nodal coordinates of each element
            inodes = self._elems.get_elem_nodes(ielem)
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
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

        # Divide the stresses by the shape function weights
        for jcol in jcols:
            values = xtable.get_col_values(None, jcol)
            xtable.set_col_values(None, jcol, values / tbwts)

        # Convert the table back to the original class
        xtable.__class__ = cls_

    def _get_stiffness(self, params, globdat):
        xtable = params[pn.TABLE]
        tbwts = params[pn.TABLEWEIGHTS]

        # Convert the table to an XTable and store the original class
        cls_ = xtable.__class__
        xtable.__class__ = XTable

        # Add the column of the Young's modulus to the table
        jcol = xtable.add_column('')

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

        # Divide the stresses by the shape function weights
        values = xtable.get_col_values(None, jcol)
        xtable.set_col_values(None, jcol, values / tbwts)

        # Convert the table back to the original class
        xtable.__class__ = cls_

    def _get_N_matrix(self, sfuncs):
        N = np.zeros((self._rank, self._dofcount))
        for i in range(self._rank):
            N[i, i::self._rank] = sfuncs.transpose()
        return N

    def _get_B_matrix(self, grads):
        B = np.zeros((self._strcount, self._dofcount))
        if self._rank == 1:
            B = grads.transpose()
        elif self._rank == 2:
            for inode in range(self._shape.node_count()):
                i = 2 * inode
                gi = grads[inode, :]
                B[0:3, i:(i + 2)] = [
                    [gi[0], 0.],
                    [0., gi[1]],
                    [gi[1], gi[0]]
                ]
        elif self._rank == 3:
            B = np.zeros((6, self._dofcount))
            for inode in range(self._shape.node_count()):
                i = 3 * inode
                gi = grads[inode, :]
                B[0:6, i:(i + 3)] = [
                    [gi[0], 0., 0.],
                    [0., gi[1], 0.],
                    [0., 0., gi[2]],
                    [gi[1], gi[0], 0.],
                    [0., gi[2], gi[1]],
                    [gi[2], 0., gi[0]]
                ]
        return B


def declare(factory):
    factory.declare_model('Solid', SolidModel)
