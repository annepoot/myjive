import numpy as np

from jive.fem.names import Actions as act
from jive.fem.names import ParamNames as pn
from jive.fem.names import GlobNames as gn
from jive.model.model import Model
from jive.util.xtable import XTable, to_xtable

ELEMENTS = 'elements'
YOUNG = 'young'
RHO = 'rho'
THICKNESS = 'thickness'
POISSON = 'poisson'
SHAPE = 'shape'
TYPE = 'type'
INTSCHEME = 'intScheme'
STATE = 'state'
DOFTYPES = ['dx', 'dy', 'dz']
PE_STATE = 'plane_strain'
PS_STATE = 'plane_stress'


class ElasticModel(Model):
    def take_action(self, action, params, globdat):
        showmsg = True

        if action == act.GETMATRIX0:
            self._get_matrix(params, globdat)
        elif action == act.GETMATRIX2:
            self._get_mass_matrix(params, globdat)
        elif action == act.GETEXTFORCE:
            self._get_body_force(params, globdat)
        elif action == act.GETTABLE:
            if 'stress' in params[pn.TABLENAME]:
                self._get_stresses(params, globdat)
            elif 'strain' in params[pn.TABLENAME]:
                self._get_strains(params, globdat)
            else:
                showmsg = False
        else:
            showmsg = False

        verbose = params.get(pn.VERBOSE, True)

        if showmsg and verbose:
            print('ElasticModel taking action', action)

    def configure(self, props, globdat):
        # This function gets only the core values from props

        # Get basic parameter values
        self._young = float(props[YOUNG])
        self._rho = float(props.get(RHO,0))

        if globdat[gn.MESHRANK] > 1:
            self._poisson = float(props[POISSON])
        if globdat[gn.MESHRANK] == 2:
            self._thickness = float(props.get(THICKNESS,1))
            self._state = props[STATE]
            if not self._state in (PE_STATE, PS_STATE):
                raise RuntimeError('ElasticModel: state in 2d should be plane_strain or plane_stress')

        # Get shape and element info
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(props[SHAPE][TYPE], props[SHAPE][INTSCHEME])
        egroup = globdat[gn.EGROUPS][props[ELEMENTS]]
        self._elems = egroup.get_elements()
        self._ielems = egroup.get_indices()
        self._nodes = self._elems.get_nodes()

        # Make sure the shape rank and mesh rank are identitcal
        if self._shape.global_rank() != globdat[gn.MESHRANK]:
            raise RuntimeError('ElasticModel: Shape rank must agree with mesh rank')

        # The rest of the configuration happens in configure_noprops
        self._configure_noprops(globdat)

    def _configure_noprops(self, globdat):
        # This function gets additional info from self and globdat
        # It has been split off from configure() to allow it to be used in inherited classes as well.

        # Get basic dimensionality info
        self._rank = self._shape.global_rank()
        self._ipcount = self._shape.ipoint_count()
        self._dofcount = self._rank * self._shape.node_count()
        self._strcount = self._rank * (self._rank + 1) // 2   # 1-->1, 2-->3, 3-->6

        # Create a new dof for every node and dof type
        nodes = np.unique([node for elem in self._elems for node in elem.get_nodes()])
        for doftype in DOFTYPES[0:self._rank]:
            globdat[gn.DOFSPACE].add_type(doftype)
            for node in nodes:
                globdat[gn.DOFSPACE].add_dof(node, doftype)

    def _get_matrix(self, params, globdat):

        for elem in self._elems:
            # Get the nodal coordinates of each element
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = self._nodes.get_some_coords(inodes)

            # Get the gradients, weights and coordinates of each integration point
            grads, weights = self._shape.get_shape_gradients(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            # Reset the element stiffness matrix
            elmat = np.zeros((self._dofcount, self._dofcount))

            for ip in range(self._ipcount):
                # Get the B and D matrices for each integration point
                B = self._get_B_matrix(grads[:,:,ip])
                D = self._get_D_matrix(ipcoords[:,ip])

                # Compute the element stiffness matrix
                elmat += weights[ip] * np.matmul(np.transpose(B), np.matmul(D, B))

            # Add the element stiffness matrix to the global stiffness matrix
            params[pn.MATRIX0][np.ix_(idofs, idofs)] += elmat

    def _get_mass_matrix(self, params, globdat):

        for elem in self._elems:
            # Get the nodal coordinates of each element
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = self._nodes.get_some_coords(inodes)

            # Get the shape functions, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            weights = self._shape.get_integration_weights(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            # Reset the element mass matrix
            elmat = np.zeros((self._dofcount, self._dofcount))

            for ip in range(self._ipcount):
                # Get the N and M matrices for each integration point
                N = self._get_N_matrix(sfuncs[:,ip])
                M = self._get_M_matrix(ipcoords[:,ip])

                # Compute the element mass matrix
                elmat += weights[ip] * np.matmul(np.transpose(N), np.matmul(M, N))

            # Add the element mass matrix to the global mass matrix
            params[pn.MATRIX2][np.ix_(idofs, idofs)] += elmat

    def _get_body_force(self, params, globdat):
        if self._rank == 2:

            for elem in self._elems:
                # Get the nodal coordinates of each element
                inodes = elem.get_nodes()
                idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
                coords = self._nodes.get_some_coords(inodes)

                # Get the shape functions, weights and coordinates of each integration point
                sfuncs = self._shape.get_shape_functions()
                weights = self._shape.get_integration_weights(coords)
                ipcoords = self._shape.get_global_integration_points(coords)

                # Reset the element force vector
                elfor = np.zeros(self._dofcount)

                for ip in range(self._ipcount):
                    # Get the N matrix and b vector for each integration point
                    N = self._get_N_matrix(sfuncs[:,ip])
                    b = self._get_b_vector(ipcoords[:,ip])

                    # Compute the element force vector
                    elfor += weights[ip] * np.matmul(np.transpose(N), b)

                # Add the element force vector to the global force vector
                params[pn.EXTFORCE][idofs] += elfor

        else:
            pass

    def _get_strains(self, params, globdat):
        table = params[pn.TABLE]
        tbwts = params[pn.TABLEWEIGHTS]

        # Convert the table to an XTable and store the original class
        xtable = to_xtable(table)

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

        for elem in self._elems:
            # Get the nodal coordinates of each element
            inodes = elem.get_nodes()
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
                B = self._get_B_matrix(grads[:,:,ip])

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
        params[pn.TABLE] = xtable.to_table()

    def _get_stresses(self, params, globdat):
        table = params[pn.TABLE]
        tbwts = params[pn.TABLEWEIGHTS]

        # Convert the table to an XTable and store the original class
        xtable = to_xtable(table)

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

        for elem in self._elems:
            # Get the nodal coordinates of each element
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0:self._rank])
            coords = self._nodes.get_some_coords(inodes)

            # Get the shape functions, gradients, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            grads, weights = self._shape.get_shape_gradients(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            # Get the nodal displacements
            eldisp = disp[idofs]

            # Reset the element stress matrix and weights
            elsig = np.zeros((self._shape.node_count(), self._strcount))
            elwts = np.zeros(self._shape.node_count())

            for ip in range(self._ipcount):
                # Get the B and D matrices for each integration point
                B = self._get_B_matrix(grads[:,:,ip])
                D = self._get_D_matrix(ipcoords[:,ip])

                if self._rank == 2:
                    D /= self._thickness

                # Get the strain and stress of the element in the integration point
                strain = np.matmul(B, eldisp)
                stress = np.matmul(D, strain)

                # Compute the element stress and weights
                elsig += np.outer(sfuncs[:, ip], stress)
                elwts += sfuncs[:, ip].flatten()

            # Add the element weights to the global weights
            tbwts[inodes] += elwts

            # Add the element stresses to the global stresses
            xtable.add_block(inodes, jcols, elsig)

        # Convert the table back to the original class
        params[pn.TABLE] = xtable.to_table()

    def _get_N_matrix(self, sfuncs):
        N = np.zeros((self._rank, self._dofcount))
        for i in range(self._rank):
            N[i,i::self._rank] = sfuncs.transpose()
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

    def _get_D_matrix(self, ipcoords):
        D = np.zeros((self._strcount, self._strcount))
        E = self._young
        if self._rank == 1:
            D[[0]] = E
            return D
        nu = self._poisson
        g = 0.5 * E / (1. + nu)

        if self._rank == 3:
            a = E * (1. - nu) / ((1. + nu) * (1. - 2. * nu))
            b = E * nu / ((1. + nu) * (1. - 2. * nu))
            D[:,:] = [
                [a, b, b, .0, .0, .0],
                [b, a, b, .0, .0, .0],
                [b, b, a, .0, .0, .0],
                [.0, .0, .0, g, .0, .0],
                [.0, .0, .0, .0, g, .0],
                [.0, .0, .0, .0, .0, g]
            ]

        elif self._rank == 2:
            if self._state == PE_STATE:
                a = E * (1. - nu) / ((1. + nu) * (1. - 2. * nu))
                b = E * nu / ((1. + nu) * (1. - 2. * nu))
                D[:,:] = [
                    [a, b, .0],
                    [b, a, .0],
                    [.0, .0, g]
                ]
            else:
                assert (self._state == PS_STATE)
                a = E / (1. - nu * nu)
                b = a * nu
                D[:,:] = [
                    [a, b, .0],
                    [b, a, .0],
                    [.0, .0, g]
                ]

            D *= self._thickness

        return D

    def _get_M_matrix(self, ipcoords):
        M = self._rho * np.identity(self._rank)

        if self._rank == 2:
            M *= self._thickness

        return M

    def _get_b_vector(self, ipcoords):

        if self._rank == 3:
            gravity = np.array([0, -1, 0])
            b = self._rho * gravity

        elif self._rank == 2:
            gravity = np.array([0, -1])
            b = self._rho * self._thickness * gravity

        else:
            raise RuntimeError('ElasticModel: self weight can only be computed in 2d or 3d')

        return b

def declare(factory):
    factory.declare_model('Elastic', ElasticModel)
