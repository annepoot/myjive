import numpy as np

from myjive.names import GlobNames as gn
from myjive.model.model import Model
from myjive.util import Table, to_xtable
from myjive.util.proputils import mandatory_argument, mandatory_dict, optional_argument

TYPE = "type"
INTSCHEME = "intScheme"
DOFTYPES = ["dx", "dy", "dz"]
PE_STATE = "plane_strain"
PS_STATE = "plane_stress"

__all__ = ["ElasticModel"]


class ElasticModel(Model):
    def GETMATRIX0(self, K, globdat, **kwargs):
        K = self._get_matrix(K, globdat, **kwargs)
        return K

    def GETMATRIX2(self, M, globdat, **kwargs):
        M = self._get_mass_matrix(M, globdat, **kwargs)
        return M

    def GETEXTFORCE(self, f_ext, globdat, **kwargs):
        f_ext = self._get_body_force(f_ext, globdat, **kwargs)
        return f_ext

    def GETTABLE(self, name, table, tbwts, globdat, **kwargs):
        if "stress" in name:
            table, tbwts = self._get_stresses(table, tbwts, globdat, **kwargs)
        elif "strain" in name:
            table, tbwts = self._get_strains(table, tbwts, globdat, **kwargs)
        return table, tbwts

    def configure(self, globdat, **props):

        # Get props
        shapeprops = mandatory_dict(
            self, props, "shape", mandatory_keys=[TYPE, INTSCHEME]
        )
        elements = mandatory_argument(self, props, "elements")
        self._young = mandatory_argument(self, props, "young")
        self._rho = optional_argument(self, props, "rho", default=0.0)

        if globdat[gn.MESHRANK] > 1:
            self._poisson = mandatory_argument(self, props, "poisson")
        if globdat[gn.MESHRANK] == 2:
            self._thickness = optional_argument(self, props, "thickness", default=1.0)
            self._state = mandatory_argument(self, props, "state")
            if not self._state in (PE_STATE, PS_STATE):
                raise RuntimeError(
                    "ElasticModel: state in 2d should be plane_strain or plane_stress"
                )

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

        # The rest of the configuration happens in configure_noprops
        self._configure_noprops(globdat)

    def _configure_noprops(self, globdat):
        # This function gets additional info from self and globdat
        # It has been split off from configure() to allow it to be used in inherited classes as well.

        # Get basic dimensionality info
        self._rank = self._shape.global_rank()
        self._ipcount = self._shape.ipoint_count()
        self._dofcount = self._rank * self._shape.node_count()
        self._strcount = self._rank * (self._rank + 1) // 2  # 1-->1, 2-->3, 3-->6

        # Create a new dof for every node and dof type
        nodes = np.unique([node for elem in self._elems for node in elem.get_nodes()])
        for doftype in DOFTYPES[0 : self._rank]:
            globdat[gn.DOFSPACE].add_type(doftype)
            for node in nodes:
                globdat[gn.DOFSPACE].add_dof(node, doftype)

    def _get_matrix(self, K, globdat):
        if K is None:
            dc = globdat[gn.DOFSPACE].dof_count()
            K = np.zeros((dc, dc))

        for elem in self._elems:
            # Get the nodal coordinates of each element
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0 : self._rank])
            coords = self._nodes.get_some_coords(inodes)

            # Get the gradients, weights and coordinates of each integration point
            grads, weights = self._shape.get_shape_gradients(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            # Reset the element stiffness matrix
            elmat = np.zeros((self._dofcount, self._dofcount))

            for ip in range(self._ipcount):
                # Get the B and D matrices for each integration point
                B_elem = self._get_B_matrix(grads[:, :, ip])
                D_elem = self._get_D_matrix(ipcoords[:, ip])

                # Compute the element stiffness matrix
                elmat += weights[ip] * np.matmul(
                    np.transpose(B_elem), np.matmul(D_elem, B_elem)
                )

            # Add the element stiffness matrix to the global stiffness matrix
            K[np.ix_(idofs, idofs)] += elmat

        return K

    def _get_mass_matrix(self, M, globdat):
        if M is None:
            dc = globdat[gn.DOFSPACE].dof_count()
            M = np.zeros((dc, dc))

        for elem in self._elems:
            # Get the nodal coordinates of each element
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0 : self._rank])
            coords = self._nodes.get_some_coords(inodes)

            # Get the shape functions, weights and coordinates of each integration point
            sfuncs = self._shape.get_shape_functions()
            weights = self._shape.get_integration_weights(coords)
            ipcoords = self._shape.get_global_integration_points(coords)

            # Reset the element mass matrix
            elmat = np.zeros((self._dofcount, self._dofcount))

            for ip in range(self._ipcount):
                # Get the N and M matrices for each integration point
                N_elem = self._get_N_matrix(sfuncs[:, ip])
                M_elem = self._get_M_matrix(ipcoords[:, ip])

                # Compute the element mass matrix
                elmat += weights[ip] * np.matmul(
                    np.transpose(N_elem), np.matmul(M_elem, N_elem)
                )

            # Add the element mass matrix to the global mass matrix
            M[np.ix_(idofs, idofs)] += elmat

        return M

    def _get_body_force(self, f_ext, globdat):
        if f_ext is None:
            f_ext = np.zeros(globdat[gn.DOFSPACE].dof_count())

        if self._rank == 2:
            for elem in self._elems:
                # Get the nodal coordinates of each element
                inodes = elem.get_nodes()
                idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0 : self._rank])
                coords = self._nodes.get_some_coords(inodes)

                # Get the shape functions, weights and coordinates of each integration point
                sfuncs = self._shape.get_shape_functions()
                weights = self._shape.get_integration_weights(coords)
                ipcoords = self._shape.get_global_integration_points(coords)

                # Reset the element force vector
                elfor = np.zeros(self._dofcount)

                for ip in range(self._ipcount):
                    # Get the N matrix and b vector for each integration point
                    N = self._get_N_matrix(sfuncs[:, ip])
                    b = self._get_b_vector(ipcoords[:, ip])

                    # Compute the element force vector
                    elfor += weights[ip] * np.matmul(np.transpose(N), b)

                # Add the element force vector to the global force vector
                f_ext[idofs] += elfor

        return f_ext

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

        # Add the columns of all stress components to the table
        if self._rank == 1:
            jcols = xtable.add_columns(["xx"])
        elif self._rank == 2:
            jcols = xtable.add_columns(["xx", "yy", "xy"])
        elif self._rank == 3:
            jcols = xtable.add_columns(["xx", "yy", "zz", "xy", "yz", "zx"])

        for elem in self._elems:
            # Get the nodal coordinates of each element
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0 : self._rank])
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

        for elem in self._elems:
            # Get the nodal coordinates of each element
            inodes = elem.get_nodes()
            idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0 : self._rank])
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
                B = self._get_B_matrix(grads[:, :, ip])
                D = self._get_D_matrix(ipcoords[:, ip])

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
        table = xtable.to_table()

        return table, tbwts

    def _get_N_matrix(self, sfuncs):
        N = np.zeros((self._rank, self._dofcount))
        for i in range(self._rank):
            N[i, i :: self._rank] = sfuncs.transpose()
        return N

    def _get_B_matrix(self, grads):
        B = np.zeros((self._strcount, self._dofcount))
        if self._rank == 1:
            B = grads.transpose()
        elif self._rank == 2:
            for inode in range(self._shape.node_count()):
                i = 2 * inode
                gi = grads[inode, :]
                B[0:3, i : (i + 2)] = [[gi[0], 0.0], [0.0, gi[1]], [gi[1], gi[0]]]
        elif self._rank == 3:
            B = np.zeros((6, self._dofcount))
            for inode in range(self._shape.node_count()):
                i = 3 * inode
                gi = grads[inode, :]
                B[0:6, i : (i + 3)] = [
                    [gi[0], 0.0, 0.0],
                    [0.0, gi[1], 0.0],
                    [0.0, 0.0, gi[2]],
                    [gi[1], gi[0], 0.0],
                    [0.0, gi[2], gi[1]],
                    [gi[2], 0.0, gi[0]],
                ]
        return B

    def _get_D_matrix(self, ipcoords):
        D = np.zeros((self._strcount, self._strcount))
        E = self._young
        if self._rank == 1:
            D[[0]] = E
            return D
        nu = self._poisson
        g = 0.5 * E / (1.0 + nu)

        if self._rank == 3:
            a = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))
            b = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
            D[:, :] = [
                [a, b, b, 0.0, 0.0, 0.0],
                [b, a, b, 0.0, 0.0, 0.0],
                [b, b, a, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, g, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, g, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, g],
            ]

        elif self._rank == 2:
            if self._state == PE_STATE:
                a = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))
                b = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
                D[:, :] = [[a, b, 0.0], [b, a, 0.0], [0.0, 0.0, g]]
            else:
                assert self._state == PS_STATE
                a = E / (1.0 - nu * nu)
                b = a * nu
                D[:, :] = [[a, b, 0.0], [b, a, 0.0], [0.0, 0.0, g]]

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
            raise RuntimeError(
                "ElasticModel: self weight can only be computed in 2d or 3d"
            )

        return b
