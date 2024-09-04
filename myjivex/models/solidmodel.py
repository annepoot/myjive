import numpy as np
from numba import njit

from myjive.names import GlobNames as gn
from myjive.model.model import Model
from ..materials import new_material
import myjive.util.proputils as pu
from myjive.util import to_xtable
from myjive.util.proputils import check_dict, split_off_type

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
        if "displacement" in name or "solution" in name:
            table, tbwts = self._get_solution_by_node(table, tbwts, globdat, **kwargs)
        if "strain" in name:
            table, tbwts = self._get_strain_by_node(table, tbwts, globdat, **kwargs)
        elif "stress" in name:
            table, tbwts = self._get_stress_by_node(table, tbwts, globdat, **kwargs)
        elif "stiffness" in name:
            table, tbwts = self._get_stiffness_by_node(table, tbwts, globdat)
        elif "size" in name:
            table, tbwts = self._get_elem_size_by_node(table, tbwts, globdat)
        return table, tbwts

    def GETELEMTABLE(self, name, table, globdat, **kwargs):
        if "strain" in name:
            table = self._get_strain_by_elem(table, globdat, **kwargs)
        elif "stress" in name:
            table = self._get_stress_by_elem(table, globdat, **kwargs)
        elif "stiffness" in name:
            table = self._get_stiffness_by_elem(table, globdat, **kwargs)
        elif "size" in name:
            table = self._get_elem_size_by_elem(table, globdat, **kwargs)
        return table

    @Model.save_config
    def configure(self, globdat, *, shape, elements, material, thickness=1.0):
        # Validate input arguments
        check_dict(self, shape, [TYPE, INTSCHEME])
        check_dict(self, material, [TYPE])
        self._thickness = thickness

        # Configure the material
        mattype, matprops = split_off_type(material)
        self._mat = new_material(mattype, "material")
        self._mat.configure(globdat, **matprops)
        self._config["material"] = self._mat.get_config()

        # Get shape and element info
        self._shape = globdat[gn.SHAPEFACTORY].get_shape(shape[TYPE], shape[INTSCHEME])
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
                B_elem = self._get_B_matrix(grads[ip])

                if not unit_matrix:
                    D_elem = self._mat.stiff_at_point(ipcoords[ip])

                # Compute the element stiffness matrix
                elmat += weights[ip] * B_elem.T @ D_elem @ B_elem

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
                N_elem = self._get_N_matrix(sfuncs[ip])

                if not unit_matrix:
                    M_elem = self._mat.mass_at_point(ipcoords[ip])

                # Compute the element mass matrix
                elmat += weights[ip] * N_elem.T @ M_elem @ N_elem

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
                B_elem = self._get_B_matrix(grads[ip])

                # Compute the element strain and weights
                for i in range(self._strcount):
                    elbmat[i * node_count : (i + 1) * node_count, :] += np.outer(
                        sfuncs[ip], B_elem[i, :]
                    )
                    elwts[i * node_count : (i + 1) * node_count] += sfuncs[ip].flatten()

            # Get the node index vector
            node_idx = np.zeros(node_count * self._strcount, dtype=int)
            for i in range(self._strcount):
                node_idx[i * node_count : (i + 1) * node_count] = inodes + nc * i

            # Add the element strain matrix to the global strain matrix
            B[np.ix_(node_idx, idofs)] += elbmat

            # Add the element weights to the global weights
            wts[node_idx] += elwts

        return B, wts

    def _get_solution_by_node(self, table, tbwts, globdat, solution=None):
        disp = globdat[gn.STATE0] if solution is None else solution

        comps = self._get_solution_comps()
        # table, tbwts = self._fill_by_node(table, tbwts, comps, func, globdat, disp=disp)

        xtable = to_xtable(table)
        jcols = xtable.add_columns(comps)

        for inode, node in enumerate(self._nodes):
            idofs = globdat[gn.DOFSPACE].get_dofs([inode], DOFTYPES[0 : self._rank])
            sol = disp[idofs]
            tbwts[inode] += 1
            xtable.add_row_values(inode, jcols, sol)

        table = xtable.to_table()
        return table, tbwts

    def _get_strain_by_node(self, table, tbwts, globdat, solution=None):
        disp = globdat[gn.STATE0] if solution is None else solution

        comps = self._get_gradient_comps()
        func = self._get_node_strain
        table, tbwts = self._fill_by_node(table, tbwts, comps, func, globdat, disp=disp)
        return table, tbwts

    def _get_stress_by_node(self, table, tbwts, globdat, solution=None):
        disp = globdat[gn.STATE0] if solution is None else solution

        comps = self._get_gradient_comps()
        func = self._get_node_stress
        table, tbwts = self._fill_by_node(table, tbwts, comps, func, globdat, disp=disp)
        return table, tbwts

    def _get_stiffness_by_node(self, table, tbwts, globdat):
        comps = [""]
        func = self._get_node_stiffness
        table, tbwts = self._fill_by_node(table, tbwts, comps, func, globdat)
        return table, tbwts

    def _get_elem_size_by_node(self, table, tbwts, globdat):
        comps = [""]
        func = self._get_node_size
        table, tbwts = self._fill_by_node(table, tbwts, comps, func, globdat)
        return table, tbwts

    def _get_strain_by_elem(self, table, globdat, solution=None):
        disp = globdat[gn.STATE0] if solution is None else solution

        comps = self._get_gradient_comps()
        func = self._get_elem_strain
        table = self._fill_by_elem(table, comps, func, globdat, disp=disp)
        return table

    def _get_stress_by_elem(self, table, globdat, solution=None):
        disp = globdat[gn.STATE0] if solution is None else solution

        comps = self._get_gradient_comps()
        func = self._get_elem_stress
        table = self._fill_by_elem(table, comps, func, globdat, disp=disp)
        return table

    def _get_stiffness_by_elem(self, table, globdat, solution=None):
        if self._ipcount != 1:
            raise ValueError("element stiffness only works with a single Gauss point")

        comps = [""]
        func = self._get_elem_stiffness
        table = self._fill_by_elem(table, comps, func, globdat)
        return table

    def _get_elem_size_by_elem(self, table, globdat, solution=None):
        comps = [""]
        func = self._get_elem_size
        table = self._fill_by_elem(table, comps, func, globdat)
        return table

    def _fill_by_elem(self, table, comps, func, globdat, **fargs):
        xtable = to_xtable(table)
        jcols = xtable.add_columns(comps)

        for ielem in self._ielems:
            strain = func(ielem, globdat, **fargs)
            xtable.set_row_values(ielem, jcols, strain)

        table = xtable.to_table()
        return table

    def _fill_by_node(self, table, tbwts, comps, func, globdat, **fargs):
        xtable = to_xtable(table)
        jcols = xtable.add_columns(comps)

        for ielem in self._ielems:
            inodes = self._elems.get_elem_nodes(ielem)
            field, weights = func(ielem, inodes, globdat, **fargs)

            tbwts[inodes] += weights
            xtable.add_block(inodes, jcols, field)

        table = xtable.to_table()
        return table, tbwts

    def _get_elem_strain(self, ielem, globdat, disp=None):
        inodes = self._elems.get_elem_nodes(ielem)
        idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0 : self._rank])
        coords = self._nodes.get_some_coords(inodes)
        grads, _ = self._shape.get_shape_gradients(coords)

        eldisp = disp[idofs]
        strains = np.zeros((self._ipcount, self._strcount))
        for ip in range(self._ipcount):
            strains[ip] = self._get_ip_strain(ip, grads, eldisp)
        if not np.allclose(strains, strains[0]):
            raise RuntimeError("Non-constant strain over element")

        return strains[0]

    def _get_elem_stress(self, ielem, globdat, disp=None):
        inodes = self._elems.get_elem_nodes(ielem)
        idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0 : self._rank])
        coords = self._nodes.get_some_coords(inodes)
        grads, _ = self._shape.get_shape_gradients(coords)
        ipcoords = self._shape.get_global_integration_points(coords)

        eldisp = disp[idofs]
        stresses = np.zeros((self._ipcount, self._strcount))
        for ip in range(self._ipcount):
            stresses[ip] = self._get_ip_stress(ip, ipcoords, grads, eldisp)
        if not np.allclose(stresses, stresses[0]):
            raise RuntimeError("Non-constant stresses over element")

        return stresses[0]

    def _get_elem_stiffness(self, ielem, globdat):
        inodes = self._elems.get_elem_nodes(ielem)
        coords = self._nodes.get_some_coords(inodes)
        grads, _ = self._shape.get_shape_gradients(coords)
        ipcoords = self._shape.get_global_integration_points(coords)

        stiffness = self._mat._get_E(ipcoords[0])
        return stiffness

    def _get_elem_size(self, ielem, globdat):
        inodes = self._elems.get_elem_nodes(ielem)
        coords = self._nodes.get_some_coords(inodes)

        max_edge = 0.0
        for i, inode in enumerate(inodes):
            for j, jnode in enumerate(inodes):
                icoords = coords[i]
                jcoords = coords[j]
                edge = np.sqrt(np.sum((icoords - jcoords) ** 2))
                if edge > max_edge:
                    max_edge = edge

        return max_edge

    def _get_node_strain(self, ielem, inodes, globdat, disp=None):
        idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0 : self._rank])
        coords = self._nodes.get_some_coords(inodes)
        sfuncs = self._shape.get_shape_functions()
        grads, weights = self._shape.get_shape_gradients(coords)

        eldisp = disp[idofs]
        eleps = np.zeros((self._shape.node_count(), self._strcount))
        elwts = np.zeros(self._shape.node_count())

        for ip in range(self._ipcount):
            strain = self._get_ip_strain(ip, grads, eldisp)
            eleps += np.outer(sfuncs[ip], strain)
            elwts += sfuncs[ip]

        return eleps, elwts

    def _get_node_stress(self, ielem, inodes, globdat, disp=None):
        idofs = globdat[gn.DOFSPACE].get_dofs(inodes, DOFTYPES[0 : self._rank])
        coords = self._nodes.get_some_coords(inodes)
        sfuncs = self._shape.get_shape_functions()
        grads, weights = self._shape.get_shape_gradients(coords)
        ipcoords = self._shape.get_global_integration_points(coords)

        if self._rank == 2:
            weights *= self._thickness

        eldisp = disp[idofs]
        elsig = np.zeros((self._shape.node_count(), self._strcount))
        elwts = np.zeros(self._shape.node_count())

        for ip in range(self._ipcount):
            stress = self._get_ip_stress(ip, ipcoords, grads, eldisp)
            elsig += np.outer(sfuncs[ip], stress)
            elwts += sfuncs[ip].flatten()

        return elsig, elwts

    def _get_node_stiffness(self, ielem, inodes, globdat):
        coords = self._nodes.get_some_coords(inodes)
        sfuncs = self._shape.get_shape_functions()
        ipcoords = self._shape.get_global_integration_points(coords)

        elyoung = np.zeros((self._shape.node_count(), 1))
        elwts = np.zeros(self._shape.node_count())

        for ip in range(self._ipcount):
            E = self._mat._get_E(ipcoords[ip])
            elyoung[:, 0] += E * sfuncs[ip]
            elwts += sfuncs[ip].flatten()

        return elyoung, elwts

    def _get_node_size(self, ielem, inodes, globdat):
        max_edge = self._get_elem_size(ielem, globdat)
        elsize = max_edge * np.ones((self._shape.node_count(), 1))
        elwts = np.ones(self._shape.node_count())
        return elsize, elwts

    def _get_solution_comps(self):
        if self._rank == 1:
            comps = ["dx"]
        elif self._rank == 2:
            comps = ["dx", "dy"]
        elif self._rank == 3:
            comps = ["dx", "dy", "dz"]
        return comps

    def _get_gradient_comps(self):
        if self._rank == 1:
            comps = ["xx"]
        elif self._rank == 2:
            comps = ["xx", "yy", "xy"]
        elif self._rank == 3:
            comps = ["xx", "yy", "zz", "xy", "yz", "zx"]
        return comps

    def _get_N_matrix(self, sfuncs):
        return self._get_N_matrix_jit(sfuncs, self._dofcount, self._rank)

    @staticmethod
    @njit
    def _get_N_matrix_jit(sfuncs, _dofcount, _rank):
        N = np.zeros((_rank, _dofcount))
        for i in range(_rank):
            N[i, i::_rank] = sfuncs
        return N

    def _get_B_matrix(self, grads):
        return self._get_B_matrix_jit(
            grads, self._strcount, self._dofcount, self._shape.node_count(), self._rank
        )

    @staticmethod
    @njit
    def _get_B_matrix_jit(grads, _strcount, _dofcount, _nodecount, _rank):
        B = np.zeros((_strcount, _dofcount))
        if _rank == 1:
            B = grads
        elif _rank == 2:
            for inode in range(_nodecount):
                i = 2 * inode
                gi = grads[:, inode]
                B[0, i + 0] = gi[0]
                B[1, i + 1] = gi[1]
                B[2, i + 0] = gi[1]
                B[2, i + 1] = gi[0]
        elif _rank == 3:
            for inode in range(_nodecount):
                i = 3 * inode
                gi = grads[:, inode]
                B[0, i + 0] = gi[0]
                B[1, i + 1] = gi[1]
                B[2, i + 2] = gi[2]
                B[3, i + 0] = gi[1]
                B[3, i + 1] = gi[0]
                B[4, i + 1] = gi[2]
                B[4, i + 2] = gi[1]
                B[5, i + 0] = gi[2]
                B[5, i + 2] = gi[0]
        return B

    def _get_ip_strain(self, ip, grads, eldisp):
        B = self._get_B_matrix(grads[ip])
        strain = np.matmul(B, eldisp)
        return strain

    def _get_ip_stress(self, ip, ipcoords, grads, eldisp):
        strain = self._get_ip_strain(ip, grads, eldisp)
        stress = self._mat.stress_at_point(strain, ipcoords[ip])
        return stress
