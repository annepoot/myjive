import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

from myjive.app import Module
from myjive.names import GlobNames as gn
from myjive.util import Table, to_xtable
from myjive.util.proputils import check_value


__all__ = ["ViewModule"]


class ViewModule(Module):
    @Module.save_config
    def configure(
        self,
        globdat,
        *,
        plotType,
        tables=[],
        elemTables=[],
        comps=[],
        elemComps=[],
        scale=0.0,
        line={},
        fill={},
        colorbar={},
        save={},
        figure={},
        axes={}
    ):
        # Validate input arguments
        check_value(self, plotType, ["node", "elem"])
        self._plottype = plotType
        self._ntables = tables
        self._etables = elemTables
        self._ncomps = comps
        self._ecomps = elemComps
        self._scale = scale
        self._lineprops = line
        self._fillprops = fill
        self._cbarprops = colorbar
        self._saveprops = save
        self._figprops = figure
        self._axprops = axes

    def init(self, globdat):
        pass

    def run(self, globdat):
        return "ok"

    def shutdown(self, globdat):
        # Go through the tables and comps
        for name, comp in zip(self._ntables, self._ncomps):
            if self._plottype == "node":
                self._write_node_table(name, globdat)
            elif self._plottype == "elem":
                self._write_elem_table(name, globdat)
            else:
                assert False

            self._generate_plot(name, comp, globdat)

    def _write_node_table(self, name, globdat):
        if name not in globdat[gn.TABLES]:
            nodecount = len(globdat[gn.NSET])
            table = Table(size=nodecount)
            tbwts = np.zeros(nodecount)

            for model in self.get_relevant_models("GETTABLE", globdat[gn.MODELS]):
                table, tbwts = model.GETTABLE(name, table, tbwts, globdat)

            to_xtable(table)

            for jcol in range(table.column_count()):
                values = table.get_col_values(None, jcol)
                table.set_col_values(None, jcol, values / tbwts)

            table.to_table()
            globdat[gn.TABLES][name] = table

    def _write_elem_table(self, name, globdat):
        if name not in globdat[gn.ELEMTABLES]:
            elemcount = len(globdat[gn.ESET])
            table = Table(size=elemcount)

            for model in self.get_relevant_models("GETELEMTABLE", globdat[gn.MODELS]):
                table = model.GETELEMTABLE(name, table, globdat)

            globdat[gn.ELEMTABLES][name] = table

    def _generate_plot(self, name, comp, globdat):
        # Get the necessary info from globdat
        nodes = globdat[gn.NSET]
        elems = globdat[gn.ESET]
        dofs = globdat[gn.DOFSPACE]
        types = dofs.get_types()
        shape = globdat[gn.MESHSHAPE]
        disp = globdat[gn.STATE0]

        if self._plottype == "node":
            array = globdat[gn.TABLES][name][comp]
            if len(array) != len(nodes):
                raise ValueError("Mismatch between table size and nodeset size")
        elif self._plottype == "elem":
            array = globdat[gn.ELEMTABLES][name][comp]
            if len(array) != len(elems):
                raise ValueError("Mismatch between table size and elementset size")
        else:
            assert False

        x = np.zeros(len(nodes))
        y = np.zeros(len(nodes))

        el = self._triangulate_elements(elems, shape)

        for inode, coords in enumerate(nodes):
            x[inode] = coords[0]
            y[inode] = coords[1]

        dx = np.copy(x)
        dy = np.copy(y)

        if self._scale > 0 and "dx" in types and "dy" in types:
            for n in range(len(nodes)):
                idofs = dofs.get_dofs([n], ["dx", "dy"])
                du = disp[idofs]
                dx[n] += self._scale * du[0]
                dy[n] += self._scale * du[1]

        fig, ax = plt.subplots(**self._figprops)
        ax.set_axis_off()
        ax.set_aspect("equal", adjustable="datalim")

        triang = Triangulation(dx, dy, el)

        if self._plottype == "node":
            mappable = ax.tricontourf(triang, array, **self._fillprops)
        elif self._plottype == "elem":
            mappable = ax.tripcolor(triang, array, **self._fillprops)

        if "show" in self._cbarprops:
            if self._cbarprops["show"]:
                self._cbarprops.pop("show")
                cmin, cmax = mappable.get_clim()
                ticks = np.linspace(cmin, cmax, 5, endpoint=True)
                cbar = plt.colorbar(mappable, ticks=ticks, **self._cbarprops)
                cbar.formatter.set_powerlimits((0, 0))
                self._cbarprops["show"] = True

        if "linewidth" in self._lineprops:
            ax.triplot(triang, **self._lineprops)

        ax.set(**self._axprops)

        if "fname" in self._saveprops:
            fname = self._saveprops["fname"]
            dirname = os.path.dirname(fname)
            if len(dirname) > 0 and not os.path.isdir(dirname):
                os.makedirs(dirname)
            plt.savefig(fname, *self._saveprops)

        plt.show()

    def _triangulate_elements(self, elems, shape):
        if shape == "Triangle3":
            nelem = len(elems)
        elif shape == "Triangle6":
            nelem = len(elems) * 4
        elif shape == "Quad4":
            nelem = len(elems) * 2
        elif shape == "Quad9":
            nelem = len(elems) * 8
        else:
            raise ValueError("ViewModule only supports triangles for now")

        el = np.zeros((nelem, 3), dtype=int)

        for e, inodes in enumerate(elems):
            if shape == "Triangle3":
                el[e, :] = inodes
            elif shape == "Triangle6":
                el[4 * e + 0, :] = inodes[[0, 3, 5]]
                el[4 * e + 1, :] = inodes[[1, 4, 3]]
                el[4 * e + 2, :] = inodes[[2, 5, 4]]
                el[4 * e + 3, :] = inodes[[3, 4, 5]]
            elif shape == "Quad4":
                el[2 * e + 0, :] = inodes[[0, 1, 3]]
                el[2 * e + 1, :] = inodes[[1, 2, 3]]
            elif shape == "Quad9":
                el[8 * e + 0, :] = inodes[[0, 4, 7]]
                el[8 * e + 1, :] = inodes[[4, 8, 7]]
                el[8 * e + 2, :] = inodes[[1, 8, 4]]
                el[8 * e + 3, :] = inodes[[1, 5, 8]]
                el[8 * e + 4, :] = inodes[[2, 5, 6]]
                el[8 * e + 5, :] = inodes[[5, 6, 8]]
                el[8 * e + 6, :] = inodes[[3, 8, 6]]
                el[8 * e + 7, :] = inodes[[3, 7, 8]]

        return el
