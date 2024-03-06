import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

from jive.app import Module
from jive.names import GlobNames as gn
from jive.names import ParamNames as pn
from jive.names import Actions as act
from jive.util import Table, to_xtable

LINEWIDTH = "linewidth"
PLOT = "plot"
NCOLORS = "ncolors"
DEFORM = "deform"
SOLUTION = "solution"

__all__ = ["ViewModule"]


class ViewModule(Module):
    def init(self, props, globdat):
        self._scale = 0.0
        self._linewidth = 0.2
        self._pname = ""
        self._solutionname = ""
        self._ncolors = 100

        if self._name in props:
            myprops = props.get(self._name)

            if LINEWIDTH in myprops:
                self._linewidth = float(myprops[LINEWIDTH])
            if PLOT in myprops:
                self._pname = myprops[PLOT]
            if SOLUTION in myprops:
                self._solutionname = myprops[SOLUTION]
            if DEFORM in myprops:
                self._scale = myprops[DEFORM]
            if NCOLORS in myprops:
                self._ncolors = int(myprops[NCOLORS])

    def run(self, globdat):
        nodes = globdat[gn.NSET]
        elems = globdat[gn.ESET]

        if self._solutionname == "":
            self._solution = globdat[gn.STATE0]
        else:
            if "[" in self._solutionname:
                name = self._solutionname.split("[")[0]
                comp = int(self._solutionname.split("[")[1].split("]")[0])
                self._solution = globdat[name][:, comp]
            else:
                self._solution = globdat[self._solutionname]

        dofs = globdat[gn.DOFSPACE]
        types = dofs.get_types()

        assert (
            globdat[gn.MESHSHAPE] == "Triangle3"
        ), "ViewModule only supports triangles for now"

        x = np.zeros(len(nodes))
        y = np.zeros(len(nodes))
        el = np.zeros((len(elems), 3), dtype=int)

        for n, node in enumerate(nodes):
            coords = node.get_coords()

            x[n] = coords[0]
            y[n] = coords[1]

        for e, elem in enumerate(elems):
            inodes = elem.get_nodes()

            el[e, :] = inodes

        dx = np.copy(x)
        dy = np.copy(y)

        for n in range(len(nodes)):
            idofs = dofs.get_dofs([n], types)
            du = self._solution[idofs]

            if len(idofs) == 2:
                dx[n] += self._scale * du[0]
                dy[n] += self._scale * du[1]

        plt.figure()
        ax = plt.gca()
        plt.ion()
        plt.cla()
        plt.axis("equal")
        plt.axis("off")

        triang = tri.Triangulation(dx, dy, el)

        if self._pname != "":
            z = np.zeros(len(nodes))
            if "solution" in self._pname:
                comp = self._pname.split("[")[1].split("]")[0]
                assert comp in types, "Invalid DOF type: %s" % comp

                for n in range(len(nodes)):
                    z[n] = self._solution[dofs.get_dof(n, comp)]
            else:
                name = self._pname.split("[")[0]
                comp = self._pname.split("[")[1].split("]")[0]
                self._write_table(name, globdat)
                table = globdat[gn.TABLES][name]
                assert comp in table, "Invalid component: %s" % comp

                for n in range(len(nodes)):
                    z[n] = table[comp][n]

            plt.tricontourf(
                triang, z, levels=np.linspace(z.min(), z.max(), self._ncolors)
            )

            ticks = np.linspace(z.min(), z.max(), 5, endpoint=True)
            plt.colorbar(ticks=ticks)

        plt.triplot(triang, "k-", linewidth=self._linewidth)
        plt.show(block=False)

        return "ok"

    def shutdown(self, globdat):
        pass

    def _write_table(self, name, globdat):
        nodecount = len(globdat[gn.NSET])
        models = globdat[gn.MODELS]

        if gn.TABLES not in globdat:
            globdat[gn.TABLES] = {}

        params = {}
        params[pn.TABLE] = Table(size=nodecount)
        params[pn.TABLENAME] = name
        params[pn.TABLEWEIGHTS] = np.zeros(nodecount)
        params[pn.SOLUTION] = self._solution

        for model in self.get_relevant_models("GETTABLE", models):
            model.GETTABLE(params, globdat)

        to_xtable(params[pn.TABLE])

        for jcol in range(params[pn.TABLE].column_count()):
            values = params[pn.TABLE].get_col_values(None, jcol)
            params[pn.TABLE].set_col_values(
                None, jcol, values / params[pn.TABLEWEIGHTS]
            )

        params[pn.TABLE].to_table()
        globdat[gn.TABLES][name] = params[pn.TABLE]
