import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

from module import *
from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

LINEWIDTH ='linewidth'
PLOT = 'plot'
NCOLORS = 'ncolors'
DEFORM = 'deform'

class ViewModule(Module):

    def init(self, props, globdat):
        self._scale = 0.0
        self._linewidth = 0.2
        self._pname = ''
        self._ncolors = 100

        if self._name in props:
            myprops = props.get(self._name)

            if LINEWIDTH in myprops:
                self._linewidth = float(myprops[LINEWIDTH])
            if PLOT in myprops:
                self._pname = myprops[PLOT]
            if DEFORM in myprops:
                self._scale = myprops[DEFORM]
            if NCOLORS in myprops:
                self._ncolors = int(myprops[NCOLORS])

        self._modelname = myprops.get(gn.MODEL, gn.MODEL)

    def run(self, globdat):

        nodes = globdat[gn.NSET]
        elems = globdat[gn.ESET]
        disp = globdat[gn.STATE0]
        dofs = globdat[gn.DOFSPACE]
        types = dofs.get_types()

        assert globdat[gn.MESHSHAPE] == 'Triangle3', 'ViewModule only supports triangles for now'

        x = np.zeros(len(nodes))
        y = np.zeros(len(nodes))
        el= np.zeros((len(elems),3),dtype=int)

        for n, node in enumerate(nodes):
            coords = node.get_coords()

            x[n] = coords[0]
            y[n] = coords[1]

        for e, elem in enumerate(elems):
            inodes = elem.get_nodes()

            el[e,:] = inodes

        dx = np.copy(x)
        dy = np.copy(y)

        for n in range(len(nodes)):
            idofs = dofs.get_dofs([n], types)
            du = disp[idofs]

            dx[n] += self._scale * du[0]
            dy[n] += self._scale * du[1]

        plt.figure()
        ax = plt.gca()
        plt.ion()
        plt.cla()
        plt.axis('equal')
        plt.axis('off')

        triang = tri.Triangulation (dx, dy, el)

        if self._pname != '':
            z = np.zeros(len(nodes))
            if 'solution' in self._pname:
                comp = self._pname.split('[')[1].split(']')[0]
                assert comp in types, 'Invalid DOF type: %s' % comp

                for n in range(len(nodes)):
                    z[n] = disp[dofs.get_dof(n,comp)]
            else:
                name = self._pname.split('[')[0]
                comp = self._pname.split('[')[1].split(']')[0]
                self._write_table(name, globdat)
                table = globdat[gn.TABLES][name]
                assert comp in table, 'Invalid component: %s' % comp

                for n in range(len(nodes)):
                    z[n] = table[comp][n]

            plt.tricontourf(triang,z,levels=np.linspace(z.min(),z.max(),self._ncolors))

            ticks = np.linspace(z.min(),z.max(),5,endpoint=True)
            plt.colorbar(ticks=ticks)

        plt.triplot(triang,'k-',linewidth=self._linewidth)
        plt.show(block=False)

        return 'ok'

    def shutdown(self, globdat):
        pass

    def _write_table(self, name, globdat):
        nodes = globdat[gn.NSET]
        model = globdat[self._modelname]

        globdat[gn.TABLES] = {}

        params = {}
        params[pn.TABLE] = {}
        params[pn.TABLENAME] = name
        params[pn.TABLEWEIGHTS] = np.zeros(len(nodes))

        model.take_action(act.GETTABLE, params, globdat)

        table = params[pn.TABLE]
        tbwts = params[pn.TABLEWEIGHTS]

        for col in table.values():
            col /= params[pn.TABLEWEIGHTS]

        globdat[gn.TABLES][name] = params[pn.TABLE]


def declare(factory):
    factory.declare_module('View', ViewModule)
