import numpy as np
import matplotlib.pyplot as plt

from module import *
from names import GlobNames as gn

LINEWIDTH = 'linewidth'
DEFORM = 'deform'

class FrameViewModule(Module):

    def init(self, props, globdat):
        self._scale = 0.
        self._linewidth = 0.5
        self._shape = props['model']['frame']['shape']['type']
        if self._name in props:
            myprops = props.get(self._name)
            if LINEWIDTH in myprops:
                self._linewidth = float(myprops[LINEWIDTH])
            if DEFORM in myprops:
                self._scale = float(myprops[DEFORM])

    def run(self, globdat):
        nodes = globdat[gn.NSET]
        elems = globdat[gn.ESET]
        disp = globdat[gn.STATE0]
        dofs = globdat[gn.DOFSPACE]
        types = dofs.get_types()

        dx = np.zeros(len(nodes))
        dy = np.zeros(len(nodes))

        for n, node in enumerate(nodes):
            coords = node.get_coords()

            dx[n] = coords[0] + self._scale * disp[dofs.get_dof(n,'dx')]
            dy[n] = coords[1] + self._scale * disp[dofs.get_dof(n,'dy')]

        plt.figure()
        plt.ion()
        plt.cla()
        plt.axis('equal')
        plt.axis('off')

        for elem in elems:
            inodes = elem.get_nodes()
            x = dx[inodes]
            y = dy[inodes]
            plt.plot(x, y, 'k-o', linewidth=self._linewidth)

        plt.show(block=False)

        return 'ok'

    def shutdown(self, globdat):
        pass


def declare(factory):
    factory.declare_module('FrameView', FrameViewModule)
