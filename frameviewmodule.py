import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from module import *
from names import GlobNames as gn

LINEWIDTH = 'linewidth'
DEFORM = 'deform'
INTERACTIVE = 'interactive'


class FrameViewModule(Module):

    def init(self, props, globdat):
        myprops = props[self._name]
        self._shape = props['model']['frame']['shape']['type']
        self._linewidth = float(myprops.get(LINEWIDTH, 0.5))
        self._scale = float(myprops.get(DEFORM, 0.5))
        self._interactive = bool(eval(myprops.get(INTERACTIVE, 'True')))
        self._step = 0
        self._storeHistory = self._interactive

    def run(self, globdat):
        if self._storeHistory:
            if self._step == 0:
                if gn.HISTORY in globdat:
                    self._storeHistory = False
                else:
                    globdat[gn.HISTORY] = np.array([globdat[gn.STATE0]])
            else:
                globdat[gn.HISTORY] = np.vstack((globdat[gn.HISTORY], globdat[gn.STATE0]))

        self._step += 1
        return 'ok'

    def shutdown(self, globdat):
        if not self._interactive:
            nodes = globdat[gn.NSET]
            elems = globdat[gn.ESET]
            disp = globdat[gn.STATE0]
            dofs = globdat[gn.DOFSPACE]

            dx = np.zeros(len(nodes))
            dy = np.zeros(len(nodes))

            for n, node in enumerate(nodes):
                coords = node.get_coords()
                dx[n] = coords[0] + self._scale * disp[dofs.get_dof(n, 'dx')]
                dy[n] = coords[1] + self._scale * disp[dofs.get_dof(n, 'dy')]

            plt.figure()
            plt.subplots_adjust(left=0.1, bottom=0., right=0.9, top=1.)
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

        elif self._interactive:
            if gn.HISTORY not in globdat:
                raise RuntimeError(gn.HISTORY, 'has not been defined')

            # Initial plot
            nodes = globdat[gn.NSET]
            elems = globdat[gn.ESET]
            disp0 = globdat[gn.HISTORY][0, :]
            dofs = globdat[gn.DOFSPACE]

            dx = np.zeros(len(nodes))
            dy = np.zeros(len(nodes))

            for n, node in enumerate(nodes):
                coords = node.get_coords()
                dx[n] = coords[0] + self._scale * disp0[dofs.get_dof(n, 'dx')]
                dy[n] = coords[1] + self._scale * disp0[dofs.get_dof(n, 'dy')]

            fig, ax = plt.subplots()
            plt.subplots_adjust(left=0.1, bottom=0., right=0.9, top=1.)
            plt.ion()
            plt.cla()
            plt.axis('equal')
            plt.axis('off')

            for i, elem in enumerate(elems):
                inodes = elem.get_nodes()
                if i == 0:
                    x = dx[inodes]
                    y = dy[inodes]
                else:
                    x = np.hstack((x, dx[inodes]))
                    y = np.hstack((y, dy[inodes]))

            line, = plt.plot(x, y, 'k-o', linewidth=self._linewidth)

            # Slider axes
            axcolor = 'lightgoldenrodyellow'
            axstep = plt.axes([0.1, 0.1, 0.8, 0.03], facecolor=axcolor)

            # Create slider object
            s_step = Slider(ax=axstep, label='Step', valmin=1.,
                            valmax=float(len(globdat[gn.HISTORY])), valinit=1., valstep=1.)

            # Slider function. Updates drawing with new x,y coordinates
            def update(val):
                step = int(s_step.val)
                disp_s = globdat[gn.HISTORY][step - 1, :]
                for n_s, node_s in enumerate(nodes):
                    coords_s = node_s.get_coords()
                    dx[n_s] = coords_s[0] + self._scale * disp_s[dofs.get_dof(n_s, 'dx')]
                    dy[n_s] = coords_s[1] + self._scale * disp_s[dofs.get_dof(n_s, 'dy')]

                for j, elem_s in enumerate(elems):
                    inodes_s = elem_s.get_nodes()
                    if j == 0:
                        x_s = dx[inodes_s]
                        y_s = dy[inodes_s]
                    else:
                        x_s = np.hstack((x_s, dx[inodes_s]))
                        y_s = np.hstack((y_s, dy[inodes_s]))

                line.set_xdata(x_s)
                line.set_ydata(y_s)
                fig.canvas.draw_idle()

            s_step.on_changed(update)
            plt.show()

            if gn.SLIDERS not in globdat:
                gn.SLIDERS = []

            gn.SLIDERS.append(s_step)

def declare(factory):
    factory.declare_module('FrameView', FrameViewModule)
