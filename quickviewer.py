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

def QuickViewer(array, globdat, comp=1, ax=None, linewidth=0.2, scale=0.0, ncolors=100, title=None, fname=None):

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

        if len(idofs) == 2:
          dx[n] += scale * du[0]
          dy[n] += scale * du[1]

    if ax is None:
        fig = plt.figure()
        ax = plt.gca()
    else:
        fig = ax.get_figure()

    plt.ion()
    ax.cla()
    ax.set_axis_off()
    ax.set_aspect('equal', adjustable='datalim')

    triang = tri.Triangulation (dx, dy, el)

    z = np.zeros(len(nodes))

    for n, node in enumerate(nodes):
        idofs = dofs.get_dofs([n], types)
        z[n] = array[idofs[comp]]

    mappable = ax.tricontourf(triang,z,levels=np.linspace(z.min(),z.max(),ncolors))
    ticks = np.linspace(z.min(),z.max(),5,endpoint=True)
    plt.colorbar(mappable, ticks=ticks,ax=ax)
    ax.triplot(triang,'k-',linewidth=linewidth)

    if not title is None:
        ax.set_title(title)

    if not fname is None:
        plt.savefig(fname, dpi=300)

    if ax is None:
        plt.show(block=False)
