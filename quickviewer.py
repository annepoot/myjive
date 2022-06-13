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

    shape = globdat[gn.MESHSHAPE]

    assert  shape == 'Triangle3' or shape == 'Triangle6', 'ViewModule only supports triangles for now'

    x = np.zeros(len(nodes))
    y = np.zeros(len(nodes))

    if shape == 'Triangle3':
        el= np.zeros((len(elems),3),dtype=int)
    elif shape == 'Triangle6':
        el= np.zeros((len(elems)*4,3),dtype=int)

    for n, node in enumerate(nodes):
        coords = node.get_coords()

        x[n] = coords[0]
        y[n] = coords[1]

    for e, elem in enumerate(elems):
        inodes = elem.get_nodes()

        if shape == 'Triangle3':
            el[e,:] = inodes
        elif shape == 'Triangle6':
            el[4*e+0,:] = inodes[[0,3,5]]
            el[4*e+1,:] = inodes[[1,4,3]]
            el[4*e+2,:] = inodes[[2,5,4]]
            el[4*e+3,:] = inodes[[3,4,5]]

    dx = np.copy(x)
    dy = np.copy(y)

    for n in range(len(nodes)):
        idofs = dofs.get_dofs([n], types)
        du = array[idofs]

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

    if not fname is None:
        if fname[-4:] == '.pdf':
            # Make sure the contour plot is rendered correctly as a pdf
            for contour in [mappable]:
                for c in contour.collections:
                    c.set_edgecolor("face")
                    c.set_linewidth(0.00000000000000001)

    if not title is None:
        ax.set_title(title)

    if not fname is None:
        plt.savefig(fname, dpi=300)

    if ax is None:
        plt.show(block=False)
