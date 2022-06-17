import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.axes_grid.inset_locator import inset_axes

from module import *
from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

LINEWIDTH ='linewidth'
PLOT = 'plot'
NCOLORS = 'ncolors'
DEFORM = 'deform'
minimum = 0
maximum = 0

def QuickViewer(array, globdat, comp=1, ax=None, linewidth=0.2, scale=0.0, alpha=1., alpha_fac=3., line_fac=3., ncolors=100, mincolor=None, maxcolor=None, title=None, fname=None):

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

    topx = []
    topy = []
    bottomx = []
    bottomy = []
    leftx = []
    lefty = []
    rightx = []
    righty = []

    for n in range(len(nodes)):
        idofs = dofs.get_dofs([n], types)
        du = array[idofs]

        if len(idofs) == 2:
          dx[n] += scale * du[0]
          dy[n] += scale * du[1]

        if np.isclose(y[n], np.max(y)):
            topx.append(dx[n])
            topy.append(dy[n])
        if np.isclose(y[n], np.min(y)):
            bottomx.append(dx[n])
            bottomy.append(dy[n])
        if np.isclose(x[n], np.max(x)):
            rightx.append(dx[n])
            righty.append(dy[n])
        if np.isclose(x[n], np.min(x)):
            leftx.append(dx[n])
            lefty.append(dy[n])

    topx, topy = (list(t) for t in zip(*sorted(zip(topx, topy))))
    bottomx, bottomy = (list(t) for t in zip(*sorted(zip(bottomx, bottomy))))
    rightx, righty = (list(t) for t in zip(*sorted(zip(rightx, righty))))
    leftx, lefty = (list(t) for t in zip(*sorted(zip(leftx, lefty))))

    if ax is None:
        fig = plt.figure()
        ax = plt.gca()
    else:
        fig = ax.get_figure()
        ax_inset = inset_axes(ax, width='100%', height='100%', loc=10)
        ax_inset.sharex(ax)
        ax_inset.sharey(ax)
        ax = ax_inset

    plt.ion()
    ax.cla()
    ax.set_axis_off()
    ax.set_aspect('equal', adjustable='box')
    ax.patch.set_alpha(0.00)

    triang = tri.Triangulation (dx, dy, el)

    z = np.zeros(len(nodes))

    for n, node in enumerate(nodes):
        idofs = dofs.get_dofs([n], types)
        z[n] = array[idofs[comp]]

    # alpha=0.05
    if mincolor is None:
        mincolor = z.min()
    if maxcolor is None:
        maxcolor = z.max()
    levels = np.linspace(mincolor, maxcolor, ncolors)

    alpha_fac = min(alpha_fac*alpha, 1)
    line_fac = min(line_fac*linewidth, 1)

    mappable = ax.tricontourf(triang,z,levels=levels,alpha=alpha)
    # ticks = np.linspace(z.min(),z.max(),5,endpoint=True)
    # plt.colorbar(mappable, ticks=ticks,ax=ax)
    ax.triplot(triang,'k-',linewidth=linewidth, alpha=alpha_fac)
    ax.plot(topx, topy, 'k-', linewidth=line_fac, alpha=alpha_fac)
    ax.plot(bottomx, bottomy, 'k-', linewidth=line_fac, alpha=alpha_fac)
    ax.plot(rightx, righty, 'k-', linewidth=line_fac, alpha=alpha_fac)
    ax.plot(leftx, lefty, 'k-', linewidth=line_fac, alpha=alpha_fac)

    # if not fname is None:
    #     if fname[-4:] == '.pdf':
    # Make sure the contour plot is rendered correctly as a pdf
    for contour in [mappable]:
        for c in contour.collections:
            c.set_edgecolor("face")
            c.set_linewidth(0)

    if not title is None:
        ax.set_title(title)

    if not fname is None:
        plt.savefig(fname, dpi=300)

    if ax is None:
        plt.show(block=False)
