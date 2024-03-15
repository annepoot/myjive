import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.axes_grid.inset_locator import inset_axes

from myjive.names import GlobNames as gn

LINEWIDTH = "linewidth"
PLOT = "plot"
NCOLORS = "ncolors"
DEFORM = "deform"

__all__ = ["QuickViewer"]


def QuickViewer(array, globdat, **kwargs):
    # Get all possible key word arguments
    comp = kwargs.get("comp", None)
    ax = kwargs.get("ax", None)
    inset = bool(kwargs.get("inset", False))
    scale = float(kwargs.get("scale", 0.0))
    linewidth = kwargs.get("linewidth", None)
    alpha = float(kwargs.get("alpha", 1.0))
    boundarywidth = kwargs.get("boundarywidth", None)
    linealpha = float(kwargs.get("linealpha", alpha))
    colorbar = bool(kwargs.get("colorbar", True))
    tickformat = kwargs.get("tickformat", None)
    cmap = plt.get_cmap(kwargs.get("colormap", "viridis"))
    ncolors = int(kwargs.get("ncolors", 100))
    mincolor = kwargs.get("mincolor", None)
    maxcolor = kwargs.get("maxcolor", None)
    title = kwargs.get("title", None)
    fname = kwargs.get("fname", None)
    pdf = bool(kwargs.get("pdf", False))
    dpi = float(kwargs.get("dpi", 300))
    figsize = kwargs.get("figsize", None)
    meshonly = bool(kwargs.get("meshonly", False))

    # Get the necessary info from globdat
    nodes = globdat[gn.NSET]
    elems = globdat[gn.ESET]
    dofs = globdat[gn.DOFSPACE]
    types = dofs.get_types()
    shape = globdat[gn.MESHSHAPE]

    # Set the component to the y-component, if it exists
    if comp is None and len(array) > len(nodes):
        comp = min(len(types) - 1, 1)

    x = np.zeros(len(nodes))
    y = np.zeros(len(nodes))

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

    for n, node in enumerate(nodes):
        coords = node.get_coords()

        x[n] = coords[0]
        y[n] = coords[1]

    for e, elem in enumerate(elems):
        inodes = elem.get_nodes()

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

    dx = np.copy(x)
    dy = np.copy(y)

    if scale > 0:
        for n in range(len(nodes)):
            idofs = dofs.get_dofs([n], types)
            du = array[idofs]

            if len(idofs) == 2:
                dx[n] += scale * du[0]
                dy[n] += scale * du[1]

    # !!! This should be moved to an appropriate module once BoundaryShapes have been implemented
    if boundarywidth is not None:
        topx = []
        topy = []
        bottomx = []
        bottomy = []
        leftx = []
        lefty = []
        rightx = []
        righty = []

        for n in range(len(nodes)):
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

    no_ax = ax is None

    if no_ax:
        plt.figure(figsize=figsize)
        ax = plt.gca()
    else:
        if inset:
            ax_inset = inset_axes(ax, width="100%", height="100%", loc=10)
            ax_inset.sharex(ax)
            ax_inset.sharey(ax)
            ax = ax_inset

    plt.ion()
    ax.cla()
    ax.set_axis_off()
    if inset:
        ax.set_aspect("equal", adjustable="box")
        ax.patch.set_alpha(0.0)
    else:
        ax.set_aspect("equal", adjustable="datalim")

    triang = tri.Triangulation(dx, dy, el)

    if comp is None:
        assert len(array) == len(nodes)
        z = array
    else:
        z = np.zeros(len(nodes))

        for n, node in enumerate(nodes):
            idofs = dofs.get_dofs([n], types)
            z[n] = array[idofs[comp]]

    if mincolor is None:
        mincolor = z.min()
    if maxcolor is None:
        maxcolor = z.max()

    levels = np.linspace(mincolor, maxcolor, ncolors)

    if not meshonly:
        mappable = ax.tricontourf(triang, z, levels=levels, alpha=alpha, cmap=cmap)

    if colorbar:
        ticks = np.linspace(mincolor, maxcolor, 5, endpoint=True)
        cbar = plt.colorbar(mappable, ticks=ticks, ax=ax, format=tickformat)
        cbar.formatter.set_powerlimits((0, 0))

    if linewidth is not None:
        ax.triplot(triang, "k-", linewidth=linewidth, alpha=linealpha)

    if boundarywidth is not None:
        ax.plot(topx, topy, "k-", linewidth=boundarywidth, alpha=linealpha)
        ax.plot(bottomx, bottomy, "k-", linewidth=boundarywidth, alpha=linealpha)
        ax.plot(rightx, righty, "k-", linewidth=boundarywidth, alpha=linealpha)
        ax.plot(leftx, lefty, "k-", linewidth=boundarywidth, alpha=linealpha)

    # Make sure the contour plot is rendered correctly as a pdf
    if pdf:
        for contour in [mappable]:
            for c in contour.collections:
                c.set_edgecolor("face")
                c.set_linewidth(1e-10)

    if not title is None:
        ax.set_title(title)

    if not fname is None:
        dirname = os.path.dirname(fname)
        if len(dirname) > 0:
            os.makedirs(os.path.dirname(fname), exist_ok=True)
        plt.savefig(fname, dpi=dpi, bbox_inches="tight")

    if no_ax:
        plt.show(block=True)
