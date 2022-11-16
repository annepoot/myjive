import sys
from pathlib import Path
sys.path.append('../../../')

import numpy as np
import scipy.sparse.linalg as spspla
import matplotlib.pyplot as plt
import matplotlib

import jive.util.proputils as pu
from jive.app import main

import solverutils as su
import gputils as gpu
from quickviewer import QuickViewer

matplotlib.use('Agg')

props = pu.parse_file('../beam.pro')

coarse_mesh = '../meshes/beam_medium.msh'
fine_mesh = '../meshes/beam_fine2.msh'

props['init']['mesh']['file'] = coarse_mesh
globdat_c = main.jive(props)
u_c = globdat_c['state0']

gpprops = pu.parse_file('../gpbeam.pro')
gpprops['gpinit']['mesh']['file'] = fine_mesh
gpprops['gpinit']['coarseMesh']['file'] = coarse_mesh
globdat_gp = main.jive(gpprops)
Phi = globdat_gp['Phi']

def get_res_hist(u_hist, uref):
    res_hist = []
    for u in u_hist:
        res = np.linalg.norm(u-uref)
        res_hist.append(res)
    return res_hist

props['init']['mesh']['file'] = fine_mesh
globdat = main.jive(props)

K = globdat['matrix0']
f = globdat['extForce']
c = globdat['constraints']
K, f = c.constrain(K, f)

L = gpu.incomplete_cholesky(K)

Lf = spspla.spsolve(L, f)
LTLf = spspla.spsolve(L.T, Lf)

u0 = Phi @ u_c

preconditioners = [None, 'diag', 'ichol']

for P in preconditioners:
    Pstring = 'none' if P is None else P

    u1_hist = su.preconditioned_conjugate_gradient(K, f, x0=None, P=P, get_history=True)
    u2_hist = su.preconditioned_conjugate_gradient(K, f, x0=u0, P=P, get_history=True)

    uref = spspla.spsolve(K, f)

    res1_hist = get_res_hist(u1_hist, uref)
    res2_hist = get_res_hist(u2_hist, uref)

    modes = ['solutions', 'residuals', 'differences']

    for mode in modes:
        for i, u2 in enumerate(u2_hist):
            fig, (ax1, ax2) = plt.subplots(2, 1)
            ax1.plot(res1_hist, label=r'$u0 = 0$, $P = {}$'.format(Pstring))
            ax1.plot(res2_hist, label=r'$u0 = \Phi u_c$, $P = {}$'.format(Pstring))
            ax1.axvline(i, color='k')
            ax1.set_title('plateau investigation: {}'.format(mode))
            ax1.set_ylabel('residual')
            ax1.legend(loc='upper right')

            if mode == 'solutions':
                QuickViewer(u2, globdat, ax=ax2, scale=10, mincolor=-0.05, maxcolor=0, colormap='viridis_r', tickformat='%.5f')
            elif mode == 'residuals':
                QuickViewer(u2-uref, globdat, ax=ax2, tickformat='%.5f')
            elif mode == 'differences':
                if i < len(u2_hist) - 1:
                    QuickViewer(u2 - u2_hist[i+1], globdat, ax=ax2, tickformat='%.5f')

            path = 'frames/{}/{}'.format(mode, Pstring)
            fname = path+'/frame-{:03d}.png'.format(i)

            print('Generating frame: {}'.format(fname))

            Path(path).mkdir(parents=True, exist_ok=True)

            plt.savefig(fname=fname, dpi=300)
            plt.close(fig)
