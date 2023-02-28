import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
import numpy as np
from jive.app import main
import jive.util.proputils as pu
from jive.solver.constrainer import Constrainer
from quickviewer import QuickViewer
from copy import deepcopy

props = pu.parse_file('beam.pro')

props_c = {}
props_c['init'] = deepcopy(props['gpinit'])
props_c['init']['type'] = 'Init'
props_c['solver'] = deepcopy(props['gpsolver'])
props_c['solver']['type'] = 'Linsolve'
props_c['model'] = deepcopy(props['model'])
props_c['model']['models'] = '[ solid, load, diri ]'
props_c['init']['mesh']['file'] = 'meshes/beam_coarse.msh'

globdat_c = main.jive(props_c)
u_coarse = globdat_c['state0']

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
u = globdat['state0']
f = globdat['extForce']
c = globdat['constraints']

conman = Constrainer(c, K)
Kc = conman.get_output_matrix()
fc = conman.get_rhs(f)

Phi = globdat['Phi']
Phic = globdat['Phic']

f_prior = globdat['f_prior']
u_prior = globdat['u_prior']
f_post = globdat['f_post']
u_post = globdat['u_post']

Sigma_prior = globdat['var_u_prior']
Sigma_post = globdat['var_u_post']

std_f_prior = globdat['std_f_prior']
std_u_prior = globdat['std_u_prior']
std_f_post = globdat['std_f_post']
std_u_post = globdat['std_u_post']

samples_f_prior = globdat['samples_f_prior']
samples_u_prior = globdat['samples_u_prior']
samples_f_post = globdat['samples_f_post']
samples_u_post = globdat['samples_u_post']

contraction = Sigma_post @ np.linalg.solve(Sigma_prior, u)

QuickViewer(u, globdat, title=r'$u_f$')

QuickViewer(u_post, globdat, title=r'$u^*$')

QuickViewer(u - u_post, globdat, title=r'$u_f - u^*$')

QuickViewer(contraction, globdat, title=r'$\Sigma^* \Sigma^{-1} u_f$')

def get_bottom_values(u, globdat):
    dofs = globdat['dofSpace']
    nodes = globdat['nodeSet']

    x_bottom = []
    u_bottom = []

    y_min = nodes[0].get_coords()[1]

    for node in nodes:
        y = node.get_coords()[1]
        if y < y_min:
            y_min = y

    for n, node in enumerate(nodes):
        coords = node.get_coords()

        # Check if the node in located on the bottom row
        if np.isclose(coords[1], y_min):
            x_bottom.append(coords[0])
            u_bottom.append(u[dofs.get_dof(n, 'dy')])

    # Sort the x_bottom and u_bottom simultaneously
    x_bottom, u_bottom = [list(v) for v in zip(*sorted(zip(x_bottom, u_bottom)))]

    return np.array(x_bottom), np.array(u_bottom)

xf, uf = get_bottom_values(u, globdat)
xf, mf = get_bottom_values(u_post, globdat)
xc, uc = get_bottom_values(u_coarse, globdat_c)
xf, cf = get_bottom_values(contraction, globdat)

fig, (ax1, ax2) = plt.subplots(nrows = 2, figsize=(6,6), tight_layout=True)
ax1.plot(xc, uc, label='coarse solution')
ax1.plot(xc, uc, label='coarse solution')
ax1.plot(xf, uf, label='fine solution')
ax1.legend(loc='lower center')
ax2.plot(xf, uf - mf, label=r'$u_f - u^*$')
ax2.plot(xf, cf, label=r'$\Sigma^* \Sigma^{-1} u$')
ax2.legend(loc='lower center')
ax1.set_title('Vertical displacement at bottom beam')
ax2.set_title('Discretization error vs posterior contraction')
# plt.savefig('img/beam-contraction.pdf')
plt.show()
