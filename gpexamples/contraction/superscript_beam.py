import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
import numpy as np
import main
import proputils as pu
from quickviewer import QuickViewer
from copy import deepcopy

props = pu.parse_file('beam.pro')

props_c = {}
props_c['init'] = deepcopy(props['gpinit'])
props_c['init']['type'] = 'Init'
props_c['solver'] = deepcopy(props['gpsolver'])
props_c['solver']['type'] = 'Solver'
props_c['model'] = deepcopy(props['model'])
props_c['model']['models'] = '[ solid, diri ]'
props_c['init']['mesh']['file'] = 'beam_coarse.msh'

globdat_c = main.jive(props_c)
u_coarse = globdat_c['state0']

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
u = globdat['state0']
f = globdat['extForce']
c = globdat['constraints']

Kc, fc = c.constrain(K, f)

Phi = globdat['Phi']
Phic = globdat['Phic']
f_prior = globdat['f_prior']
u_prior = globdat['u_prior']
f_post = globdat['f_post']
u_post = globdat['u_post']

Sigma_f_prior = globdat['var_f_prior']
Sigma_f_post = globdat['var_f_post']
Sigma_u_prior = globdat['var_u_prior']
Sigma_u_post = globdat['var_u_post']

sig_f_prior = Sigma_f_prior.diagonal()
sig_f_post = Sigma_f_post.diagonal()
sig_u_prior = Sigma_u_prior.diagonal()
sig_u_post = Sigma_u_post.diagonal()

std_f_prior = np.sqrt(sig_f_prior)
std_f_post = np.sqrt(sig_f_post)
std_u_prior = np.sqrt(sig_u_prior)
std_u_post = np.sqrt(sig_u_post)

samples_u_prior = globdat['samples_u_prior']
samples_f_prior = globdat['samples_f_prior']
samples_u_post = globdat['samples_u_post']
samples_f_post = globdat['samples_f_post']

def get_error_rel(u, u_coarse, Phi):
    error = Phi @ u_coarse - u
    error_rel = np.zeros(u.size)

    for i in range(len(u)):
        if np.isclose(u[i], 0):
            if np.isclose(error[i], 0):
                error_rel[i] = 0
            else:
                error_rel[i] = np.nan
        else:
            error_rel[i] = error[i] / u[i]

    return abs(error_rel)

error = abs(Phi @ u_coarse - u)
error_rel = get_error_rel(u, u_coarse, Phi)

cont_hadamard = std_u_post / std_u_prior

QuickViewer(error_rel, globdat, mincolor=0, maxcolor=1, title=r'$|u_c - u_f| / u_f')

QuickViewer(cont_hadamard, globdat, mincolor=0, maxcolor=0.1, title=r'$\sigma_{prior}/\sigma_{post}$')

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
xc, uc = get_bottom_values(u_coarse, globdat_c)
xf, ef = get_bottom_values(error_rel, globdat)
xf, cf = get_bottom_values(cont_hadamard, globdat)

fig, (ax1, ax2) = plt.subplots(nrows = 2, figsize=(6,6), tight_layout=True)
# ax1.plot(xf, u_post, label='posterior mean')
# ax1.plot(xf, u_prior, label='prior mean')
# ax1.plot(xf, samples_u_post, color='gray', linewidth=0.2)
# ax1.plot(xf, samples_u_prior, color='gray', linewidth=0.2)
# ax1.fill_between(xf, u_post - 2*std_u_post, u_post + 2*std_u_post, alpha=0.3)
# ax1.fill_between(xf, u_prior - 2*std_u_prior, u_prior + 2*std_u_prior, alpha=0.3)
ax1.plot(xc, uc, label='coarse solution')
ax1.plot(xf, uf, label='fine solution')
ax1.legend(loc='upper center')
ax2.plot(xf, ef, label=r'$(u_c - u_f)/u_f$')
ax2.plot(xf, cf, label=r'$\sigma_{prior} / \sigma_{post}$')
ax2.set_ylim(0, 0.8)
ax2.legend(loc='upper center')
ax1.set_title('Vertical displacement at bottom beam')
ax2.set_title('Relative error and inverse posterior contraction')
plt.savefig('img/beam-contraction.pdf')
plt.show()
