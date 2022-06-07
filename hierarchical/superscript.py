import sys
sys.path.append('../')

from math import exp
import matplotlib.pyplot as plt
import numpy as np
import main
import proputils as pu
import testutils as tu
from quickviewer import QuickViewer
from copy import deepcopy

props = pu.parse_file('beam.pro')

props_c = {}
props_c['init'] = deepcopy(props['init'])
props_c['solver'] = deepcopy(props['solver'])
props_c['femodel'] = deepcopy(props['femodel'])
props_c['init']['mesh']['file'] = 'beam_coarse.msh'

globdat_c = main.jive(props_c)
u_coarse = globdat_c['state0']

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
u = globdat['state0']

f_prior = globdat['f_prior']
u_prior = globdat['u_prior']
f_post = globdat['f_post']
u_post = globdat['u_post']
std_f_prior = np.sqrt(globdat['var_f_prior'])
std_f_post = np.sqrt(globdat['var_f_post'])
std_u_prior = np.sqrt(globdat['var_u_prior'])
std_u_post = np.sqrt(globdat['var_u_post'])

samples_u_prior = globdat['samples_u_prior']
samples_f_prior = globdat['samples_f_prior']
samples_u_post = globdat['samples_u_post']
samples_f_post = globdat['samples_f_post']

Phi = globdat['Phi']

err = abs(u - Phi @ globdat_c['state0'])

QuickViewer(u_post, globdat, title=r'Posterior mean diplacement ($\bar u$)')

QuickViewer(u, globdat, 1, title=r'Exact displacement ($u$)')

fig, (ax1, ax2) = plt.subplots(2, 1, tight_layout=True)
QuickViewer(err, globdat, ax=ax1, title=r'Discretization error ($|u_f - u_c|$)')
QuickViewer(std_u_post, globdat, ax=ax2, title=r'Posterior standard deviation ($\sqrt{\bar \Sigma_{ii}}$)')
# plt.savefig(dpi=450, fname='img/'+props['init']['mesh']['file'].replace('.msh','').replace('beam_', ''))
plt.show()

QuickViewer(std_u_post, globdat, title=r'Posterior standard deviation ($\sqrt{\bar \Sigma_{ii}}$)')

for i, sample in enumerate(samples_u_prior.T):

    QuickViewer(sample, globdat, 1, scale=10.0, title=r'Prior samples from $u$ (sample {})'.format(i+1))

for i, sample in enumerate(samples_u_post.T):

    QuickViewer(sample, globdat, 1, scale=10.0, title=r'Posterior samples from $u$ (sample {})'.format(i+1))

fine_list = ['post', 'coarse', 'medium', 'fine', 'fine2']
x_dict = {}
u_dict = {}

for fineness in fine_list:

    if fineness != 'post':
        pro = {}
        pro['init'] = deepcopy(props['init'])
        pro['femodel'] = deepcopy(props['femodel'])
        pro['solver'] = deepcopy(props['solver'])
        pro['init']['mesh']['file'] = 'beam_' + fineness + '.msh'

        glob = main.jive(pro)

        dofs = glob['dofSpace']
        elems = glob['elemSet']
        nodes = glob['nodeSet']
        u = glob['state0']

    else:
        dofs = globdat['dofSpace']
        elems = globdat['elemSet']
        nodes = globdat['nodeSet']
        u = globdat['u_post']
        std_u_post = np.sqrt(globdat['var_u_post'])
        std_u_bottom = []

    x_bottom = []
    u_bottom = []

    for n, node in enumerate(nodes):
        coords = node.get_coords()

        # Check if the node in located on the bottom row
        if np.isclose(coords[1], 0):
            x_bottom.append(coords[0])
            u_bottom.append(u[dofs.get_dof(n, 'dy')])

            if fineness == 'post':
                std_u_bottom.append(std_u_post[dofs.get_dof(n, 'dy')])

    if fineness == 'post':
        x_bottom,u_bottom,std_u_bottom =[list(v) for v in zip(*sorted(zip(x_bottom,u_bottom,std_u_bottom)))]
    else:
        x_bottom,u_bottom =[list(v) for v in zip(*sorted(zip(x_bottom,u_bottom)))]

    x_dict[fineness] = x_bottom
    u_dict[fineness] = u_bottom

plt.figure()

for fineness in fine_list:
    plt.plot(x_dict[fineness], u_dict[fineness], label=fineness)
    if fineness == 'post':
        u_bar = np.array(u_dict[fineness])
        std_u_bottom = np.array(std_u_bottom)
        std_u_bottom[0] = std_u_bottom[-1] = 0
        plt.fill_between(x_dict[fineness], u_bar - 2*std_u_bottom, u_bar + 2*std_u_bottom, alpha=0.3)

plt.legend()
plt.show()
