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

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
u = globdat['state0']
# phi = globdat['phi']
# phi_sub = globdat['phi_sub']
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

phi = globdat['phi']

err = abs(u_post - u)

QuickViewer(u_post, globdat, title=r'Posterior mean diplacement ($\bar u$)')

QuickViewer(u, globdat, 1, title=r'Exact displacement ($u$)')

QuickViewer(err, globdat, title=r'Posterior mean error ($|\bar u - u|$)')

QuickViewer(std_u_post, globdat, title=r'Posterior standard deviation ($\sqrt{\bar \Sigma_{ii}}$)')

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
