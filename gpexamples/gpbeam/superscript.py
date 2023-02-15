import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
import numpy as np
from jive.app import main
import jive.util.proputils as pu
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
strain_xx_c = globdat_c['tables']['strain']['xx']
strain_yy_c = globdat_c['tables']['strain']['yy']
strain_c = np.append(strain_xx_c, strain_yy_c)

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
u = globdat['state0']
strain_xx = globdat['tables']['strain']['xx']
strain_yy = globdat['tables']['strain']['yy']
strain = np.append(strain_xx, strain_yy)

f_prior = globdat['f_prior']
u_prior = globdat['u_prior']
f_post = globdat['f_post']
u_post = globdat['u_post']

std_f_prior = globdat['std_f_prior']
std_u_prior = globdat['std_u_prior']
std_f_post = globdat['std_f_post']
std_u_post = globdat['std_u_post']

samples_f_prior = globdat['samples_f_prior']
samples_u_prior = globdat['samples_u_prior']
samples_f_post = globdat['samples_f_post']
samples_u_post = globdat['samples_u_post']

Phi = globdat['Phi']

err = abs(u - Phi @ u_coarse)
err_grad = abs(strain - Phi @ strain_c)

QuickViewer(u_post, globdat, title=r'Posterior mean diplacement ($\bar u$)')

QuickViewer(u, globdat, title=r'Exact displacement ($u$)')

fig, (ax1, ax2) = plt.subplots(2, 1, tight_layout=True)
QuickViewer(err, globdat, ax=ax1, title=r'Discretization error ($|u_f - u_c|$)')
QuickViewer(std_u_post, globdat, ax=ax2, title=r'Posterior standard deviation ($\sqrt{\bar \Sigma_{ii}}$)')
# plt.savefig(fname='img/'+props['init']['mesh']['file'].replace('.msh','').replace('beam_', '')+'.pdf')
plt.show()

fig, (ax1, ax2) = plt.subplots(2, 1, tight_layout=True)
QuickViewer(err_grad, globdat, ax=ax1, comp=0, title=r'Discretization error ($|\varepsilon_f^{xx} - \varepsilon_c^{xx}|$)')
QuickViewer(std_u_post, globdat, ax=ax2, title=r'Posterior standard deviation ($\sqrt{\bar \Sigma_{ii}}$)')
# plt.savefig(fname='img/'+props['init']['mesh']['file'].replace('.msh','').replace('beam_', '')+'.pdf')
plt.show()

fig, (ax1, ax2) = plt.subplots(2, 1, tight_layout=True)
QuickViewer(err_grad, globdat, ax=ax1, comp=1, title=r'Discretization error ($|\varepsilon_f^{yy} - \varepsilon_c^{yy}|$)')
QuickViewer(std_u_post, globdat, ax=ax2, title=r'Posterior standard deviation ($\sqrt{\bar \Sigma_{ii}}$)')
# plt.savefig(fname='img/'+props['init']['mesh']['file'].replace('.msh','').replace('beam_', '')+'.pdf')
plt.show()

QuickViewer(std_u_post, globdat, title=r'Posterior standard deviation ($\sqrt{\bar \Sigma_{ii}}$)')

for i, sample in enumerate(samples_u_prior.T[:3]):

    QuickViewer(sample, globdat, scale=10.0, title=r'Prior samples from $u$ (sample {})'.format(i+1))

for i, sample in enumerate(samples_u_post.T[:3]):

    QuickViewer(sample, globdat, scale=10.0, title=r'Posterior samples from $u$ (sample {})'.format(i+1))

fine_list = ['post', 'coarse', 'medium', 'fine', 'fine2']
x_dict = {}
u_dict = {}

for fineness in fine_list:

    if fineness != 'post':
        pro = deepcopy(props_c)
        pro['init']['mesh']['file'] = 'meshes/beam_' + fineness + '.msh'
        pro['solver']['type'] = 'Linsolve'

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
        std_u_post = np.sqrt(globdat['var_u_post'].diagonal())
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
