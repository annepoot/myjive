import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
import numpy as np
from jive.app import main
import jive.util.proputils as pu
from quickviewer import QuickViewer
from copy import deepcopy
from plotutils import create_dat

props = pu.parse_file('bridge.pro')

props_c = {}
props_c['init'] = deepcopy(props['gpinit'])
props_c['init']['type'] = 'Init'
props_c['solver'] = deepcopy(props['gpsolver'])
props_c['solver']['type'] = 'Linsolve'
props_c['model'] = deepcopy(props['model'])
props_c['model']['models'] = '[ solid, load, diri ]'
props_c['init']['mesh']['file'] = 'meshes/bridge-q4-r0.msh'

globdat_c = main.jive(props_c)
u_coarse = globdat_c['state0']
eps_xx_c = globdat_c['tables']['strain']['xx']
eps_yy_c = globdat_c['tables']['strain']['yy']
eps_c = np.append(eps_xx_c, eps_yy_c)

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
u = globdat['state0']
eps_xx = globdat['tables']['strain']['xx']
eps_yy = globdat['tables']['strain']['yy']
eps = np.append(eps_xx, eps_yy)

mean = globdat['gp']['mean']
u_prior = mean['prior']['state0']
f_prior = mean['prior']['extForce']
u_post = mean['posterior']['state0']
f_post = mean['posterior']['extForce']

std = globdat['gp']['std']
std_u_prior = std['prior']['state0']
std_f_prior = std['prior']['extForce']
std_u_post = std['posterior']['state0']
std_f_post = std['posterior']['extForce']

samples = globdat['gp']['samples']
samples_u_prior = samples['prior']['state0']
samples_f_prior = samples['prior']['extForce']
samples_u_post = samples['posterior']['state0']
samples_f_post = samples['posterior']['extForce']

samples_eps_xx_prior = samples['prior']['strain']['xx']
samples_eps_yy_prior = samples['prior']['strain']['yy']
samples_eps_xx_post = samples['posterior']['strain']['xx']
samples_eps_yy_post = samples['posterior']['strain']['yy']

eps_xx_prior = np.mean(samples_eps_xx_prior, axis=1)
eps_yy_prior = np.mean(samples_eps_yy_prior, axis=1)
eps_xx_post = np.mean(samples_eps_xx_post, axis=1)
eps_yy_post = np.mean(samples_eps_yy_post, axis=1)

eps_prior = np.append(eps_xx_prior, eps_yy_prior)
eps_post = np.append(eps_xx_post, eps_yy_post)

std_eps_xx_prior = np.std(samples_eps_xx_prior, axis=1)
std_eps_yy_prior = np.std(samples_eps_yy_prior, axis=1)
std_eps_xx_post = np.std(samples_eps_xx_post, axis=1)
std_eps_yy_post = np.std(samples_eps_yy_post, axis=1)

std_eps_prior = np.append(std_eps_xx_prior, std_eps_yy_prior)
std_eps_post = np.append(std_eps_xx_post, std_eps_yy_post)

Phi = globdat['Phi']

err = abs(u - Phi @ u_coarse)
err_grad = abs(eps - Phi @ eps_c)

plt.figure(figsize=(8,2), tight_layout=True)
ax = plt.gca()
QuickViewer(u_post, globdat, ax=ax, colorbar=False, pdf=True, scale=10.0, title=r'Posterior mean displacement ($\mu_u$)')
plt.savefig('img/posterior-mean-state0.pdf')
plt.show()

plt.figure(figsize=(8,2), tight_layout=True)
ax = plt.gca()
QuickViewer(std_u_post, globdat, ax=ax, colorbar=False, pdf=True, title=r'Posterior standard deviation ($\sigma_u$)')
plt.savefig('img/posterior-std-state0.pdf')
plt.show()

plt.figure(figsize=(8,2), tight_layout=True)
ax = plt.gca()
QuickViewer(std_eps_post, globdat, ax=ax, colorbar=False, pdf=True, title=r'Posterior standard deviation ($\sigma_\varepsilon$)')
plt.savefig('img/posterior-std-strain.pdf')
plt.show()

for i, sample in enumerate(samples_u_prior.T[:3]):
    plt.figure(figsize=(8,2), tight_layout=True)
    ax = plt.gca()
    QuickViewer(sample, globdat, ax=ax, colorbar=False, pdf=True, scale=10.0, title=r'Prior samples from $u$ (sample {})'.format(i+1))
    plt.savefig('img/samples/prior-state0-{}.pdf'.format(i+1))
    plt.show()

for i, sample in enumerate(samples_f_prior.T[:3]):
    plt.figure(figsize=(8,2), tight_layout=True)
    ax = plt.gca()
    QuickViewer(sample, globdat, ax=ax, colorbar=False, pdf=True, title=r'Prior samples from $f$ (sample {})'.format(i+1))
    plt.savefig('img/samples/prior-extForce-{}.pdf'.format(i+1))
    plt.show()

for i, sample in enumerate(samples_u_post.T[:3]):
    plt.figure(figsize=(8,2), tight_layout=True)
    ax = plt.gca()
    QuickViewer(sample, globdat, ax=ax, colorbar=False, pdf=True, scale=10.0, title=r'Posterior samples from $u$ (sample {})'.format(i+1))
    plt.savefig('img/samples/posterior-state0-{}.pdf'.format(i+1))
    plt.show()

for i, sample in enumerate(samples_f_post.T[:3]):
    plt.figure(figsize=(8,2), tight_layout=True)
    ax = plt.gca()
    QuickViewer(sample, globdat, ax=ax, colorbar=False, pdf=True, title=r'Posterior samples from $f$ (sample {})'.format(i+1))
    plt.savefig('img/samples/posterior-extForce-{}.pdf'.format(i+1))
    plt.show()

fine_list = ['r0', 'r1', 'r2', 'r3']
x_dict = {}
u_dict = {}

for fineness in fine_list:

    if fineness != 'post':
        pro = deepcopy(props_c)
        pro['init']['mesh']['file'] = 'meshes/bridge-q4-' + fineness + '.msh'
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
        u = globdat['gp']['mean']['posterior']['state0']
        std_u_post = globdat['gp']['std']['posterior']['state0']
        std_u_bottom = []

    x_bottom = []
    u_bottom = []

    for n, node in enumerate(nodes):
        coords = node.get_coords()

        # Check if the node in located on the bottom row
        if np.isclose(coords[1], 0.5):
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
    h_inv = 2**int(fineness[-1])
    if h_inv == 1:
        label = r'h = 1'
    else:
        label = r'$h = \frac{1}{' + str(h_inv) + '}$'
    
    plt.plot(x_dict[fineness], u_dict[fineness], label=label)
    if fineness == 'post':
        u_bar = np.array(u_dict[fineness])
        std_u_bottom = np.array(std_u_bottom)
        std_u_bottom[0] = std_u_bottom[-1] = 0
        plt.fill_between(x_dict[fineness], u_bar - 2*std_u_bottom, u_bar + 2*std_u_bottom, alpha=0.3)

    create_dat([x_dict[fineness],u_dict[fineness]], ['x', 'u_top'], 'results/' + fineness + '-top')
    
plt.legend()
plt.show()
