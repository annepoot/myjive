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
props_c['init'] = deepcopy(props['gpinit'])
props_c['init']['type'] = 'Init'
props_c['solver'] = deepcopy(props['gpsolver'])
props_c['solver']['type'] = 'Solver'
props_c['model'] = deepcopy(props['model'])
props_c['model']['models'] = '[ elastic, diri, neum ]'
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

# QuickViewer(u_post, globdat, title=r'Posterior mean diplacement ($\bar u$)')

# QuickViewer(u, globdat, 1, title=r'Exact displacement ($u$)')

# QuickViewer(err, globdat, 1, title=r'Discretization error ($|u_f - u_c|$)')

# QuickViewer(std_u_post, globdat, title=r'Posterior standard deviation ($\sqrt{\bar \Sigma_{ii}}$)')

fig, ax = plt.subplots()
# ax.cla()
ax.set_axis_off()
# ax.set_aspect('equal', adjustable='datalim')

insert_mean = int(samples_u_prior.shape[1] * 0.8)
for i, sample in enumerate(samples_u_prior.T):

    QuickViewer(sample, globdat, 1, ax=ax, alpha=0.05, mincolor=np.min(samples_u_prior), maxcolor=np.max(samples_u_prior), scale=10.)

    if i == insert_mean:
        QuickViewer(u_prior, globdat, 1, ax=ax, line_fac=1, alpha_fac=1, mincolor=np.min(samples_u_prior), maxcolor=np.max(samples_u_prior), scale=10.)

# plt.title(r'Prior samples from $u$')
plt.show()

# for i, sample in enumerate(samples_u_post.T):

#     QuickViewer(sample, globdat, 1, title=r'Posterior samples from $u$ (sample {})'.format(i+1))
