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

props = pu.parse_file('poisson.pro')

props_f = {}
props_f['init'] = deepcopy(props['init'])
props_f['solver'] = deepcopy(props['solver'])
props_f['femodel'] = deepcopy(props['femodel'])

globdat_f = main.jive(props_f)
u_fine = globdat_f['state0']

props_c = deepcopy(props_f)
props_c['init']['mesh']['file'] = 'tri3mesh.msh'
props_c['femodel']['poisson']['shape']['type'] = 'Triangle3'

globdat_c = main.jive(props_c)
u_coarse = globdat_c['state0']

QuickViewer(u_coarse, globdat_c, comp=0, title=r'Coarse displacement field ($u_c$)')
QuickViewer(u_fine, globdat_f, comp=0, title=r'Fine displacement field ($u_f$)')

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

QuickViewer(u_post, globdat, 0, title=r'Posterior mean diplacement ($\bar u$)')

QuickViewer(globdat_c['state0'], globdat_c, 0, title=r'Coarse solution ($u_c$)')

QuickViewer(Phi @ globdat_c['state0'], globdat, 0, title=r'Projected coarse solution ($u_c$)')

QuickViewer(u, globdat, 0, title=r'Exact displacement ($u$)')

QuickViewer(err, globdat, 0, title=r'Discretization error ($|u_f - u_c|$)')

QuickViewer(std_u_post, globdat, 0, title=r'Posterior standard deviation ($\sqrt{\bar \Sigma_{ii}}$)')

for i, sample in enumerate(samples_u_prior.T):

    QuickViewer(sample, globdat, 0, title=r'Prior samples from $u$ (sample {})'.format(i+1))

for i, sample in enumerate(samples_u_post.T):

    QuickViewer(sample, globdat, 0, title=r'Posterior samples from $u$ (sample {})'.format(i+1))
