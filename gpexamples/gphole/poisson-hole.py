import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
import numpy as np
from jive.app import main
import jive.util.proputils as pu
from quickviewer import QuickViewer
from copy import deepcopy

props = pu.parse_file('poisson-hole.pro')

props_c = {}
props_c['init'] = deepcopy(props['gpinit'])
props_c['init']['type'] = 'Init'
props_c['solver'] = deepcopy(props['gpsolver'])
props_c['solver']['type'] = 'Linsolve'
props_c['model'] = deepcopy(props['model'])
props_c['model']['models'] = '[ poisson, diri ]'
props_c['init']['mesh']['file'] = 'meshes/q4/hole-9.msh'
props_c['model']['poisson']['shape']['type'] = 'Quad4'

globdat_c = main.jive(props_c)
u_c = globdat_c['state0']

globdat = main.jive(props)
u = globdat['state0']

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

err = abs(u - Phi @ u_c)

fig, (ax1, ax2) = plt.subplots(2, 1, tight_layout=True)
QuickViewer(err, globdat, comp=0, ax=ax1, title=r'Discretization error ($u$) ($|u_f - u_c|$)')
QuickViewer(std_u_post, globdat, comp=0, ax=ax2, title=r'Posterior standard deviation ($u$) ($\sqrt{\bar \Sigma_{ii}}$)')
# plt.savefig(fname='img/'+props['init']['mesh']['file'].replace('.msh','').replace('beam_', '')+'.pdf')
plt.show()
