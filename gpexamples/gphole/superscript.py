import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
import numpy as np
from jive.app import main
import jive.util.proputils as pu
from quickviewer import QuickViewer
from copy import deepcopy

props = pu.parse_file('cantilever-hole.pro')

props_c = {}
props_c['init'] = deepcopy(props['gpinit'])
props_c['init']['type'] = 'Init'
props_c['solver'] = deepcopy(props['gpsolver'])
props_c['solver']['type'] = 'Linsolve'
props_c['model'] = deepcopy(props['model'])
props_c['model']['models'] = '[ solid, diri ]'
props_c['init']['mesh']['file'] = 'meshes/q4/hole-9.msh'
props_c['model']['solid']['shape']['type'] = 'Quad4'

globdat_c = main.jive(props_c)
u_c = globdat_c['state0']
strain_xx_c = globdat_c['tables']['strain']['xx']
strain_yy_c = globdat_c['tables']['strain']['yy']
strain_c = np.append(strain_xx_c, strain_yy_c)

globdat = main.jive(props)
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

err = abs(u - Phi @ u_c)
err_grad = abs(strain - Phi @ strain_c)

fig, (ax1, ax2) = plt.subplots(2, 1, tight_layout=True)
QuickViewer(err, globdat, comp=0, ax=ax1, title=r'Discretization error ($x$) ($|u_f - u_c|$)')
QuickViewer(std_u_post, globdat, comp=0, ax=ax2, title=r'Posterior standard deviation ($x$) ($\sqrt{\bar \Sigma_{ii}}$)')
# plt.savefig(fname='img/'+props['init']['mesh']['file'].replace('.msh','').replace('beam_', '')+'.pdf')
plt.show()

fig, (ax1, ax2) = plt.subplots(2, 1, tight_layout=True)
QuickViewer(err, globdat, comp=1, ax=ax1, title=r'Discretization error ($y$) ($|u_f - u_c|$)')
QuickViewer(std_u_post, globdat, comp=1, ax=ax2, title=r'Posterior standard deviation ($y$) ($\sqrt{\bar \Sigma_{ii}}$)')
# plt.savefig(fname='img/'+props['init']['mesh']['file'].replace('.msh','').replace('beam_', '')+'.pdf')
plt.show()
