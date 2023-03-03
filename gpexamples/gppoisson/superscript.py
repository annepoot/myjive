import sys
sys.path.append('../../')

from jive.app import main
import jive.util.proputils as pu
from quickviewer import QuickViewer
from copy import deepcopy

props = pu.parse_file('poisson.pro')

props_f = {}
props_f['init'] = deepcopy(props['gpinit'])
props_f['init']['type'] = 'Init'
props_f['solver'] = deepcopy(props['gpsolver'])
props_f['solver']['type'] = 'Linsolve'
props_f['model'] = deepcopy(props['model'])
props_f['model']['models'] = props['model']['models'].replace('gp,', '')

globdat_f = main.jive(props_f)
u_fine = globdat_f['state0']

props_c = deepcopy(props_f)
props_c['init']['mesh']['file'] = 'tri3mesh.msh'
props_c['model']['poisson']['shape']['type'] = 'Triangle3'
props_c['model']['load']['shape']['type'] = 'Triangle3'

globdat_c = main.jive(props_c)
u_coarse = globdat_c['state0']

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
u = globdat['state0']

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

Phi = globdat['Phi']

err = abs(u - Phi @ globdat_c['state0'])

QuickViewer(u_post, globdat, title=r'Posterior mean diplacement ($\bar u$)')

QuickViewer(globdat_c['state0'], globdat_c, title=r'Coarse solution ($u_c$)')

QuickViewer(u, globdat, title=r'Fine solution ($u_f$)')

QuickViewer(err, globdat, title=r'Discretization error ($|u_f - u_c|$)')

QuickViewer(std_u_post, globdat, title=r'Posterior standard deviation ($\sqrt{\bar \Sigma_{ii}}$)')
