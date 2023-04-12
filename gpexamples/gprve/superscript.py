import sys
sys.path.append('../../')

import numpy as np
from jive.app import main
import jive.util.proputils as pu
from quickviewer import QuickViewer
from copy import deepcopy

props = pu.parse_file('gprve.pro')

props_c = {}
props_c['init'] = deepcopy(props['gpinit'])
props_c['init']['type'] = 'Init'
props_c['solver'] = deepcopy(props['gpsolver'])
props_c['solver']['type'] = 'Linsolve'
props_c['model'] = deepcopy(props['model'])
props_c['model']['models'] = '[ solid, diri ]'
props_c['init']['mesh']['file'] = 'meshes-void/nfib-16_r0.msh'

globdat_c = main.jive(props_c)
u_coarse = globdat_c['state0']
strain_xx_c = globdat_c['tables']['strain']['xx']
strain_yy_c = globdat_c['tables']['strain']['yy']
strain_c = np.append(strain_xx_c, strain_yy_c)

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
f = globdat['extForce']
u = globdat['state0']

strain_xx = globdat['tables']['strain']['xx']
strain_yy = globdat['tables']['strain']['yy']
strain = np.append(strain_xx, strain_yy)

loglikelihood = globdat['gp']['logLikelihood']

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

samples_eps_prior = samples['prior']['strain']
std_eps_xx_prior = np.std(samples_eps_prior['xx'], axis=1)
std_eps_yy_prior = np.std(samples_eps_prior['yy'], axis=1)

samples_eps_post = samples['posterior']['strain']
std_eps_xx_post = np.std(samples_eps_post['xx'], axis=1)
std_eps_yy_post = np.std(samples_eps_post['yy'], axis=1)

eps_xx_post = np.mean(samples_eps_post['xx'], axis=1)
eps_yy_post = np.mean(samples_eps_post['yy'], axis=1)

Phi = globdat['Phi']

err = abs(u - Phi @ u_coarse)
err_grad = abs(strain - Phi @ strain_c)

QuickViewer(u, globdat, comp=0, title=r'Displacement field ($u_x$)')
QuickViewer(u, globdat, comp=1, title=r'Displacement field ($u_y$)')
QuickViewer(strain_xx, globdat, title=r'Strain field ($\varepsilon_{xx}$)')
QuickViewer(strain_yy, globdat, title=r'Strain field ($\varepsilon_{yy}$)')

QuickViewer(err, globdat, comp=0, title=r'Discretization error ($u_x^f - \Phi u_x^c$)')
QuickViewer(err, globdat, comp=1, title=r'Discretization error ($u_x^f - \Phi u_x^c$)')
QuickViewer(err_grad, globdat, comp=0, title=r'Discretization error ($\varepsilon^f_{xx} - \Phi \varepsilon^c_{xx}$)')
QuickViewer(err_grad, globdat, comp=1, title=r'Discretization error ($\varepsilon^f_{yy} - \Phi \varepsilon^c_{yy}$)')

QuickViewer(std_u_prior, globdat, comp=0, title=r'Prior standard deviation ($\sigma_{u_x}$)')
QuickViewer(std_u_prior, globdat, comp=1, title=r'Prior standard deviation ($\sigma_{u_y}$)')
QuickViewer(std_u_post, globdat, comp=0, title=r'Posterior standard deviation ($\sigma_{u_x}$)')
QuickViewer(std_u_post, globdat, comp=1, title=r'Posterior standard deviation ($\sigma_{u_y}$)')
QuickViewer(std_eps_xx_prior, globdat, title=r'Prior standard deviation ($\sigma_{\varepsilon_{xx}}$)')
QuickViewer(std_eps_yy_prior, globdat, title=r'Prior standard deviation ($\sigma_{\varepsilon_{yy}}$)')
QuickViewer(std_eps_xx_post, globdat, title=r'Posterior standard deviation ($\sigma_{\varepsilon_{xx}}$)')
QuickViewer(std_eps_yy_post, globdat, title=r'Posterior standard deviation ($\sigma_{\varepsilon_{yy}}$)')
