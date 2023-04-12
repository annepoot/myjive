import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
import numpy as np
from jive.app import main
import jive.util.proputils as pu
from quickviewer import QuickViewer
from copy import deepcopy

props = pu.parse_file('singularity.pro')

props_c = {}
props_c['init'] = deepcopy(props['gpinit'])
props_c['init']['type'] = 'Init'
props_c['solver'] = deepcopy(props['gpsolver'])
props_c['solver']['type'] = 'Linsolve'
props_c['model'] = deepcopy(props['model'])
props_c['model']['models'] = '[ solid, diri, neum ]'
props_c['init']['mesh']['file'] = 'meshes/singularity-r1.msh'

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

QuickViewer(u, globdat, comp=0, pdf=True, title=r'Displacement field ($u_x$)', fname='img/core-plots/state0-x.pdf')
QuickViewer(u, globdat, comp=1, pdf=True, title=r'Displacement field ($u_y$)', fname='img/core-plots/state0-y.pdf')
QuickViewer(strain_xx, globdat, pdf=True, title=r'Strain field ($\varepsilon_{xx}$)', fname='img/core-plots/strain-xx.pdf')
QuickViewer(strain_yy, globdat, pdf=True, title=r'Strain field ($\varepsilon_{yy}$)', fname='img/core-plots/strain-yy.pdf')

QuickViewer(err, globdat, comp=0, pdf=True, title=r'Discretization error ($u_x^f - \Phi u_x^c$)', fname='img/core-plots/error_state0-x.pdf')
QuickViewer(err, globdat, comp=1, pdf=True, title=r'Discretization error ($u_x^f - \Phi u_x^c$)', fname='img/core-plots/error_state0-y.pdf')
QuickViewer(err_grad, globdat, comp=0, pdf=True, title=r'Discretization error ($\varepsilon^f_{xx} - \Phi \varepsilon^c_{xx}$)', fname='img/core-plots/error_strain-xx.pdf')
QuickViewer(err_grad, globdat, comp=1, pdf=True, title=r'Discretization error ($\varepsilon^f_{yy} - \Phi \varepsilon^c_{yy}$)', fname='img/core-plots/error_strain-yy.pdf')

QuickViewer(std_u_prior, globdat, comp=0, pdf=True, title=r'Prior standard deviation ($\sigma_{u_x}$)', fname='img/std_state0-x_prior.pdf')
QuickViewer(std_u_prior, globdat, comp=1, pdf=True, title=r'Prior standard deviation ($\sigma_{u_y}$)', fname='img/std_state0-y_prior.pdf')
QuickViewer(std_u_post, globdat, comp=0, pdf=True, title=r'Posterior standard deviation ($\sigma_{u_x}$)', fname='img/std_state0-x_posterior.pdf')
QuickViewer(std_u_post, globdat, comp=1, pdf=True, title=r'Posterior standard deviation ($\sigma_{u_y}$)', fname='img/std_state0-y_posterior.pdf')
QuickViewer(std_eps_xx_prior, globdat, pdf=True, title=r'Prior standard deviation ($\sigma_{\varepsilon_{xx}}$)', fname='img/std_strain-xx_prior.pdf')
QuickViewer(std_eps_yy_prior, globdat, pdf=True, title=r'Prior standard deviation ($\sigma_{\varepsilon_{yy}}$)', fname='img/std_strain-yy_prior.pdf')
QuickViewer(std_eps_xx_post, globdat, pdf=True, title=r'Posterior standard deviation ($\sigma_{\varepsilon_{xx}}$)', fname='img/std_strain-xx_posterior.pdf')
QuickViewer(std_eps_yy_post, globdat, pdf=True, title=r'Posterior standard deviation ($\sigma_{\varepsilon_{yy}}$)', fname='img/std_strain-yy_posterior.pdf')
