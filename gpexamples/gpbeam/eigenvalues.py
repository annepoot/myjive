import sys
sys.path.append('../../')

import numpy as np
from jive.app import main
import jive.util.proputils as pu
from quickviewer import QuickViewer
from copy import deepcopy

props = pu.parse_file('beam.pro')
props['model']['gp']['prior']['func'] = 'K'
props['model']['gp']['prior']['hyperparams'] = {}

props_c = {}
props_c['init'] = deepcopy(props['gpinit'])
props_c['init']['type'] = 'Init'
props_c['solver'] = deepcopy(props['gpsolver'])
props_c['solver']['type'] = 'Linsolve'
props_c['model'] = deepcopy(props['model'])
props_c['model']['models'] = '[ solid, diri, load ]'
props_c['init']['mesh']['file'] = 'meshes/beam_coarse.msh'

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

cov = globdat['gp']['covariance']
cov_u_prior = cov['prior']['state0']
cov_f_prior = cov['prior']['extForce']
cov_u_post = cov['posterior']['state0']
cov_f_post = cov['posterior']['extForce']

samples = globdat['gp']['samples']
samples_u_prior = samples['prior']['state0']
samples_f_prior = samples['prior']['extForce']
samples_u_post = samples['posterior']['state0']
samples_f_post = samples['posterior']['extForce']

Phi = globdat['Phi']

err = abs(u - Phi @ u_coarse)
err_grad = abs(strain - Phi @ strain_c)

QuickViewer(u, globdat, comp=1, dpi=600, title=r'displacement ($u_x$)')
QuickViewer(err, globdat, comp=1, dpi=600, title=r'displacement error ($\Delta_{u_x}$)')

QuickViewer(u_post, globdat, comp=1, dpi=600, title=r'posterior mean displacement ($m^*$)')
QuickViewer(std_u_post, globdat, comp=1, dpi=600, title=r'posterior std displacement ($\sigma^*_{u_{x}}$)')
QuickViewer(cov_u_post @ f, globdat, comp=1, dpi=600, title=r'error estimate ($\Sigma^* f$)')

dc = len(u)
pdNoise = 1e-4

cov_u_prior += pdNoise**2 * np.identity(dc)
cov_u_post += pdNoise**2 * np.identity(dc)

l, Q = np.linalg.eigh(cov_u_post)

# newl = l * (Q.T @ f)**2
# newl = (l * (Q.T @ f))**2
newl = l * abs(Q.T @ f)
newcov = Q @ np.diag(newl) @ Q.T
newvar = newcov.diagonal()
newstd = np.sqrt(newvar)

QuickViewer(newvar, globdat, comp=1, dpi=600, title=r'rescaled posterior var ($\sigma_{u_x}$)')
QuickViewer(newstd, globdat, comp=1, dpi=600, title=r'rescaled posterior std ($\sigma^*_{u_x}$)')

nsample = samples_u_post.shape[1]
sample_mean = np.mean(samples_u_post, axis=1)
A = samples_u_post - np.tile(sample_mean, (nsample,1)).T
sample_cov = 1/(nsample-1) * A @ A.T
sample_std = np.sqrt(sample_cov.diagonal())

U, d, VT = np.linalg.svd(A, full_matrices=False)
l = d**2

# newl_sample = l * (U.T @ f)**2
# newl_sample = (l * (U.T @ f))**2
newl_sample = l * abs(U.T @ f)
newcov_sample = 1/(nsample-1) * U @ np.diag(newl_sample) @ U.T
newvar_sample = newcov_sample.diagonal()
newstd_sample = np.sqrt(newvar_sample)

QuickViewer(newstd_sample, globdat, comp=1, dpi=600, title=r'rescaled sample posterior std ($\sigma^*_{u_x}$)')
