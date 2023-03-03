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

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
f = globdat['extForce']
u = globdat['state0']

mean = globdat['gp']['mean']
u_prior = mean['prior']['state0']
f_prior = mean['prior']['extForce']
u_post = mean['posterior']['state0']
f_post = mean['posterior']['extForce']

covariance = globdat['gp']['covariance']
Sigma_prior = covariance['prior']['state0']
Sigma_post = covariance['posterior']['state0']

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

err = abs(u - Phi @ u_coarse)

correction = Sigma_post @ np.linalg.solve(Sigma_prior, u)

QuickViewer(u, globdat, title=r'Exact displacement ($u$)')

QuickViewer(u_post, globdat, title=r'Posterior mean diplacement ($\bar u$)')

fig, (ax1, ax2) = plt.subplots(2, 1, tight_layout=True)
QuickViewer(u - u_post, globdat, ax=ax1, title=r'Discretization error ($u_f - \bar u|$)')
QuickViewer(correction, globdat, ax=ax2, title=r'Correction ($\Sigma^* \Sigma^{-1} u$)')
plt.show()

QuickViewer(std_u_post, globdat, title=r'Posterior standard deviation ($\sqrt{\bar \Sigma_{ii}}$)')

mean_error = []
correction_error = []
total_error= []

sample_count = (2**(np.arange(1, np.log2(10000), step=0.3))).astype(int)
for i in sample_count:
    subsample = samples_u_post[:,:i]

    sample_mean = np.mean(subsample, axis=1)
    diff = subsample - np.tile(sample_mean, (i,1)).T
    sample_covariance = diff @ diff.T / (i-1)

    correction = sample_covariance @ np.linalg.solve(Sigma_prior, u)

    estimate = sample_mean + correction

    mean_error.append(np.sqrt(np.mean((sample_mean - u_post)**2)))
    correction_error.append(np.sqrt(np.mean((correction - (u - u_post))**2)))
    total_error.append(np.sqrt(np.mean((estimate - u)**2)))

expected_mean_error = 1e-4 / np.sqrt(sample_count)

plt.figure()
plt.loglog(sample_count, mean_error, label='mean sample error')
plt.loglog(sample_count, correction_error, label='correction sample error')
plt.loglog(sample_count, total_error, label='total error')
plt.loglog(sample_count, expected_mean_error, label=r'$O(1/\sqrt{n})$')
plt.legend()
plt.show()
