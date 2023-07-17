import sys
sys.path.append('../../')

import numpy as np
from jive.app import main
import jive.util.proputils as pu
from quickviewer import QuickViewer

props = pu.parse_file('singularity.pro')

gamma_list   = [ 1e-2, 1e1, 1e3, 1e1, 1e3 ]
epsilon_list = [ 1e-3, 1e0, 1e2, 1e-3, 1e0 ]

for g_rel, e_rel in zip(gamma_list, epsilon_list):
    alpha = 1.
    gamma = alpha * g_rel
    epsilon = alpha * e_rel

    print('Running the model for gamma={:.1e}, epsilon={:.1e}'.format(gamma, epsilon))
    props['model']['gp']['prior']['func'] = 'alpha**2 * M + gamma**2 * F'
    props['model']['gp']['prior']['hyperparams']['alpha'] = str(alpha)
    props['model']['gp']['prior']['hyperparams']['gamma'] = str(gamma)
    props['model']['gp']['obsNoise'] = epsilon

    globdat = main.jive(props)

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

    samples_eps_prior = samples['prior']['strain']
    std_eps_xx_prior = np.std(samples_eps_prior['xx'], axis=1)
    std_eps_yy_prior = np.std(samples_eps_prior['yy'], axis=1)

    samples_eps_post = samples['posterior']['strain']
    std_eps_xx_post = np.std(samples_eps_post['xx'], axis=1)
    std_eps_yy_post = np.std(samples_eps_post['yy'], axis=1)

    fname = 'img/systematic-plots/alpha-{:.0e}_gamma-{:.0e}_epsilon-{:.0e}_'.format(alpha, gamma, epsilon)

    QuickViewer(u_post, globdat, comp=0, dpi=600, colorbar=False, figsize=(6,6), fname=fname+'mean_state0-x_posterior.png')
    QuickViewer(u_post, globdat, comp=1, dpi=600, colorbar=False, figsize=(6,6), fname=fname+'mean_state0-y_posterior.png')

    QuickViewer(std_u_prior, globdat, comp=0, dpi=600, colorbar=False, figsize=(6,6), fname=fname+'std_state0-x_prior.png')
    QuickViewer(std_u_prior, globdat, comp=1, dpi=600, colorbar=False, figsize=(6,6), fname=fname+'std_state0-y_prior.png')
    QuickViewer(std_u_post, globdat, comp=0, dpi=600, colorbar=False, figsize=(6,6), fname=fname+'std_state0-x_posterior.png')
    QuickViewer(std_u_post, globdat, comp=1, dpi=600, colorbar=False, figsize=(6,6), fname=fname+'std_state0-y_posterior.png')

    QuickViewer(std_eps_xx_prior, globdat, dpi=600, colorbar=False, figsize=(6,6), fname=fname+'std_strain-xx_prior.png')
    QuickViewer(std_eps_yy_prior, globdat, dpi=600, colorbar=False, figsize=(6,6), fname=fname+'std_strain-yy_prior.png')
    QuickViewer(std_eps_xx_post, globdat, dpi=600, colorbar=False, figsize=(6,6), fname=fname+'std_strain-xx_posterior.png')
    QuickViewer(std_eps_yy_post, globdat, dpi=600, colorbar=False, figsize=(6,6), fname=fname+'std_strain-yy_posterior.png')
