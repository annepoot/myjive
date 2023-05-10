import sys
sys.path.append('../../')

import numpy as np
from jive.app import main
import jive.util.proputils as pu
from quickviewer import QuickViewer

props = pu.parse_file('rve.pro')

beta_rel = np.logspace(-6,4,11)
epsilon_rel = np.logspace(-4,4,9)

for e_rel in epsilon_rel:
    for b_rel in beta_rel:
        alpha = 1.
        beta = alpha * b_rel
        epsilon = alpha * e_rel

        print('Running the model for beta={:.1e}, epsilon={:.1e}'.format(beta, epsilon))
        props['model']['gp']['prior']['func'] = 'alpha**2 * M'
        props['model']['gp']['prior']['hyperparams']['alpha'] = str(alpha)
        props['model']['gp']['obsNoise'] = epsilon
        props['model']['gp']['boundary']['covs'] = '[{},{}]'.format(beta,beta)

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

        fname = 'img/systematic-plots/alpha-{:.0e}_beta-{:.0e}_epsilon-{:.0e}_'.format(alpha, beta, epsilon)

        QuickViewer(u_post, globdat, comp=0, pdf=True, colorbar=False, figsize=(6,6), fname=fname+'mean_state0-x_posterior.pdf')
        QuickViewer(u_post, globdat, comp=1, pdf=True, colorbar=False, figsize=(6,6), fname=fname+'mean_state0-y_posterior.pdf')

        QuickViewer(std_u_prior, globdat, comp=0, pdf=True, colorbar=False, figsize=(6,6), fname=fname+'std_state0-x_prior.pdf')
        QuickViewer(std_u_prior, globdat, comp=1, pdf=True, colorbar=False, figsize=(6,6), fname=fname+'std_state0-y_prior.pdf')
        QuickViewer(std_u_post, globdat, comp=0, pdf=True, colorbar=False, figsize=(6,6), fname=fname+'std_state0-x_posterior.pdf')
        QuickViewer(std_u_post, globdat, comp=1, pdf=True, colorbar=False, figsize=(6,6), fname=fname+'std_state0-y_posterior.pdf')

        QuickViewer(std_eps_xx_prior, globdat, pdf=True, colorbar=False, figsize=(6,6), fname=fname+'std_strain-xx_prior.pdf')
        QuickViewer(std_eps_yy_prior, globdat, pdf=True, colorbar=False, figsize=(6,6), fname=fname+'std_strain-yy_prior.pdf')
        QuickViewer(std_eps_xx_post, globdat, pdf=True, colorbar=False, figsize=(6,6), fname=fname+'std_strain-xx_posterior.pdf')
        QuickViewer(std_eps_yy_post, globdat, pdf=True, colorbar=False, figsize=(6,6), fname=fname+'std_strain-yy_posterior.pdf')
