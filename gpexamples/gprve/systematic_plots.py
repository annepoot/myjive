import sys
sys.path.append('../../')

import numpy as np
from jive.app import main
import jive.util.proputils as pu
from quickviewer import QuickViewer

props = pu.parse_file('rve.pro')

beta_rel = np.logspace(-4,0,5)
epsilon_rel = np.logspace(-3,1,5)

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

        samples = globdat['gp']['samples']

        samples_eps_prior = samples['prior']['strain']
        std_eps_xx_prior = np.std(samples_eps_prior['xx'], axis=1)
        std_eps_yy_prior = np.std(samples_eps_prior['yy'], axis=1)

        samples_eps_post = samples['posterior']['strain']
        std_eps_xx_post = np.std(samples_eps_post['xx'], axis=1)
        std_eps_yy_post = np.std(samples_eps_post['yy'], axis=1)

        fname = 'img/systematic-plots/alpha-{:.0e}_beta-{:.0e}_epsilon-{:.0e}_'.format(alpha, beta, epsilon)

        QuickViewer(std_eps_xx_prior, globdat, pdf=True, title=r'Prior standard deviation ($\sigma_{\varepsilon_{xx}}$)', fname=fname+'std_strain-xx_prior.pdf')
        QuickViewer(std_eps_yy_prior, globdat, pdf=True, title=r'Prior standard deviation ($\sigma_{\varepsilon_{yy}}$)', fname=fname+'std_strain-yy_prior.pdf')
        QuickViewer(std_eps_xx_post, globdat, pdf=True, title=r'Posterior standard deviation ($\sigma_{\varepsilon_{xx}}$)', fname=fname+'std_strain-xx_posterior.pdf')
        QuickViewer(std_eps_xx_post, globdat, pdf=True, title=r'Posterior standard deviation ($\sigma_{\varepsilon_{yy}}$)', fname=fname+'std_strain-yy_posterior.pdf')
