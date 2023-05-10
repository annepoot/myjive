import sys
sys.path.append('../../')

import numpy as np
from jive.app import main
import jive.util.proputils as pu
from plotutils import create_dat

props = pu.parse_file('contour.pro')

globdat = main.jive(props)
u = globdat['state0']
u_post = globdat['gp']['mean']['posterior']['state0']

beta_rel = np.logspace(-6,4,21)
epsilon_rel = np.logspace(-4,4,17)

alphas = []
betas = []
betas_rel = []
epsilons = []
epsilons_rel = []
errors = []
errors_rel = []

for i, e_rel in enumerate(epsilon_rel):
    for j, b_rel in enumerate(beta_rel):
        alpha = 1.
        beta = alpha * b_rel
        epsilon = alpha * e_rel

        print('Running the model for beta={:.1e}, epsilon={:.1e}'.format(beta, epsilon))
        props['model']['gp']['prior']['func'] = 'alpha**2 * M'
        props['model']['gp']['prior']['hyperparams']['alpha'] = str(alpha)
        props['model']['gp']['boundary']['covs'] = '['+str(beta)+','+str(beta)+']'
        props['model']['gp']['obsNoise'] = epsilon

        globdat = main.jive(props)
        u = globdat['state0']
        u_post = globdat['gp']['mean']['posterior']['state0']

        error = np.linalg.norm(u-u_post)
        error_rel = error / np.linalg.norm(u)

        alphas.append(alpha)
        betas.append(beta)
        betas_rel.append(b_rel)
        epsilons.append(epsilon)
        epsilons_rel.append(e_rel)
        errors.append(error)
        errors_rel.append(error_rel)

create_dat(data=[alphas, betas, betas_rel, epsilons, epsilons_rel, errors, errors_rel],
           headers=['alpha', 'beta', 'beta_rel', 'epsilon', 'epsilon_rel', 'error', 'error_rel'],
           fname='output/contour.dat')
