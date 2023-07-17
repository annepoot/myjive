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

gamma_rel = np.logspace(-4,6,21)
epsilon_rel = np.logspace(-4,4,17)

alphas = []
gammas = []
gammas_rel = []
epsilons = []
epsilons_rel = []
errors = []
errors_rel = []

for i, e_rel in enumerate(epsilon_rel):
    for j, g_rel in enumerate(gamma_rel):
        alpha = 1.
        gamma = alpha * g_rel
        epsilon = alpha * e_rel

        print('Running the model for gamma={:.1e}, epsilon={:.1e}'.format(gamma, epsilon))
        props['model']['gp']['prior']['func'] = 'alpha**2 * M + gamma**2 * F'
        props['model']['gp']['prior']['hyperparams']['alpha'] = str(alpha)
        props['model']['gp']['prior']['hyperparams']['gamma'] = str(gamma)
        props['model']['gp']['obsNoise'] = epsilon

        globdat = main.jive(props)
        u = globdat['state0']
        u_post = globdat['gp']['mean']['posterior']['state0']

        error = np.linalg.norm(u-u_post)
        error_rel = error / np.linalg.norm(u)

        alphas.append(alpha)
        gammas.append(gamma)
        gammas_rel.append(g_rel)
        epsilons.append(epsilon)
        epsilons_rel.append(e_rel)
        errors.append(error)
        errors_rel.append(error_rel)

create_dat(data=[alphas, gammas, gammas_rel, epsilons, epsilons_rel, errors, errors_rel],
           headers=['alpha', 'gamma', 'gamma_rel', 'epsilon', 'epsilon_rel', 'error', 'error_rel'],
           fname='output/contour.dat')
