import sys
sys.path.append('../')

import matplotlib.pyplot as plt
import numpy as np
import main
import proputils as pu
import testutils as tu
from quickviewer import QuickViewer
from copy import deepcopy
from matplotlib.ticker import ScalarFormatter


props = pu.parse_file('beam.pro')

fig, ax = plt.subplots()

# Make data.
alpha = np.logspace(-1, 1, 20)
beta = np.round(np.linspace(0, 0.1, 5), 4)

for b in beta:
    likelihood = np.zeros_like(alpha)

    for i, a in enumerate(alpha):
        props['model']['gp']['prior']['hyperparams']['alpha'] = a
        props['model']['gp']['prior']['hyperparams']['beta'] = b
        globdat = main.jive(props)
        likelihood[i] = globdat['logLikelihood']

    plt.plot(alpha, likelihood, label=r'$\beta = {}$'.format(b))

ax.set_xlabel(r'$\alpha$')
ax.set_xscale('log')
ax.get_xaxis().set_major_formatter(ScalarFormatter())
ax.set_ylabel(r'log likelihood')
plt.legend()
plt.savefig('img/alpha_beta_likelihood.pdf')
plt.show()
