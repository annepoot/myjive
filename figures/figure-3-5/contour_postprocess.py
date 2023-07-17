import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

data = np.genfromtxt('output/contour.dat', dtype=float, names=True, usecols=('gamma_rel', 'epsilon_rel', 'error_rel'))

gamma_rel = data['gamma_rel']
epsilon_rel = data['epsilon_rel']
error_rel = data['error_rel']

n_gamma = len(np.unique(gamma_rel))
n_epsilon = len(np.unique(epsilon_rel))
shape = (n_epsilon, n_gamma)

Gamma_rel = np.reshape(gamma_rel, shape)
Epsilon_rel = np.reshape(epsilon_rel, shape)
Error_rel = np.reshape(error_rel, shape)

gamma_list   = [ 1e-2, 1e1, 1e3, 1e1, 1e3 ]
epsilon_list = [ 1e-3, 1e0, 1e2, 1e-3, 1e0 ]
label_list   = [ 1, 2, 3, 4, 5 ]

fig,ax=plt.subplots(1,1)
cp = ax.contour(Gamma_rel, Epsilon_rel, Error_rel, levels=[1e-16,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0],norm = LogNorm())
fig.colorbar(cp,format='%.0e')
ax.set_title(r'Phase diagram of the effects of $\alpha$, $\gamma$ and $\varepsilon$ on the posterior')
ax.set_xlabel(r'$\gamma / \alpha$')
ax.set_ylabel(r'$\varepsilon / \alpha$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim((1e-4,1e4))
ax.set_ylim((1e-4,1e4))
ax.plot(gamma_list, epsilon_list, 'k.')
for i, label in enumerate(label_list):
    ax.annotate(label, (gamma_list[i], epsilon_list[i]))
plt.savefig('output/contour.pdf')
plt.show()
