import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

data = np.genfromtxt('output/contour.dat', dtype=float, names=True, usecols=('beta_rel', 'epsilon_rel', 'error_rel'))

beta_rel = data['beta_rel']
epsilon_rel = data['epsilon_rel']
error_rel = data['error_rel']

n_beta = len(np.unique(beta_rel))
n_epsilon = len(np.unique(epsilon_rel))
shape = (n_epsilon, n_beta)

Beta_rel = np.reshape(beta_rel, shape)
Epsilon_rel = np.reshape(epsilon_rel, shape)
Error_rel = np.reshape(error_rel, shape)

beta_interest    = [1e-1,1e3 ,1e5 ,1e2 ,1e0, 1e4]
epsilon_interest = [1e-3,1e-3,1e-2,1e-1,1e1, 1e1]
label_interest = [1,2,3,4,5,6]

fig,ax=plt.subplots(1,1)
cp = ax.contour(Epsilon_rel, Beta_rel, Error_rel, levels=[1e-10,1e-8,1e-6,1e-4,1e-2,1e-1],norm = LogNorm())
fig.colorbar(cp,format='%.0e') # Add a colorbar to a plot
ax.set_title(r'Phase diagram of the effects of $\alpha$, $\beta$ and $\varepsilon$ on the posterior')
ax.set_xlabel(r'$\varepsilon / \alpha$')
ax.set_ylabel(r'$\beta / \alpha$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim((1e-4,1e2))
ax.set_ylim((1e-4,1e4))
# ax.plot(epsilon_interest, beta_interest, 'k.')
# for i, label in enumerate(label_interest):
#     ax.annotate(label, (epsilon_interest[i], beta_interest[i]))
plt.savefig('output/contour.pdf')
plt.show()
