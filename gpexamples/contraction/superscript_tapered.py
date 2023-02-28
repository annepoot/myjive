import sys
sys.path.append('../../')

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

from jive.app import main
import jive.util.proputils as pu
from jive.solver.constrainer import Constrainer

def mesher_lin(L, n, fname='2nodebar_coarse'):
    dx = L / n
    if not '.' in fname:
        fname += '.mesh'

    with open(fname, 'w') as fmesh:
        fmesh.write('nodes (ID, x, [y], [z])\n')
        for i in range(n + 1):
            fmesh.write('%d %f\n' % (i, i * dx))
        fmesh.write('elements (node#1, node#2, [node#3, ...])\n')
        for i in range(n):
            fmesh.write('%d %d\n' % (i, i + 1))

props = pu.parse_file('tapered.pro')

props_c = {}
props_c['gpinit'] = deepcopy(props['gpinit'])
props_c['gpsolver'] = deepcopy(props['gpsolver'])
props_c['model'] = deepcopy(props['model'])
props_c['gpinit']['mesh']['file'] = '2nodebar_coarse.mesh'

globdat_c = main.jive(props_c)
u_coarse = globdat_c['state0']

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
u = globdat['state0']
f = globdat['extForce']
c = globdat['constraints']

conman = Constrainer(c, K)
Kc = conman.get_output_matrix()
fc = conman.get_rhs(f)

xf = np.linspace(0, 10, len(u))
xc = np.linspace(0, 10, len(u_coarse))

Phi = globdat['Phi']
Phic = globdat['Phic']

f_prior = globdat['f_prior']
u_prior = globdat['u_prior']
f_post = globdat['f_post']
u_post = globdat['u_post']

Sigma_prior = globdat['var_u_prior']
Sigma_post = globdat['var_u_post']

std_f_prior = globdat['std_f_prior']
std_u_prior = globdat['std_u_prior']
std_f_post = globdat['std_f_post']
std_u_post = globdat['std_u_post']

samples_f_prior = globdat['samples_f_prior']
samples_u_prior = globdat['samples_u_prior']
samples_f_post = globdat['samples_f_post']
samples_u_post = globdat['samples_u_post']

fig, (ax1, ax2) = plt.subplots(nrows = 2, figsize=(6,6), tight_layout=True)
ax1.plot(xf, u_post, label='posterior mean')
ax1.plot(xf, u_prior, label='prior mean')
ax1.plot(xf, samples_u_post, color='gray', linewidth=0.2)
ax1.plot(xf, samples_u_prior, color='gray', linewidth=0.2)
ax1.fill_between(xf, u_post - 2*std_u_post, u_post + 2*std_u_post, alpha=0.3)
ax1.fill_between(xf, u_prior - 2*std_u_prior, u_prior + 2*std_u_prior, alpha=0.3)
ax1.plot(xf, Phi @ u_coarse, label='coarse solution')
ax1.plot(xf, u, label='fine solution')
ax1.legend(loc='upper left')
ax2.plot(xf, u - u_post, label=r'$u_f - u^*$')
ax2.plot(xf, Sigma_post @ np.linalg.solve(Sigma_prior, u - u_prior), label=r'$\Sigma_{post} \Sigma_{prior}^{-1} u_f$')
ax2.legend(loc='upper center')
ax1.set_title('Prior and posterior distribution of the tapered bar problem')
ax2.set_title('Discretization error vs posterior contraction')
plt.show()
