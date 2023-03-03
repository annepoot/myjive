import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
import numpy as np
from jive.app import main
import jive.util.proputils as pu
from jive.solver.constrainer import Constrainer
from copy import deepcopy

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

step = int((len(xf)-1) / (len(xc)-1))

u_post = globdat['u_post']

print(u_coarse)
print(u[::step])
print(u_post[::step])

Phi = globdat['Phi']
Phic = globdat['Phic']

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
samples_u_prior = samples['prior']['state0']
samples_f_prior = samples['prior']['extForce']
samples_u_post = samples['posterior']['state0']
samples_f_post = samples['posterior']['extForce']

plt.figure()
plt.plot(xf, u_post, label='posterior mean')
plt.plot(xf, u_prior, label='prior mean')
plt.plot(xf, samples_u_post, color='gray', linewidth=0.2)
plt.plot(xf, samples_u_prior, color='gray', linewidth=0.2)
plt.fill_between(xf, u_post - 2*std_u_post, u_post + 2*std_u_post, alpha=0.3)
plt.fill_between(xf, u_prior - 2*std_u_prior, u_prior + 2*std_u_prior, alpha=0.3)
plt.plot(xc, u_coarse, label='coarse solution')
plt.plot(xf, u, label='fine solution')
plt.legend()
plt.show()
