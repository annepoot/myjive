import sys
sys.path.append('../')

from math import exp
import matplotlib.pyplot as plt
import numpy as np
import main
import proputils as pu
import testutils as tu
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
props_c['init'] = deepcopy(props['init'])
props_c['solver'] = deepcopy(props['solver'])
props_c['femodel'] = deepcopy(props['femodel'])
props_c['init']['mesh']['file'] = '2nodebar_coarse.mesh'

P = 1
L = 10
# EA = float(props['femodel']['bar']['EA'])
k = float(props['femodel']['bar']['k'])

globdat_c = main.jive(props_c)
Kc = globdat_c['matrix0']
uc = globdat_c['state0']

globdat = main.jive(props)
Kf = globdat['matrix0']
Mf = globdat['matrix2']
uf = globdat['state0']

K = Kf
M = Mf
u = uf

xf = np.linspace(0, 10, len(uf))
xc = np.linspace(0, 10, len(uc))

step = int((len(xf)-1) / (len(xc)-1))

u_post = globdat['u_post']

print(uc)
print(uf[::step])
print(u_post[::step])

phi = globdat['phi']
f_prior = globdat['f_prior']
u_prior = globdat['u_prior']
f_post = globdat['f_post']
u_post = globdat['u_post']
std_f_prior = np.sqrt(globdat['var_f_prior'])
std_f_post = np.sqrt(globdat['var_f_post'])
std_u_prior = np.sqrt(globdat['var_u_prior'])
std_u_post = np.sqrt(globdat['var_u_post'])

samples_u_prior = globdat['samples_u_prior']
samples_f_prior = globdat['samples_f_prior']
samples_u_post = globdat['samples_u_post']
samples_f_post = globdat['samples_f_post']

plt.figure()
plt.plot(xc, uc, label='coarse solution')
plt.plot(xf, uf, label='fine solution')
plt.plot(xf, u_post, label='posterior mean')
plt.fill_between(xf, u_post - 2*std_u_post, u_post + 2*std_u_post, alpha=0.3)
plt.legend()
plt.show()
