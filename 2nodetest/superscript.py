import sys
sys.path.append('../')

from math import exp
import matplotlib.pyplot as plt
import numpy as np
import main
import proputils as pu
import testutils as tu

def mesher_lin(L, n):
    dx = L / n
    with open('2nodebar.mesh', 'w') as fmesh:
        fmesh.write('nodes (ID, x, [y], [z])\n')
        for i in range(n + 1):
            fmesh.write('%d %f\n' % (i, i * dx))
        fmesh.write('elements (node#1, node#2, [node#3, ...])\n')
        for i in range(n):
            fmesh.write('%d %d\n' % (i, i + 1))

props = pu.parse_file('2nodebar.pro')

P = 1
L = 10
EA = float(props['femodel']['bar']['EA'])
k = float(props['femodel']['bar']['k'])

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
u = globdat['state0']
phi = globdat['phi']
phi_sub = globdat['phi_sub']
f_post = globdat['f_post']
u_post = globdat['u_post']
var_f_prior = np.sqrt(globdat['sigma_f_prior'])
var_f_post = np.sqrt(globdat['sigma_f_post'])
var_u_prior = np.sqrt(globdat['sigma_u_prior'])
var_u_post = np.sqrt(globdat['sigma_u_post'])

x = range(len(u_post))

fix, [[ax1,ax2],[ax3,ax4]] = plt.subplots(2,2)
ax1.fill_between(x, -2*var_f_prior, 2*var_f_prior, alpha=0.3)
ax2.fill_between(x, -2*var_u_prior, 2*var_u_prior, alpha=0.3)
ax3.plot(x, f_post)
ax3.fill_between(x, f_post - 2*var_f_post, f_post + 2*var_f_post, alpha=0.3)
ax3.plot(x, K @ u)
ax4.plot(x, u_post)
ax4.fill_between(x, u_post - 2*var_u_post, u_post + 2*var_u_post, alpha=0.3)
ax4.plot(x, u)
plt.show()