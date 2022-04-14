import sys
sys.path.append('../')

from math import exp
import matplotlib.pyplot as plt
import numpy as np
import main
import proputils as pu
import testutils as tu

def mesher_quad(L, n):
    dx = L / n / 2
    with open('3nodebar.mesh', 'w') as fmesh:
        fmesh.write('nodes (ID, x, [y], [z])\n')
        for i in range(2 * n + 1):
            fmesh.write('%d %f\n' % (i, i * dx))
        fmesh.write('elements (node#1, node#2, [node#3, ...])\n')
        for i in range(n):
            fmesh.write('%d %d %d\n' % (2 * i, 2 * i + 1, 2 * i + 2))

props = pu.parse_file('3nodebar.pro')

P = 1
L = 10
EA = float(props['model']['bar']['EA'])
k = float(props['model']['bar']['k'])

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
u = globdat['state0']
phi = globdat['phi']
phi_sub = globdat['phi_sub']
u_post = globdat['u_post']
var_u_post = np.sqrt(globdat['sigma_u_post'])

x = range(len(u_post))

fix, ax = plt.subplots()
ax.plot(x, u_post)
ax.fill_between(x, u_post + 2*var_u_post, u_post - 2*var_u_post, alpha=0.3)
ax.plot(x, u)
plt.show()