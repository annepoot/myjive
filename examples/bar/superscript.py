import sys, os
rootdir = os.path.abspath(os.path.join(os.getcwd(), "..", ".."))
jivedir = os.path.abspath(os.path.join(rootdir, "..", "myjive"))
sys.path.append(rootdir)
sys.path.append(jivedir)

from math import exp
import matplotlib.pyplot as plt
import numpy as np
from jive.app import main
import jive.util.proputils as pu

def mesher_lin(L, n):
    dx = L / n
    with open('bar.mesh', 'w') as fmesh:
        fmesh.write('nodes (ID, x, [y], [z])\n')
        for i in range(n + 1):
            fmesh.write('%d %f\n' % (i, i * dx))
        fmesh.write('elements (node#1, node#2, [node#3, ...])\n')
        for i in range(n):
            fmesh.write('%d %d\n' % (i, i + 1))


def mesher_quad(L, n):
    dx = L / n / 2
    with open('bar.mesh', 'w') as fmesh:
        fmesh.write('nodes (ID, x, [y], [z])\n')
        for i in range(2 * n + 1):
            fmesh.write('%d %f\n' % (i, i * dx))
        fmesh.write('elements (node#1, node#2, [node#3, ...])\n')
        for i in range(n):
            fmesh.write('%d %d %d\n' % (2 * i, 2 * i + 1, 2 * i + 2))


props = pu.parse_file('bar.pro')

P = 1
L = 10
EA = float(props['model']['bar']['EA'])
k = float(props['model']['bar']['k'])
alpha = np.sqrt(k / EA)
N_L = P*exp(-alpha*L)
energy = alpha*P*P/2/k
energyCorrected = alpha*P*P/2/k - alpha*N_L*N_L/2/k
u_L = -N_L/EA/alpha

ns = [4,8,16,32,64];
E1 = np.zeros(len(ns))
E2 = np.zeros(len(ns))

for i in range(len(ns)):
    print('\n\nrunning %d\n\n' % ns[i])

    props['model']['bar']['shape']['type'] = 'Line2'
    props['model']['bar']['shape']['intScheme'] = 'Gauss2'
    props['model']['diri']['values'] = '['+str(u_L) + ']'
    mesher_lin(L, ns[i])
    globdat = main.jive(props)
    K = globdat['matrix0']
    u = globdat['state0']
    E1[i] = 0.5 * u @ K @ u

    props['model']['bar']['shape']['type'] = 'Line3'
    props['model']['bar']['shape']['intScheme'] = 'Gauss3'
    mesher_quad(L, ns[i])
    globdat = main.jive(props)
    K = globdat['matrix0']
    u = globdat['state0']
    E2[i] = 0.5 * u @ K @ u

plt.figure()
ref = energy
plt.loglog(ns, abs(E1 - ref)/ref, label='Linear')
plt.loglog(ns, abs(E2 - ref)/ref, label='Quadratic')
# ref = energyCorrected
# plt.loglog(ns, abs(E1 - ref)/ref, label='Linear corrected')
# plt.loglog(ns, abs(E2 - ref)/ref, label='Quadratic corrected')
plt.xlabel('Number of elements')
plt.ylabel('Error')
plt.legend()
plt.show()
