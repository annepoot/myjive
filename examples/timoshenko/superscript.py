import sys, os
cwd = os.getcwd()
rootdir = os.path.join(cwd[:cwd.rfind(os.path.sep + "myjive")], "myjive")
if rootdir not in sys.path:
    sys.path.append(rootdir)

import matplotlib.pyplot as plt
import numpy as np
from jive.app import main
import jive.util.proputils as pu

def mesher(L,n):
    dx = L/n
    with open('timoshenko.mesh','w') as fmesh:
        fmesh.write('nodes (ID, x, [y], [z])\n')
        for i in range(n+1):
            fmesh.write('%d %f\n' %(i,i*dx))
        fmesh.write('elements (node#1, node#2, [node#3, ...])\n')
        for i in range(n):
            fmesh.write('%d %d\n' %(i,i+1))


props = pu.parse_file('timoshenko.pro')

P = 1
L = 10
EI = float(props['model']['timoshenko']['EI'])
GA = float(props['model']['timoshenko']['GAs'])
exact = P*L**3/3/EI + P*L/GA

ns = [1,2,4,8,16,32]
u1 = np.zeros(len(ns))
u2 = np.zeros(len(ns))

for i in range(len(ns)):
    print('\n\nrunning %d\n\n' %ns[i])
    mesher(L,ns[i])

    props['model']['timoshenko']['shape']['intScheme'] = 'Gauss1'
    globdat = main.jive(props)
    u1[i] = globdat['state0'][-1]
    props['model']['timoshenko']['shape']['intScheme'] = 'Gauss2'
    globdat = main.jive(props)
    u2[i] = globdat['state0'][-1]

plt.figure()
plt.loglog(ns,abs(u1-exact))
plt.loglog(ns,abs(u2-exact))
plt.show()
