import sys
sys.path.append('../../')

import numpy as np
import matplotlib.pyplot as plt
from sksparse import cholmod as cm
from time import perf_counter

import jive.util.proputils as pu
from jive.app import main

import gputils as gpu

props = pu.parse_file('beam.pro')

densities = ['coarse', 'medium', 'fine', 'fine2', 'fine3', 'fine4', 'fine5', 'fine6']

N_range = []
nops1_list = []
nops2_list = []
nops3_list = []
time1_list = []
time2_list = []
time3_list = []

for density in densities:
    props['init']['mesh']['file'] = 'meshes/beam_' + density + '.msh'
    globdat = main.jive(props)

    K = globdat['matrix0']
    f = globdat['extForce']
    c = globdat['constraints']
    K, f = c.constrain(K, f)

    N = K.shape[0]
    N_range.append(N)

    time1 = perf_counter()
    chol = cm.cholesky(K, ordering_method='amd')
    time1 = perf_counter() - time1
    time1_list.append(time1)

    time2 = perf_counter()
    chol = cm.analyze(K, ordering_method='amd')
    time2 = perf_counter() - time2
    time2_list.append(time2)

    time3 = perf_counter()
    chol = chol.cholesky_inplace(K)
    time3 = perf_counter() - time3
    time3_list.append(time3)

    print(N, time1, time2, time3)

N_range = np.array(N_range)
O_n_fit = N_range * 10**-3

plt.figure()
plt.loglog(N_range[:len(time1_list)], time1_list)
plt.loglog(N_range[:len(time2_list)], time2_list)
plt.loglog(N_range[:len(time3_list)], time3_list)
plt.loglog(N_range, O_n_fit)
plt.show()
