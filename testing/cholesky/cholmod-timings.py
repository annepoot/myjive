import sys
sys.path.append('../../')

import numpy as np
import matplotlib.pyplot as plt
from sksparse import cholmod as cm
from time import perf_counter

import jive.util.proputils as pu
from jive.app import main

props = pu.parse_file('beam.pro')

densities = ['coarse', 'medium', 'fine', 'fine2', 'fine3', 'fine4', 'fine5', 'fine6', 'fine7', 'fine8']

N_range = []
time1_list = []
time2_list = []
time3_list = []
time4_list = []
time5_list = []
time6_list = []

with open('runtimes.csv', 'w') as file:
    file.write('N,simplicial_total,simplicial_analyze,simplicial_inplace,supernodal_total,supernodal_analyze,supernodal_inplace\n')

for density in densities:
    props['init']['mesh']['file'] = 'meshes/beam_' + density + '.msh'
    globdat = main.jive(props)

    K = globdat['matrix0']
    f = globdat['extForce']
    c = globdat['constraints']
    K, f = c.constrain(K, f)

    K = K.tocsc()

    N = K.shape[0]
    N_range.append(N)

    time1 = perf_counter()
    chol = cm.cholesky(K, mode='simplicial', ordering_method='amd')
    time1 = perf_counter() - time1
    time1_list.append(time1)

    time2 = perf_counter()
    chol = cm.analyze(K, mode='simplicial', ordering_method='amd')
    time2 = perf_counter() - time2
    time2_list.append(time2)

    time3 = perf_counter()
    chol = chol.cholesky_inplace(K)
    time3 = perf_counter() - time3
    time3_list.append(time3)

    time4 = perf_counter()
    chol = cm.cholesky(K, mode='supernodal', ordering_method='amd')
    time4 = perf_counter() - time4
    time4_list.append(time4)

    time5 = perf_counter()
    chol = cm.analyze(K, mode='supernodal', ordering_method='amd')
    time5 = perf_counter() - time5
    time5_list.append(time5)

    time6 = perf_counter()
    chol = chol.cholesky_inplace(K)
    time6 = perf_counter() - time6
    time6_list.append(time6)

    with open('runtimes.csv', 'a') as file:
        file.write('{},{},{},{},{},{},{}\n'.format(N, time1, time2, time3, time4, time5, time6))

    print(N, time1, time2, time3)

N_range = np.array(N_range)
O_simplical_fit = N_range**1.6 * 10**-8
O_supernodal_fit = N_range**1 * 10**-4

plt.figure()
plt.loglog(N_range[:len(time1_list)], time1_list)
plt.loglog(N_range[:len(time2_list)], time2_list)
plt.loglog(N_range[:len(time3_list)], time3_list)
plt.loglog(N_range[:len(time4_list)], time4_list)
plt.loglog(N_range[:len(time5_list)], time5_list)
plt.loglog(N_range[:len(time6_list)], time6_list)
plt.loglog(N_range, O_simplical_fit)
plt.loglog(N_range, O_supernodal_fit)
plt.show()
