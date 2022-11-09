import sys
sys.path.append('../../')

import numpy as np
from time import perf_counter
import matplotlib.pyplot as plt
import gputils as gpu
import jive.util.proputils as pu
from jive.app import main

props = pu.parse_file('beam.pro')

densities = ['coarse', 'medium', 'fine', 'fine2', 'fine3', 'fine4']

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

    if N < 1000:
        K_full = K.toarray()
        time1 = perf_counter()
        L1, n_ops1 = gpu.cholesky(K_full, get_nops=True)
        time1 = perf_counter() - time1
        nops1_list.append(n_ops1)
        time1_list.append(time1)

        time2 = perf_counter()
        L2, n_ops2 = gpu.sparse_cholesky(K, get_nops=True)
        time2 = perf_counter() - time2
        nops2_list.append(n_ops2)
        time2_list.append(time2)

    time3 = perf_counter()
    L3, n_ops3 = gpu.incomplete_cholesky(K, get_nops=True)
    time3 = perf_counter() - time3
    nops3_list.append(n_ops3)
    time3_list.append(time3)

    print(N, n_ops1, n_ops2, n_ops3, time1, time2, time3)

N_range = np.array(N_range)
O_n_fit = N_range * 10**-3

plt.figure()
plt.loglog(N_range[:len(time1_list)], time1_list)
plt.loglog(N_range[:len(time2_list)], time2_list)
plt.loglog(N_range[:len(time3_list)], time3_list)
plt.loglog(N_range, O_n_fit)
plt.show()
