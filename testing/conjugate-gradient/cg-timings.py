import sys
sys.path.append('../../')

import numpy as np
from time import perf_counter
import matplotlib.pyplot as plt
import solverutils as su
import gputils as gpu
import jive.util.proputils as pu

from jive.app import main
from jive.solver.constrainer import Constrainer

props = pu.parse_file('beam.pro')

densities = ['coarse', 'medium', 'fine', 'fine2', 'fine3', 'fine4']

N_range = []
nops1_list = []
nops2_list = []
nops3_list = []
time1_list = []
time2_list = []
time3_list = []
time4_list = []

for density in densities:
    props['init']['mesh']['file'] = 'meshes/beam_' + density + '.msh'
    globdat = main.jive(props)

    K = globdat['matrix0']
    f = globdat['extForce']
    c = globdat['constraints']

    conman = Constrainer(c, K)
    K = conman.get_output_matrix()
    f = conman.get_rhs(f)

    L = gpu.incomplete_cholesky(K)

    N = K.shape[0]
    N_range.append(N)

    print('time1')
    time1 = perf_counter()
    u1 = su.conjugate_gradient(K, f)
    time1 = perf_counter() - time1
    time1_list.append(time1)

    print('time2')
    time2 = perf_counter()
    u2 = su.preconditioned_conjugate_gradient(K, f, P=None)
    time2 = perf_counter() - time2
    time2_list.append(time2)

    print('time3')
    time3 = perf_counter()
    u3 = su.preconditioned_conjugate_gradient(K, f, P='diag')
    time3 = perf_counter() - time3
    time3_list.append(time3)

    if N < 5000:
        print('time4')
        time4 = perf_counter()
        u4 = su.preconditioned_conjugate_gradient(K, f, P='ichol', L=L)
        time4 = perf_counter() - time4
        time4_list.append(time4)

    print(N, time1, time2, time3, time4)

N_range = np.array(N_range)

plt.figure()
plt.loglog(N_range[:len(time1_list)], time1_list)
plt.loglog(N_range[:len(time2_list)], time2_list)
plt.loglog(N_range[:len(time3_list)], time3_list)
plt.loglog(N_range[:len(time4_list)], time4_list)
plt.show()
