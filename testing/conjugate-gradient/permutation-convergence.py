import sys
sys.path.append('../../')

import numpy as np
import scipy.sparse.linalg as spspla
import matplotlib.pyplot as plt
from time import perf_counter

from jive.app import main
import jive.util.proputils as pu
from jive.solver.constrainer import Constrainer

import solverutils as su
import gputils as gpu

props = pu.parse_file('beam.pro')

coarse_mesh = 'meshes/beam_medium.msh'
fine_mesh = 'meshes/beam_fine3.msh'

props['init']['mesh']['file'] = coarse_mesh
globdat_c = main.jive(props)
u_c = globdat_c['state0']

gpprops = pu.parse_file('gpbeam.pro')
gpprops['gpinit']['mesh']['file'] = fine_mesh
gpprops['gpinit']['coarseMesh']['file'] = coarse_mesh
globdat_gp = main.jive(gpprops)
Phi = globdat_gp['Phi']

def get_res_hist(u_hist, uref):
    res_hist = []
    for u in u_hist:
        res = np.linalg.norm(u-uref)
        res_hist.append(res)
    return res_hist

props['init']['mesh']['file'] = fine_mesh
globdat = main.jive(props)

K = globdat['matrix0']
f = globdat['extForce']
c = globdat['constraints']

conman = Constrainer(c, K)
K = conman.get_output_matrix()
f = conman.get_rhs(f)

L = gpu.incomplete_cholesky(K)

reorder = su.get_reorder(K)
P = su.get_reorder_matrix(reorder)
K_r = su.reorder_matrix(K, P)
f_r = su.reorder_vector(f, P)

plt.figure()
for u0_string in ['None', 'Phi @ u_c']:

    u0 = eval(u0_string)

    if u0 is None:
        u0_r = None
    else:
        u0_r = su.reorder_vector(u0, P)

    N = K.shape[0]

    print('time1')
    time1 = perf_counter()
    u1_hist = su.preconditioned_conjugate_gradient(K, f, x0=u0, P=None, get_history=True)
    time1 = perf_counter() - time1

    print('time2')
    time2 = perf_counter()
    u2_hist = su.preconditioned_conjugate_gradient(K, f, x0=u0, P='diag', get_history=True)
    time2 = perf_counter() - time2

    if N < 5000:
        print('time3')
        time3 = perf_counter()
        u3_hist = su.preconditioned_conjugate_gradient(K, f, x0=u0, P='ichol', get_history=True)
        time3 = perf_counter() - time3

        print('time4')
        time4 = perf_counter()
        u4_r_hist = su.preconditioned_conjugate_gradient(K_r, f_r, x0=u0_r, P='ichol', get_history=True)
        time4 = perf_counter() - time4

        u4_hist = []

        for u4_r in u4_r_hist:
            u4 = su.rev_reorder_vector(u4_r, P)
            u4_hist.append(u4)

    uref = spspla.spsolve(K, f)

    print(N, time1, time2, time3, time4)

    res1_hist = get_res_hist(u1_hist, uref)
    res2_hist = get_res_hist(u2_hist, uref)
    res3_hist = get_res_hist(u3_hist, uref)
    res4_hist = get_res_hist(u4_hist, uref)

    plt.plot(res1_hist, label=r'$u0 = {}$, $P = None$'.format(u0_string))
    plt.plot(res2_hist, label=r'$u0 = {}$, $P = diag$'.format(u0_string))
    plt.plot(res3_hist, label=r'$u0 = {}$, $P = ichol$'.format(u0_string))
    plt.plot(res4_hist, label=r'$u0 = {}$, $P = ichol\,(amd)$'.format(u0_string))
plt.legend()
plt.xlabel('iterations')
plt.ylabel('residual norm')
plt.title('convergence study')
plt.figure()
