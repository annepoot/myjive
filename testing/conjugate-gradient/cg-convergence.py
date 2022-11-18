import sys
sys.path.append('../../')

import numpy as np
import scipy.sparse.linalg as spspla
from time import perf_counter
import matplotlib.pyplot as plt
import solverutils as su
import gputils as gpu
import jive.util.proputils as pu
from jive.app import main

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

densities = ['coarse', 'medium', 'fine', 'fine2', 'fine3', 'fine4']

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
K, f = c.constrain(K, f)

L = gpu.incomplete_cholesky(K)

Lf = spspla.spsolve(L, f)
LTLf = spspla.spsolve(L.T, Lf)

plt.figure()
for u0_string in ['Phi @ u_c', 'None']:

    u0 = eval(u0_string)

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
        u3_hist = su.preconditioned_conjugate_gradient(K, f, x0=u0, P='ichol', L=L, get_history=True)
        time3 = perf_counter() - time3

    uref = spspla.spsolve(K, f)

    print(N, time1, time2, time3)

    res1_hist = get_res_hist(u1_hist, uref)
    res2_hist = get_res_hist(u2_hist, uref)
    res3_hist = get_res_hist(u3_hist, uref)

    plt.plot(res1_hist, label=r'$u0 = {}$, $P = None$'.format(u0_string))
    plt.plot(res2_hist, label=r'$u0 = {}$, $P = diag$'.format(u0_string))
    plt.plot(res3_hist, label=r'$u0 = {}$, $P = ichol$'.format(u0_string))
plt.legend()
plt.xlabel('iterations')
plt.ylabel('residual norm')
plt.title('convergence study')
plt.figure()
