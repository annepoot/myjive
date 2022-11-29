import sys
sys.path.append('../../')

import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spspla
from numba import njit
from time import perf_counter

import jive.util.proputils as pu
from jive.app import main
from ichol import incomplete_cholesky

@njit
def solve_lower_triangular(data, indices, indptr, b):
    x = np.zeros_like(b)

    for i in range(len(b)):
        cols = indices[indptr[i]:indptr[i+1]]
        vals = data[indptr[i]:indptr[i+1]]
        numerator = b[i]

        for j, L_ij in zip(cols[:-1], vals[:-1]):
            numerator -= L_ij * x[j]

        x[i] = numerator / vals[-1]
    return x


@njit
def solve_upper_triangular(data, indices, indptr, b):
    x = np.zeros_like(b)

    for i in range(len(b)-1, -1, -1):
        cols = indices[indptr[i]:indptr[i+1]]
        vals = data[indptr[i]:indptr[i+1]]
        numerator = b[i]

        for j, L_ij in zip(cols[1:], vals[1:]):
            numerator -= L_ij * x[j]

        x[i] = numerator / vals[0]
    return x


@njit
def solve_triangular(data, indices, indptr, b, lower=True):
    x = np.zeros_like(b)

    if lower:
        rng = range(len(b))
    else:
        rng = range(len(b)-1, -1, -1)

    for i in rng:
        if lower:
            start = indptr[i]
            end = indptr[i+1]-1
            L_ii = data[indptr[i+1]-1]
        else:
            start = indptr[i]+1
            end = indptr[i+1]
            L_ii = data[indptr[i]]

        cols = indices[start:end]
        vals = data[start:end]
        numerator = b[i]

        for j, L_ij in zip(cols, vals):
            numerator -= L_ij * x[j]
        x[i] = numerator / L_ii
    return x

props = pu.parse_file('beam.pro')
props['init']['mesh']['file'] = 'meshes/beam_fine4.msh'
props['solver']['solver'] = 'direct'
globdat = main.jive(props)

K = globdat['matrix0']
f = globdat['extForce']
u = globdat['state0']
L = spsp.tril(K, format='csr')
L.data = incomplete_cholesky(L.data, L.indices, L.indptr)
LT = L.T.tocsr()

time1 = perf_counter()
tmp =  spspla.spsolve(L, f)
u1 =  spspla.spsolve(LT, tmp)
time1 = perf_counter() - time1

time2 = perf_counter()
tmp = spspla.spsolve_triangular(L, f, lower=True)
u2 = spspla.spsolve_triangular(LT, tmp, lower=False)
time2 = perf_counter() - time2

time3 = perf_counter()
tmp = solve_lower_triangular(L.data, L.indices, L.indptr, f)
u3 = solve_upper_triangular(LT.data, LT.indices, LT.indptr, tmp)
time3 = perf_counter() - time3

time4 = perf_counter()
tmp = solve_lower_triangular(L.data, L.indices, L.indptr, f)
u4 = solve_upper_triangular(LT.data, LT.indices, LT.indptr, tmp)
time4 = perf_counter() - time4

time5 = perf_counter()
tmp = solve_triangular(L.data, L.indices, L.indptr, f, lower=True)
u5 = solve_triangular(LT.data, LT.indices, LT.indptr, tmp, lower=False)
time5 = perf_counter() - time5

time6 = perf_counter()
tmp = solve_triangular(L.data, L.indices, L.indptr, f, lower=True)
u6 = solve_triangular(LT.data, LT.indices, LT.indptr, tmp, lower=False)
time6 = perf_counter() - time6

print('Reference time (spsolve):', time1)
print('Reference time (spsolve_triangular):', time2)
print('Time including compilation:', time3)
print('Time excluding compilation:', time4)
print('Time including compilation:', time5)
print('Time excluding compilation:', time6)
print('Speedup factor:', time2 / time6)
