import sys
sys.path.append('../../')

import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spspla
from time import perf_counter
import matplotlib.pyplot as plt
import solverutils as su

# Function to quickly generate a 1D bar stiffness matrix
def get_stiffness(N):

    indices = []
    indptr = [0]
    values = []

    for i in range(N):
        if i == 0:
            indices.append(i)
            values.append(1)
        elif i == N-1:
            indices.append(i)
            values.append(1)
        elif i == 1:
            indices.extend([i, i+1])
            values.extend([2, -1])
        elif i == N-2:
            indices.extend([i-1, i])
            values.extend([-1, 2])
        else:
            indices.extend([i-1, i, i+1])
            values.extend([-1, 2, -1])
        indptr.append(len(indices))

    K = spsp.csr_array((values, indices, indptr), shape=(N,N), dtype=float)

    return K

# Function to quickly generate a force vector
def get_force(N):

    f = np.zeros(N, dtype=float)

    if N % 2 == 0:
        f[N//2] = f[N//2-1] = 5
    else:
        f[(N-1)//2] = 10

    return f

N_range = 2**np.arange(0, 20)
time1_list = []
time2_list = []
time3_list = []
time4_list = []

for N in N_range:
    K = get_stiffness(N)
    f_full = get_force(N)

    f = spsp.csr_array(f_full)
    fT = f.T

    time1 = perf_counter()
    for _ in range(10):
        u = spspla.spsolve(K, f_full)
    time1 = perf_counter() - time1
    time1_list.append(time1)

    time2 = perf_counter()
    for _ in range(10):
        u = spspla.spsolve(K, fT)
    time2 = perf_counter() - time2
    time2_list.append(time2)

    time3 = perf_counter()
    for _ in range(10):
        u = spspla.spsolve(K, f.T)
    time3 = perf_counter() - time3
    time3_list.append(time3)

    time4 = perf_counter()
    for _ in range(10):
        u = spspla.spsolve(K, f.transpose())
    time4 = perf_counter() - time4
    time4_list.append(time4)

    print(N, time1, time2, time3, time4)

plt.figure()
plt.loglog(N_range[:len(time1_list)], time1_list, label=r'$K_{s}^{-1} f_{f}$')
plt.loglog(N_range[:len(time2_list)], time2_list, label=r'$K_{s}^{-1} f_{s}$')
plt.loglog(N_range[:len(time3_list)], time3_list, label=r'$K_{s}^{-1} f_{s}.T$')
plt.loglog(N_range[:len(time4_list)], time4_list, label=r'$K_{s}^{-1} f_{s}.tranpose()$')
plt.legend()
plt.show()
