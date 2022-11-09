import sys
sys.path.append('../../')

import numpy as np
import scipy.sparse as spsp
from time import perf_counter
import matplotlib.pyplot as plt
import gputils as gpu

# Function to quickly generate a 1D bar stiffness matrix
def get_stiffness(N):

    K = np.zeros((N,N), dtype=float)
    K_elem = np.array([[1,-1],[-1,1]], dtype=float)

    for i in range(N-1):
        K[i:i+2,i:i+2] += K_elem

    K[ 0, :] = K[ :, 0] = 0
    K[-1, :] = K[ :,-1] = 0
    K[ 0, 0] = K[-1,-1] = 1

    return K

N_range = 2**np.arange(0, 15)
nops1_list = []
nops2_list = []
nops3_list = []
time1_list = []
time2_list = []
time3_list = []

for N in N_range:
    K = get_stiffness(N)

    if N < 1000:
        time1 = perf_counter()
        L, n_ops1 = gpu.cholesky(K, get_nops=True)
        time1 = perf_counter() - time1
        nops1_list.append(n_ops1)
        time1_list.append(time1)

    K = spsp.csr_array(K)

    time2 = perf_counter()
    L, n_ops2 = gpu.sparse_cholesky(K, get_nops=True)
    time2 = perf_counter() - time2
    nops2_list.append(n_ops2)
    time2_list.append(time2)

    time3 = perf_counter()
    L, n_ops3 = gpu.incomplete_cholesky(K, get_nops=True)
    time3 = perf_counter() - time3
    nops3_list.append(n_ops3)
    time3_list.append(time3)

    print(N, n_ops1, n_ops2, n_ops3, time1, time2, time3)

O_n3_fit = N_range**3 * 10**-7
O_n_fit = N_range * 10**-4

plt.figure()
plt.loglog(N_range[:len(time1_list)], time1_list)
plt.loglog(N_range[:len(time2_list)], time2_list)
plt.loglog(N_range[:len(time3_list)], time3_list)
plt.loglog(N_range, O_n3_fit)
plt.loglog(N_range, O_n_fit)
plt.show()
