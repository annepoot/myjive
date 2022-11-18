import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('runtimes.txt', delimiter=',', skip_header=1)

N = data[:,0]
simplicial_total = data[:,1]
simplicial_analyze = data[:,2]
simplicial_inplace = data[:,3]
supernodal_total = data[:,4]
supernodal_analyze = data[:,5]
supernodal_inplace = data[:,6]

O_simplical_fit = N**1.6 * 10**-8
O_supernodal_fit = N**1 * 10**-4

plt.figure()
# plt.loglog(N, simplicial_total, label='simplicial total')
plt.loglog(N, simplicial_analyze, label='simplicial analyze')
plt.loglog(N, simplicial_inplace,label='simplicial inplace')
# plt.loglog(N, supernodal_total, label='supernodal total')
plt.loglog(N, supernodal_analyze, label='supernodal analyze')
plt.loglog(N, supernodal_inplace, label='supernodal inplace')

plt.loglog(N, O_simplical_fit, label=r'simplicial fit: $O(n^{1.6})$')
plt.loglog(N, O_supernodal_fit, label=r'supernodal fit: $O(n)$')
plt.title
plt.legend(bbox_to_anchor=(1, 1))
plt.title('CHOLMOD runtimes')
plt.xlabel(r'dofs')
plt.ylabel(r'runtime')
plt.show()
