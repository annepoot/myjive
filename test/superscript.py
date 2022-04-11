import sys
sys.path.append('../')

from math import exp
import matplotlib.pyplot as plt
import numpy as np
import main
import proputils as pu
import testutils as tu

props = pu.parse_file('bar.pro')

P = 1
L = 10
EA = float(props['model']['bar']['EA'])
k = float(props['model']['bar']['k'])
alpha = np.sqrt(k / EA)
N_L = P*exp(-alpha*L)
energy = alpha*P*P/2/k
energyCorrected = alpha*P*P/2/k - alpha*N_L*N_L/2/k
u_L = -N_L/EA/alpha

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
u = globdat['state0']
