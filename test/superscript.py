import sys
sys.path.append('../')

from math import exp
import matplotlib.pyplot as plt
import numpy as np
import main
import proputils as pu
import testutils as tu

props = pu.parse_file('test.pro')

P = 1
L = 10
EA = float(props['model']['bar']['EA'])
k = float(props['model']['bar']['k'])

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
u = globdat['state0']
phi = globdat['phi']
phi_sub = globdat['phi_sub']
