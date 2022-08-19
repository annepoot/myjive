import sys
sys.path.append('../')

import matplotlib.pyplot as plt
import numpy as np
import main
import proputils as pu
import testutils as tu
from quickviewer import QuickViewer
from copy import deepcopy

props = pu.parse_file('prior-image.pro')

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
u = globdat['state0']

f_prior = globdat['f_prior']
u_prior = globdat['u_prior']
f_post = globdat['f_post']
u_post = globdat['u_post']
std_f_prior = np.sqrt(globdat['var_f_prior'])
std_f_post = np.sqrt(globdat['var_f_post'])
std_u_prior = np.sqrt(globdat['var_u_prior'])
std_u_post = np.sqrt(globdat['var_u_post'])

samples_u_prior = globdat['samples_u_prior']
samples_f_prior = globdat['samples_f_prior']
samples_u_post = globdat['samples_u_post']
samples_f_post = globdat['samples_f_post']

Phi = globdat['Phi']

QuickViewer(u_prior, globdat, scale=-1.0, title=r'Prior mean displacement')
