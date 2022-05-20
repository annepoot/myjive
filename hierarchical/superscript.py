import sys
sys.path.append('../')

from math import exp
import matplotlib.pyplot as plt
import numpy as np
import main
import proputils as pu
import testutils as tu
from quickviewer import QuickViewer

props = pu.parse_file('beam.pro')

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
u = globdat['state0']
# phi = globdat['phi']
# phi_sub = globdat['phi_sub']
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

phi = globdat['phi']

c = globdat['constraints']
cdofs, cvals = c.get_constraints()
u_post[cdofs] = 0
samples_u_prior[cdofs,:] = 0
samples_u_post[cdofs,:] = 0
std_u_post[cdofs] = 0

err = abs(u_post - u)

QuickViewer(u_post, globdat, title=r'Posterior mean diplacement ($\bar u$)')

QuickViewer(u, globdat, 1, title=r'Exact displacement ($u$)')

QuickViewer(err, globdat, title=r'Posterior mean error ($|\bar u - u|$)')

QuickViewer(std_u_post, globdat, title=r'Posterior standard deviation ($\sqrt{\bar \Sigma_{ii}}$)')
