import sys
sys.path.append('../')

from math import exp
import matplotlib.pyplot as plt
import numpy as np
import main
import proputils as pu
import testutils as tu

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
var_f_prior = np.sqrt(globdat['sigma_f_prior'])
var_f_post = np.sqrt(globdat['sigma_f_post'])
var_u_prior = np.sqrt(globdat['sigma_u_prior'])
var_u_post = np.sqrt(globdat['sigma_u_post'])

samples_u_prior = globdat['samples_u_prior']
samples_f_prior = globdat['samples_f_prior']
samples_u_post = globdat['samples_u_post']
samples_f_post = globdat['samples_f_post']

phi = globdat['phi']

x = range(len(u_post))

fix, [[ax1,ax2],[ax3,ax4]] = plt.subplots(2,2)
ax1.fill_between(x, f_prior - 2*var_f_prior, f_prior + 2*var_f_prior, alpha=0.3)
ax1.plot(x, samples_f_prior, color='gray', linewidth=0.2)
ax1.plot(x, f_prior)
ax2.fill_between(x, u_prior - 2*var_u_prior, u_prior + 2*var_u_prior, alpha=0.3)
ax2.plot(x, samples_u_prior, color='gray', linewidth=0.2)
ax2.plot(x, u_prior)
ax3.fill_between(x, f_post - 2*var_f_post, f_post + 2*var_f_post, alpha=0.3)
ax3.plot(x, samples_f_post, color='gray', linewidth=0.2)
ax3.plot(x[1:-1], (K @ u)[1:-1])
ax3.plot(x, f_post)
ax4.fill_between(x, u_post - 2*var_u_post, u_post + 2*var_u_post, alpha=0.3)
ax4.plot(x, samples_u_post, color='gray', linewidth=0.2)
ax4.plot(x, u)
ax4.plot(x, u_post)
plt.show()
