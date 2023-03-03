import sys
sys.path.append('../../')

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

from jive.app import main
import jive.util.proputils as pu

from plotutils import create_dat

# Function to generate 1D meshes
def mesher(n, L=1, fname='bar'):
    dx = L / n
    if not '.' in fname:
        fname += '.mesh'

    with open(fname, 'w') as fmesh:
        fmesh.write('nodes (ID, x, [y], [z])\n')
        for i in range(n + 1):
            fmesh.write('%d %f\n' % (i, i * dx))
        fmesh.write('elements (node#1, node#2, [node#3, ...])\n')
        for i in range(n):
            fmesh.write('%d %d\n' % (i, i + 1))

# Load the properties of the run
props = pu.parse_file('tapered.pro')

# Do a coarse run first, to obtain the coarse solution
props_c = {}
props_c['gpinit'] = deepcopy(props['gpinit'])
props_c['gpsolver'] = deepcopy(props['gpsolver'])
props_c['model'] = deepcopy(props['model'])
props_c['gpinit']['mesh']['file'] = 'bar_coarse.mesh'

# Save the coarse solution
globdat_c = main.jive(props_c)
u_coarse = globdat_c['state0']

# Now do the actual run
globdat = main.jive(props)

# Save the fine and projected coarse displacements
u = globdat['state0']
Phi = globdat['Phi']
u_coarse = Phi @ u_coarse

# Get the prior and posterior means and standard deviations
u_prior     = globdat['gp']['mean']['prior']['state0']
u_post      = globdat['gp']['mean']['posterior']['state0']
std_u_prior = globdat['gp']['std']['prior']['state0']
std_u_post  = globdat['gp']['std']['posterior']['state0']

# Get the prior and posterior samples
samples_u_prior = globdat['gp']['samples']['prior']['state0']
samples_u_post = globdat['gp']['samples']['posterior']['state0']

# Use a fine linspace for plotting
x = np.linspace(0, 1, len(u))

# Create figure 1 directly using matplotlib
plt.figure()
plt.plot(x, u_post, label='posterior mean')
plt.plot(x, u_prior, label='prior mean')
plt.plot(x, samples_u_post, color='gray', linewidth=0.2)
plt.plot(x, samples_u_prior, color='gray', linewidth=0.2)
plt.fill_between(x, u_post - 2*std_u_post, u_post + 2*std_u_post, alpha=0.3)
plt.fill_between(x, u_prior - 2*std_u_prior, u_prior + 2*std_u_prior, alpha=0.3)
plt.plot(x, u_coarse, label='coarse solution')
plt.plot(x, u, label='fine solution')
plt.legend()
plt.show()

# Create output files for latex
create_dat(data=x,
           headers='x',
           fname='output/mesh.dat')

create_dat(data=[u_prior,u_post,std_u_prior,std_u_post,u_coarse,u],
           headers=['u_prior','u_posterior','std_u_prior','std_u_posterior','u_coarse','u_fine'],
           fname='output/results.dat')

create_dat(data=samples_u_prior,
           headers='prior_sample_{}',
           fname='output/samples_prior.dat')

create_dat(data=samples_u_post,
           headers='posterior_sample_{}',
           fname='output/samples_posterior.dat')
