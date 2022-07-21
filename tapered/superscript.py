import sys
sys.path.append('../')

from math import exp
import matplotlib.pyplot as plt
import numpy as np
import main
import proputils as pu
import testutils as tu
import plotutils as plu
from copy import deepcopy

def mesher_lin(L, n, fname='2nodebar_coarse'):
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

props = pu.parse_file('tapered.pro')

props_c = {}
props_c['gpinit'] = deepcopy(props['gpinit'])
props_c['gpsolver'] = deepcopy(props['gpsolver'])
props_c['model'] = deepcopy(props['model'])
props_c['gpinit']['mesh']['file'] = '2nodebar_coarse.mesh'

globdat_c = main.jive(props_c)
u_coarse = globdat_c['state0']

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
u = globdat['state0']
f = globdat['extForce']
c = globdat['constraints']

Kc, fc = c.constrain(K, f)

xf = np.linspace(0, 10, len(u))
xc = np.linspace(0, 10, len(u_coarse))

step = int((len(xf)-1) / (len(xc)-1))

u_post = globdat['u_post']

print(u_coarse)
print(u[::step])
print(u_post[::step])

Phi = globdat['Phi']
Phic = globdat['Phic']
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

folder = '/home/anne/Storage/owncloud/phd/latex/go-no-go/img/1d-bar/'

plt.figure()
plt.plot(xf, u_post, label='posterior mean')
plt.plot(xf, u_prior, label='prior mean')
plt.plot(xf, samples_u_post, color='gray', linewidth=0.2, alpha=0.3)
plt.plot(xf, samples_u_prior, color='gray', linewidth=0.2, alpha=0.3)
plt.fill_between(xf, u_post - 2*std_u_post, u_post + 2*std_u_post, alpha=0.3)
plt.fill_between(xf, u_prior - 2*std_u_prior, u_prior + 2*std_u_prior, alpha=0.3)
plt.plot(xc, u_coarse, label='coarse solution')
plt.plot(xf, u, label='fine solution')
plt.legend()
# plt.savefig(fname=folder+'1d-bar-results.pdf')
plt.show()

plu.create_dat(data=xf,
                headers='x',
                fname=folder+'mesh')

plu.create_dat(data=[u_prior, u_post, std_u_prior, std_u_post, Phi @ u_coarse, u],
               headers=['u_prior', 'u_posterior', 'std_u_prior', 'std_u_posterior', 'u_coarse', 'u_fine'],
               fname=folder+'results')

plu.create_dat(data=samples_u_prior,
               headers='prior_sample_{}',
               fname=folder+'samples_prior')

plu.create_dat(data=samples_u_post,
               headers='posterior_sample_{}',
               fname=folder+'samples_posterior')
