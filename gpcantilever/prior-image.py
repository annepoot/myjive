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

def PrettyViewer(array, globdat, **kwargs):
    maxcolor = np.max(globdat['samples_u_prior'])

    mean = kwargs.get('mean', False)

    defaults = {
        'mincolor': np.min(globdat['samples_u_prior']),
        'maxcolor': np.max(globdat['samples_u_prior']),
        'colorbar': None,
        'colormap': 'BuPu',
        'scale': 1.0,
        'inset': True,
    }

    if mean:
        defaults['boundarywidth'] = 0.5
    else:
        defaults['boundarywidth'] = 0.35
        defaults['alpha'] = 0.1
        defaults['linealpha'] = 0.8

    defaults.update(kwargs)

    QuickViewer(array, globdat, **defaults)

fig, ax = plt.subplots(figsize=(5.2,4), tight_layout = True)
ax.set_axis_off()

for i, sample in enumerate(samples_u_prior.T):
    if np.isclose(i / samples_u_prior.shape[1], 0.8):
        PrettyViewer(u_prior, globdat, ax=ax, mean=True)
    else:
        PrettyViewer(sample, globdat, ax=ax)
PrettyViewer(u_prior, globdat, ax=ax, mean=True, alpha=0.0, linealpha=1.0, boundarywidth=0.75)

folder = '/home/anne/Storage/owncloud/phd/images/cantilever-prior/'

plt.savefig(folder + 'prior-image.pdf', transparent=True)
for dpi in [300, 600, 1200]:
    plt.savefig(folder + 'prior-image-{}dpi.png'.format(dpi), dpi=dpi, transparent=True)

plt.show()
