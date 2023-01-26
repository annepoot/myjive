import sys
sys.path.append('../../')

from jive.app import main
import jive.util.proputils as pu
from quickviewer import QuickViewer
from copy import deepcopy

props = pu.parse_file('poisson.pro')

props_f = {}
props_f['init'] = deepcopy(props['gpinit'])
props_f['init']['type'] = 'Init'
props_f['solver'] = deepcopy(props['gpsolver'])
props_f['solver']['type'] = 'Linsolve'
props_f['model'] = deepcopy(props['model'])
props_f['model']['models'] = props['model']['models'].replace('gp,', '')

globdat_f = main.jive(props_f)
u_fine = globdat_f['state0']

props_c = deepcopy(props_f)
props_c['init']['mesh']['file'] = 'tri3mesh.msh'
props_c['model']['poisson']['shape']['type'] = 'Triangle3'
props_c['model']['load']['shape']['type'] = 'Triangle3'

globdat_c = main.jive(props_c)
u_coarse = globdat_c['state0']

globdat = main.jive(props)
K = globdat['matrix0']
M = globdat['matrix2']
u = globdat['state0']

f_prior = globdat['f_prior']
u_prior = globdat['u_prior']
f_post = globdat['f_post']
u_post = globdat['u_post']

std_f_prior = globdat['std_f_prior']
std_u_prior = globdat['std_u_prior']
std_f_post = globdat['std_f_post']
std_u_post = globdat['std_u_post']

samples_f_prior = globdat['samples_f_prior']
samples_u_prior = globdat['samples_u_prior']
samples_f_post = globdat['samples_f_post']
samples_u_post = globdat['samples_u_post']

Phi = globdat['Phi']

err = abs(u - Phi @ globdat_c['state0'])

QuickViewer(u_post, globdat, title=r'Posterior mean diplacement ($\bar u$)')

QuickViewer(globdat_c['state0'], globdat_c, title=r'Coarse solution ($u_c$)')

QuickViewer(u, globdat, title=r'Fine solution ($u_f$)')

QuickViewer(err, globdat, title=r'Discretization error ($|u_f - u_c|$)')

QuickViewer(std_u_post, globdat, title=r'Posterior standard deviation ($\sqrt{\bar \Sigma_{ii}}$)')

# # Use a direct solver for reference
# props['model']['gp']['solver']['type'] = 'cholmod'
# globdat_ref = main.jive(props)

# data = []

# for preconditioner in ['id', 'diag', 'ichol']:
#     for coarse_init in [True, False]:
#         for max_iter in 2**np.arange(8):
#             print(preconditioner, coarse_init, max_iter)

#             props['model']['gp']['solver']['type'] = 'CG'
#             props['model']['gp']['solver']['allowMaxIter'] = 'True'
#             props['model']['gp']['solver']['maxIter'] = str(max_iter)
#             props['model']['gp']['preconditioner']['type'] = preconditioner
#             props['model']['gp']['coarseInit'] = str(coarse_init)

#             globdat = main.jive(props)

#             sample = globdat['samples_u_post'][:,0]
#             sample_ref = globdat_ref['samples_u_post'][:,0]
#             sample_rmse = np.sqrt(np.sum((sample-sample_ref)**2))

#             std = globdat['std_u_post']
#             std_ref = globdat_ref['std_u_post']
#             std_rmse = np.sqrt(np.sum((std-std_ref)**2))

#             data.append([max_iter, preconditioner, coarse_init, sample_rmse, std_rmse])

#             QuickViewer(globdat['samples_u_post'][:,0], globdat,
#                         fname='img/sample-post/sample-post_iterMax-{}_P-{}_u0-{}.png'.format(max_iter, preconditioner, 'uc' if coarse_init else '0'),
#                         title='Single posterior sample ($i_{{max}} = {}, P = {}, u_0 = {}$)'.format(max_iter, preconditioner, 'u_c' if coarse_init else '0'))

#             QuickViewer(globdat['std_u_post'], globdat,
#                         fname='img/std-post/std-post_iterMax-{}_P-{}_u0-{}.png'.format(max_iter, preconditioner, 'uc' if coarse_init else '0'),
#                         title='Posterior std ($i_{{max}} = {}, P = {}, u_0 = {}$)'.format(max_iter, preconditioner, 'u_c' if coarse_init else '0'))

# df = pd.DataFrame(data, columns=['maxIter', 'preconditioner', 'coarseInit', 'rmse_sample', 'rmse_std'])

# df.to_csv('rmse_data.csv')
