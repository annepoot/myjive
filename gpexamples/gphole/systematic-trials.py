import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
import numpy as np
from jive.app import main
import jive.util.proputils as pu
from quickviewer import QuickViewer
from copy import deepcopy
from plotutils import create_dat

props = pu.parse_file('systematic-trials.pro')

for load in ['forced displacement', 'point load']:
    for direction in ['horizontal', 'vertical', 'poisson']:
        for loc in [1, 3, 5, 7, 9]:
            for prior in ['M', 'K']:
                if direction == 'poisson':
                    model = 'poisson'
                else:
                    model = 'solid'

                props['model']['models'] = '[' + model + ',gp,diri,neum]'

                if load == 'forced displacement':
                    props['model']['neum']['groups'] = '[]'
                    props['model']['neum']['dofs']   = '[]'
                    props['model']['neum']['values'] = '[]'
                    if model == 'solid':
                        props['model']['diri']['groups'] = '[l,l,rm]'
                        if direction == 'horizontal':
                            props['model']['diri']['dofs']   = '[dx,dy,dx]'
                            props['model']['diri']['values'] = '[0.,0.,1.]'
                        elif direction == 'vertical':
                            props['model']['diri']['dofs']   = '[dx,dy,dy]'
                            props['model']['diri']['values'] = '[0.,0.,-1.]'
                        else:
                            raise ValueError('direction must be "horizontal" or "vertical"')
                    elif model == 'poisson':
                        props['model']['diri']['groups'] = '[l,rm]'
                        props['model']['diri']['dofs']   = '[u,u]'
                        props['model']['diri']['values'] = '[0.,1.]'
                    else:
                        raise ValueError('model must be "solid" or "poisson"')
                elif load == 'point load':
                    if model == 'solid':
                        props['model']['diri']['groups'] = '[l,l]'
                        props['model']['diri']['dofs']   = '[dx,dy]'
                        props['model']['diri']['values'] = '[0.,0.]'
                        props['model']['neum']['groups'] = '[rm]'
                        if direction == 'horizontal':
                            props['model']['neum']['dofs']   = '[dx]'
                            props['model']['neum']['values'] = '[1.]'
                        elif direction == 'vertical':
                            props['model']['neum']['dofs']   = '[dy]'
                            props['model']['neum']['values'] = '[-1.]'
                        else:
                            raise ValueError('direction must be "horizontal" or "vertical"')
                    elif model == 'poisson':
                        props['model']['diri']['groups'] = '[l]'
                        props['model']['diri']['dofs']   = '[u]'
                        props['model']['diri']['values'] = '[0.]'
                        props['model']['neum']['groups'] = '[rm]'
                        props['model']['neum']['dofs']   = '[u]'
                        props['model']['neum']['values'] = '[1.]'
                    else:
                        raise ValueError('model must be "solid" or "poisson"')
                else:
                    raise ValueError('load must be "forced displacement" or "point load"')

                props['gpinit']['mesh']['file'] = 'meshes/q9/hole-{}.msh'.format(loc)
                props['gpinit']['coarseMesh']['file'] = 'meshes/q4/hole-{}.msh'.format(loc)
                props['model']['gp']['prior']['func'] = 'alpha**2 * ' + prior

                props_c = {}
                props_c['init'] = deepcopy(props['gpinit'])
                props_c['init']['type'] = 'Init'
                props_c['solver'] = deepcopy(props['gpsolver'])
                props_c['solver']['type'] = 'Linsolve'
                props_c['model'] = deepcopy(props['model'])
                props_c['model']['models'] = '[ ' + model + ' , diri, neum ]'
                props_c['init']['mesh']['file'] = props['gpinit']['coarseMesh']['file']
                props_c['model'][model]['shape']['type'] = 'Quad4'

                props_e = {}
                props_e['gpinit'] = deepcopy(props['gpinit'])
                props_e['gpinit']['mesh']['file'] = 'meshes/q9/hole-{}-r2.msh'.format(loc)
                props_e['gpinit']['coarseMesh']['file'] = 'meshes/q9/hole-{}.msh'.format(loc)
                props_e['solver'] = deepcopy(props['gpsolver'])
                props_e['solver']['type'] = 'Linsolve'
                props_e['model'] = deepcopy(props['model'])
                props_e['model']['models'] = '[ ' + model + ' , gp, diri, neum ]'
                props_e['model']['gp']['shape']['type'] = 'Quad9'

                globdat_c = main.jive(props_c)
                u_c = globdat_c['state0']

                globdat_e = main.jive(props_e)
                u_e = globdat_e['state0']
                Phi_e = globdat_e['Phi']
                u_e = np.linalg.solve(Phi_e.T @ Phi_e, Phi_e.T @ u_e)

                globdat = main.jive(props)
                u = globdat['state0']

                f_prior = globdat['f_prior']
                u_prior = globdat['u_prior']
                f_post = globdat['f_post']
                u_post = globdat['u_post']

                std_f_prior = globdat['std_f_prior']
                std_u_prior = globdat['std_u_prior']
                std_f_post = globdat['std_f_post']
                std_u_post = globdat['std_u_post']

                Phi = globdat['Phi']
                Phi_sub = np.zeros_like(Phi)

                for i in range(Phi.shape[0]):
                    for j in range(Phi.shape[1]):
                        if np.isclose(Phi[i,j], 1.):
                            Phi_sub[i,j] = 1.

                u_c = Phi @ u_c

                err = abs(u - u_c)
                err_post = abs(u - u_post)
                err_exact = abs(u_e - u)
                err_exact_post = abs(u_e - u_post)

                std_u_post_proj = np.linalg.solve(Phi.T @ Phi, Phi.T @ std_u_post)
                std_u_post_sub = Phi_sub.T @ std_u_post

                info = '{}-prior, {}, {}, hole {}'.format(prior, load, direction, loc)
                fname = info.replace(', ', '_').replace(' ', '-')

                if model == 'solid':
                    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(8,12), tight_layout=True)
                    QuickViewer(err, globdat, comp=0, ax=ax1, title=r'Discretization error ($x$) ($|u_f - u_c|$)')
                    QuickViewer(err_exact, globdat, comp=0, ax=ax2, title=r'Discretization error ($x$) ($|u_e - u_f|$)')
                    QuickViewer(err_post, globdat, comp=0, ax=ax3, title=r'Discretization error ($x$) ($|u_f - \bar{u}|$)')
                    QuickViewer(std_u_post, globdat, comp=0, ax=ax4, title=r'Posterior std ($x$) ($\sqrt{\bar \Sigma_{ii}}$)')
                    QuickViewer(std_u_post_sub, globdat_c, comp=0, ax=ax5, title=r'Posterior std (subselected) ($x$) ($\sqrt{\bar \Sigma_{ii}}$)')
                    fig.suptitle(info)
                    plt.savefig(fname='img/'+fname+'_x.png', dpi=450)
                    plt.show()

                    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(8,12), tight_layout=True)
                    QuickViewer(err, globdat, comp=1, ax=ax1, title=r'Discretization error ($y$) ($|u_f - u_c|$)')
                    QuickViewer(err_exact, globdat, comp=1, ax=ax2, title=r'Discretization error ($y$) ($|u_e - u_f|$)')
                    QuickViewer(err_post, globdat, comp=1, ax=ax3, title=r'Discretization error ($y$) ($|u_f - \bar{u}|$)')
                    QuickViewer(std_u_post, globdat, comp=1, ax=ax4, title=r'Posterior std ($y$) ($\sqrt{\bar \Sigma_{ii}}$)')
                    QuickViewer(std_u_post_sub, globdat_c, comp=1, ax=ax5, title=r'Posterior std (subselected) ($y$) ($\sqrt{\bar \Sigma_{ii}}$)')
                    fig.suptitle(info)
                    plt.savefig(fname='img/'+fname+'_y.png', dpi=450)
                    plt.show()

                elif model == 'poisson':
                    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(8,12), tight_layout=True)
                    QuickViewer(err, globdat, comp=0, ax=ax1, title=r'Discretization error ($u$) ($|u_f - u_c|$)')
                    QuickViewer(err_exact, globdat, comp=0, ax=ax2, title=r'Discretization error ($u$) ($|u_e - u_f|$)')
                    QuickViewer(err_post, globdat, comp=0, ax=ax3, title=r'Discretization error ($u$) ($|u_f - \bar{u}|$)')
                    QuickViewer(std_u_post, globdat, comp=0, ax=ax4, title=r'Posterior std ($u$) ($\sqrt{\bar \Sigma_{ii}}$)')
                    QuickViewer(std_u_post_sub, globdat_c, comp=0, ax=ax5, title=r'Posterior std (subselected) ($u$) ($\sqrt{\bar \Sigma_{ii}}$)')
                    fig.suptitle(info)
                    plt.savefig(fname='img/'+fname+'_u.png', dpi=450)
                    plt.show()

                else:
                    raise ValueError('model must be "solid" or "poisson"')

                create_dat(data=[u,u_c,err,err_post,err_exact,err_exact_post,
                                 u_prior,u_post,std_u_prior,std_u_post],
                           headers=['u','u_c','error','error_post','error_exact','error_exact_post',
                                    'u_prior','u_post','std_u_prior','std_u_post'],
                           fname='results/'+fname+'.dat')
