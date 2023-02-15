import sys
sys.path.append('../../')

import matplotlib.pyplot as plt
import numpy as np
from jive.app import main
import jive.util.proputils as pu
from quickviewer import QuickViewer
from copy import deepcopy

props = pu.parse_file('prior-investigation.pro')

for load in ['forced displacement', 'point load', 'self weight']:
    for prior in ['M', 'K']:
        if load == 'forced displacement':
            props['model']['models'] = '[solid,gp,diri]'
            props['model']['diri']['groups'] = '[lb,lb,rb,tm]'
            props['model']['diri']['dofs']   = '[dx,dy,dy,dy]'
            props['model']['diri']['values'] = '[0.,0.,0.,-1.]'
        elif load == 'point load':
            props['model']['models'] = '[solid,gp,diri,neum]'
            props['model']['diri']['groups'] = '[lb,lb,rb]'
            props['model']['diri']['dofs']   = '[dx,dy,dy]'
            props['model']['diri']['values'] = '[0.,0.,0.]'
            props['model']['neum']['groups'] = '[tm]'
            props['model']['neum']['dofs']   = '[dy]'
            props['model']['neum']['values'] = '[-1.]'
        elif load == 'self weight':
            props['model']['models'] = '[solid,gp,diri,load]'
            props['model']['diri']['groups'] = '[lb,lb,rb]'
            props['model']['diri']['dofs']   = '[dx,dy,dy]'
            props['model']['diri']['values'] = '[0.,0.,0.]'
            props['model']['load']['dofs']   = '[dy]'
            props['model']['load']['values'] = '[-1.]'
        else:
            raise ValueError('load should be either "forced displacement" or "point load" or "self weight".')

        props['model']['gp']['prior']['func'] = 'alpha**2 * ' + prior

        props_c = {}
        props_c['init'] = deepcopy(props['gpinit'])
        props_c['init']['type'] = 'Init'
        props_c['solver'] = deepcopy(props['gpsolver'])
        props_c['solver']['type'] = 'Linsolve'
        props_c['model'] = deepcopy(props['model'])
        props_c['model']['models'] = props['model']['models'].replace(',gp,', ',')
        props_c['init']['mesh']['file'] = 'meshes/beam_coarse.msh'

        globdat_c = main.jive(props_c)
        u_c = globdat_c['state0']

        globdat = main.jive(props)
        u = globdat['state0']

        u_prior = globdat['u_prior']
        u_post = globdat['u_post']
        std_u_prior = globdat['std_u_prior']
        std_u_post = globdat['std_u_post']

        Phi = globdat['Phi']

        info = '{} prior, {}'.format(prior, load)
        fname = info.replace(', ', '_').replace(' ', '-')

        fine_list = ['post', 'coarse', 'medium', 'fine', 'fine2']
        x_dict = {}
        u_dict = {}

        for fineness in fine_list:

            if fineness != 'post':
                pro = deepcopy(props_c)
                pro['init']['mesh']['file'] = 'meshes/beam_' + fineness + '.msh'
                pro['solver']['type'] = 'Linsolve'

                glob = main.jive(pro)

                dofs = glob['dofSpace']
                elems = glob['elemSet']
                nodes = glob['nodeSet']
                u = glob['state0']

            else:
                dofs = globdat['dofSpace']
                elems = globdat['elemSet']
                nodes = globdat['nodeSet']
                u = globdat['u_post']
                std_u_post = np.sqrt(globdat['var_u_post'].diagonal())
                std_u_bottom = []

            x_bottom = []
            u_bottom = []

            for n, node in enumerate(nodes):
                coords = node.get_coords()

                # Check if the node in located on the bottom row
                if np.isclose(coords[1], 0):
                    x_bottom.append(coords[0])
                    u_bottom.append(u[dofs.get_dof(n, 'dy')])

                    if fineness == 'post':
                        std_u_bottom.append(std_u_post[dofs.get_dof(n, 'dy')])

            if fineness == 'post':
                x_bottom,u_bottom,std_u_bottom =[list(v) for v in zip(*sorted(zip(x_bottom,u_bottom,std_u_bottom)))]
            else:
                x_bottom,u_bottom =[list(v) for v in zip(*sorted(zip(x_bottom,u_bottom)))]

            x_dict[fineness] = x_bottom
            u_dict[fineness] = u_bottom

        plt.figure()

        for fineness in fine_list:
            plt.plot(x_dict[fineness], u_dict[fineness], label=fineness)
            if fineness == 'post':
                u_bar = np.array(u_dict[fineness])
                std_u_bottom = np.array(std_u_bottom)
                std_u_bottom[0] = std_u_bottom[-1] = 0
                plt.fill_between(x_dict[fineness], u_bar - 2*std_u_bottom, u_bar + 2*std_u_bottom, alpha=0.3)

        plt.title(info)
        plt.legend(loc='lower left')
        plt.savefig(fname='img/'+fname+'.png', dpi=450)
        plt.show()
