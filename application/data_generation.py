import sys

sys.path.append('../')

import numpy as np
import main
import proputils as pu
from names import GlobNames as gn

props = pu.parse_file('beam.pro')

nsamples = 10
rng = np.random.default_rng(0)

with open('damaged_beams.dat', 'w') as f:
    f.write('sample E_pure nu deteriorations intervention node x y dx dy E_true\n')

    for sample in range(nsamples):

        print('\n' + 50 * '=')
        print('\tSAMPLE {} OUT OF {}'.format(sample+1, nsamples))
        print(50 * '=' + '\n')

        E = np.exp(rng.normal(np.log(10000), 0.1))
        nu = np.exp(rng.normal(np.log(0.2), 0.1))

        n_det = int(10 * sample / nsamples)

        if n_det > 8:
            risk = 'demolition'
        elif n_det > 4:
            risk = 'maintenance'
        else:
            risk = 'unnecessary'

        props['model']['solid']['material']['E'] = E
        props['model']['solid']['material']['nu'] = nu
        props['model']['solid']['material']['deteriorations'] = n_det
        props['model']['solid']['material']['seed'] = sample

        globdat = main.jive(props)
        nodes = globdat[gn.NSET]
        u = globdat[gn.STATE0]
        dx = u[:len(u)//2]
        dy = u[len(u)//2:]
        E_true = globdat[gn.TABLES]['stiffness']['']

        for i, node in enumerate(nodes):
            coords = node.get_coords()

            f.write('{sample:d} {E_pure:f} {nu:f} {n_det:d} {risk:s} {node} {x:f} {y:f} {dx:f} {dy:f} {E_true:f}\n'.format(
                sample=sample,
                E_pure=E,
                nu=nu,
                n_det=n_det,
                risk=risk,
                node=i,
                x=coords[0],
                y=coords[1],
                dx=dx[i],
                dy=dy[i],
                E_true=E_true[i]
            ))
