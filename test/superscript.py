import sys

sys.path.append('../')

import numpy as np
import main
import proputils as pu
import testutils as tu
from names import GlobNames as gn

props = pu.parse_file('beam.pro')

H = 2
L = 10
t = float(props['model']['elastic']['thickness'])
E = float(props['model']['elastic']['young'])
EI = E*H**3*t/12
F = 1
uexact = F*L**3/48/EI

globdat = main.jive(props)
K = globdat['matrix0']
u = globdat['state0']
f = np.matmul(K, u)
bodyforces_y = f[715:]
reactions_y = f[713:715]

globdat = main.jive(props)
view = globdat[gn.MODULEFACTORY].get_module('View','view')

props['view'] = {}
props['view']['plot'] = 'stresses[stress_xy]'
props['view']['deform'] = 0
props['view']['ncolors'] = 100

view.init(props, globdat)
status = view.run(globdat)
