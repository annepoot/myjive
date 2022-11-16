import sys
sys.path.append('../../')

from sksparse import cholmod as cm

from jive.util import proputils as pu
from jive.app import main

import testutils as tu

props = pu.parse_file('beam.pro')
props['init']['mesh']['file'] = 'meshes/beam_fine.msh'
globdat = main.jive(props)

K = globdat['matrix0']
f = globdat['extForce']
c = globdat['constraints']

K, f = c.constrain(K, f)

for ordering in ['natural', 'amd', 'colamd', 'metis', 'nesdis', 'default', 'best']:
    chol = cm.cholesky(K, ordering_method=ordering)
    P = chol.P()
    L = chol.L()
    print('using {} ordering, L has {} non-zero values'.format(ordering, L.getnnz()))
    tu.showmat(L, mask_zeros=True, title=ordering + ' ordering')
