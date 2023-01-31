import sys
sys.path.append('../../')

import numpy as np
from jive.app import main
from jive.solver.constrainer import Constrainer
import jive.util.proputils as pu

props = pu.parse_file('elem.pro')

H = 1
L = 1
t = float(props['model']['solid']['thickness'])
E = float(props['model']['solid']['material']['E'])

globdat = main.jive(props)
u = globdat['state0']
f = globdat['extForce']
K = globdat['matrix0']
c = globdat['constraints']

conman = Constrainer(c, K)
Kc = conman.get_output_matrix()
