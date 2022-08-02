import sys

sys.path.append('../')

import numpy as np
import main
import proputils as pu
import pyvista

props = pu.parse_file('beam.pro')

H = 2
L = 10
t = float(props['model']['solid']['material']['thickness'])
E = float(props['model']['solid']['material']['E'])
rho = float(props['model']['solid']['material']['rho'])
q = rho * H * t
EI = E*H**3*t/12
uexact = 5*q*L**4/384/EI
weight_exact = H * L * t * rho

globdat = main.jive(props)
K = globdat['matrix0']
u = globdat['state0']
f = K @ u

bodyforces_y = f[715:]
reactions_y = f[713:715]
umid = u[globdat['dofSpace'].get_dof(3,'dy')]

print('Exact weight = ', weight_exact)
print('Weight as sum of fem body forces = ', -np.sum(bodyforces_y))
print('Vertical reactions = ', np.sum(reactions_y))

print('\n\nBody force displacement check: exact %f, computed %f\n\n' %(uexact,umid))

reader = pyvista.get_reader('stiffness1.vtu')
mesh = reader.read()
mesh.plot(scalars='stiffness', cpos='xy')
