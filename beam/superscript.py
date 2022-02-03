import sys

sys.path.append('../')

import numpy as np
import main
import proputils as pu

props = pu.parse_file('beam.pro')

H = 2
L = 10
t = float(props['model']['elastic']['thickness'])
E = float(props['model']['elastic']['young'])
EI = E*H**3*t/12
F = 1
uexact = F*L**3/48/EI

print('\n\n first run: as is\n\n')

globdat = main.jive(props)
u = globdat['state0']
uref = u[globdat['dofSpace'].get_dof(3,'dy')]

print('\n\n second run: with one roller\n\n')

props['model']['diri']['groups'] = '[lb,lb,rb]'
props['model']['diri']['dofs'] = '[dx,dy,dy]'
props['model']['diri']['values'] = '[0,0,0]'

globdat = main.jive(props)
u = globdat['state0']
urol = u[globdat['dofSpace'].get_dof(3,'dy')]

print('\n\nPoint load check: exact %f, constrained %f, roller %f\n\n' %(uexact,uref,urol))


rho = 1
q = rho * H * t
weight_exact = H * L * t * rho
props['solver']['storeMatrix'] = 'True'
props['model']['neum']['values'] = '[0.0]'
props['model']['elastic']['rho'] = str(rho)


globdat = main.jive(props)
K = globdat['matrix0']
u = globdat['state0']
f = np.matmul(K, u)
bodyforces_y = f[715:]
reactions_y = f[713:715]

print('Exact weight = ', weight_exact)
print('Weight as sum of fem body forces = ', -np.sum(bodyforces_y))
print('Vertical reactions = ', np.sum(reactions_y))

umid = u[globdat['dofSpace'].get_dof(3,'dy')]
uexact = 5*q*L**4/384/EI

print('\n\nBody force displacement check: exact %f, computed %f\n\n' %(uexact,umid))
