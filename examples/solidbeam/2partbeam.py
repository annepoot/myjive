import sys, os
rootdir = os.path.abspath(os.path.join(os.getcwd(), "..", ".."))
jivedir = os.path.abspath(os.path.join(rootdir, "..", "myjive"))
sys.path.append(rootdir)
sys.path.append(jivedir)

from jive.app import main
from jive.util import proputils as pu

props = pu.parse_file('2partbeam.pro')

H = 2
L = 10
t = float(props['model']['left']['thickness'])
E_left = float(props['model']['left']['material']['E'])
EI_left = E_left*H**3*t/12
F = 1
# uexact = F*L**3/48/EI
uexact = 1/96*F*L**3/EI_left

globdat = main.jive(props)
u = globdat['state0']
umid = u[globdat['dofSpace'].get_dof(4,'dy')]

print('\n\nPoint load check: exact %f, numerical %f\n\n' %(uexact,umid))
