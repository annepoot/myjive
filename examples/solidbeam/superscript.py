import sys, os

cwd = os.getcwd()
rootdir = os.path.join(cwd[: cwd.rfind(os.path.sep + "myjive")], "myjive")
if rootdir not in sys.path:
    sys.path.append(rootdir)

import numpy as np
from jive.app import main
from jive.util import proputils as pu
from util.quickviewer import QuickViewer

props = pu.parse_file("beam.pro")

H = 2
L = 10
t = float(props["model"]["solid"]["thickness"])
E = float(props["model"]["solid"]["material"]["E"])
EI = E * H**3 * t / 12
F = 1
uexact = F * L**3 / 48 / EI

print("\n\n first run: as is\n\n")

globdat = main.jive(props)
u = globdat["state0"]
uref = u[globdat["dofSpace"].get_dof(3, "dy")]

print("\n\n second run: with one roller\n\n")

props["model"]["diri"]["groups"] = "[lb,lb,rb]"
props["model"]["diri"]["dofs"] = "[dx,dy,dy]"
props["model"]["diri"]["values"] = "[0,0,0]"

globdat = main.jive(props)
u = globdat["state0"]
urol = u[globdat["dofSpace"].get_dof(3, "dy")]

print(
    "\n\nPoint load check: exact %f, constrained %f, roller %f\n\n"
    % (uexact, uref, urol)
)

rho = 1
q = rho * H * t
weight_exact = H * L * t * rho
props["model"]["neum"]["values"] = "[0.0]"
props["model"]["load"]["values"] = "[-" + str(rho * t) + "]"

globdat = main.jive(props)
K = globdat["matrix0"]
u = globdat["state0"]
f = K @ u
bodyforces_y = f[715:]
reactions_y = f[713:715]

print("Exact weight = ", weight_exact)
print("Weight as sum of fem body forces = ", -np.sum(bodyforces_y))
print("Vertical reactions = ", np.sum(reactions_y))

umid = u[globdat["dofSpace"].get_dof(3, "dy")]
uexact = 5 * q * L**4 / 384 / EI

print("\n\nBody force displacement check: exact %f, computed %f\n\n" % (uexact, umid))

QuickViewer(u, globdat, scale=10, title=r"$u(x)$")
