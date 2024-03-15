from jive.app import main
from jive.solver import Constrainer
import jive.util.proputils as pu
from core import declare_all as declare_core
props = pu.parse_file("elem.pro")

H = 1
L = 1
t = float(props["model"]["solid"]["thickness"])
E = float(props["model"]["solid"]["material"]["E"])

globdat = main.jive(props, extra_declares=[declare_core])
u = globdat["state0"]
f = globdat["extForce"]
K = globdat["matrix0"]
c = globdat["constraints"]

conman = Constrainer(c, K)
Kc = conman.get_output_matrix()
