from myjive.app import main
from myjive.solver import Constrainer
import myjive.util.proputils as pu
from myjivex import declare_all as declarex

props = pu.parse_file("elem.pro")

H = 1
L = 1
t = float(props["model"]["solid"]["thickness"])
E = float(props["model"]["solid"]["material"]["E"])

globdat = main.jive(props, extra_declares=[declarex])
u = globdat["state0"]
f = globdat["extForce"]
K = globdat["matrix0"]
c = globdat["constraints"]

conman = Constrainer(c, K)
Kc = conman.get_output_matrix()
