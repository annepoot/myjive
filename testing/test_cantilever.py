import pytest
import sys, os

cwd = os.getcwd()
rootdir = os.path.join(cwd[: cwd.rfind(os.path.sep + "myjive")], "myjive")
if rootdir not in sys.path:
    sys.path.append(rootdir)

import numpy as np
import jive.util.proputils as pu
from jive.app import main
from jive.solver import Constrainer


@pytest.fixture(autouse=True)
def change_test_dir(monkeypatch):
    monkeypatch.chdir(rootdir)
    monkeypatch.chdir("examples/cantilever")


@pytest.fixture
def props():
    return pu.parse_file("beam.pro")


@pytest.mark.rank2
@pytest.mark.cantilever
@pytest.mark.core
def test_cantilever(props):
    globdat = main.jive(props)

    K = globdat["matrix0"]
    u = globdat["state0"]
    f = globdat["extForce"]
    c = globdat["constraints"]

    conman = Constrainer(c, K)
    Kc = conman.get_output_matrix()
    fc = conman.get_rhs(f)

    # Check solver solution
    assert np.isclose(Kc @ u, fc).all()

    u_rt = u[globdat["dofSpace"].get_dof(2, "dy")]
    u_y = u[len(u) // 2 :]

    # Check displacement field
    assert np.isclose(u_rt, -0.5106734913602741)
    assert np.isclose(min(u_y), u_rt)
    assert np.isclose(max(u_y), 0)

    # Check force equilibrium
    assert np.isclose(sum(f), -1)
    assert np.isclose(sum(K @ u), 0)
