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
    monkeypatch.chdir("examples/solidbeam")


@pytest.fixture
def props1():
    return pu.parse_file("beam.pro")


@pytest.fixture
def props2():
    return pu.parse_file("2partbeam.pro")


@pytest.mark.rank2
@pytest.mark.solidbeam
@pytest.mark.core
def test_solidbeam(props1):
    props1["model"]["diri"]["groups"] = "[lb,lb,rb]"
    props1["model"]["diri"]["dofs"] = "[dx,dy,dy]"
    props1["model"]["diri"]["values"] = "[0,0,0]"

    globdat = main.jive(props1)

    K = globdat["matrix0"]
    u = globdat["state0"]
    f = globdat["extForce"]
    c = globdat["constraints"]

    conman = Constrainer(c, K)
    Kc = conman.get_output_matrix()
    fc = conman.get_rhs(f)

    # Check solver solution
    assert np.isclose(Kc @ u, fc).all()

    u_mid = u[globdat["dofSpace"].get_dof(3, "dy")]
    u_y = u[len(u) // 2 :]

    # Check displacement field
    assert np.isclose(u_mid, -0.018766977331516582)
    assert np.isclose(min(u_y), u_mid)
    assert np.isclose(max(u_y), 0)

    f = K @ u
    bodyforces_y = f[715:]
    reactions_y = f[713:715]

    # Check force equilibrium
    assert np.isclose(sum(reactions_y), 1)
    assert np.isclose(-sum(bodyforces_y), sum(reactions_y))


@pytest.mark.rank2
@pytest.mark.solidbeam
@pytest.mark.core
def test_2partbeam(props2):
    props2["model"]["diri"]["groups"] = "[lb,lb,rb]"
    props2["model"]["diri"]["dofs"] = "[dx,dy,dy]"
    props2["model"]["diri"]["values"] = "[0,0,0]"

    globdat = main.jive(props2)

    K = globdat["matrix0"]
    u = globdat["state0"]
    f = globdat["extForce"]
    c = globdat["constraints"]

    conman = Constrainer(c, K)
    Kc = conman.get_output_matrix()
    fc = conman.get_rhs(f)

    # Check solver solution
    assert np.isclose(Kc @ u, fc).all()

    dofs_y = globdat["dofSpace"].get_dofs(np.arange(len(globdat["nodeSet"])), ["dy"])
    u_y = u[dofs_y]

    dof_mid = globdat["dofSpace"].get_dof(4, "dy")
    u_mid = u[dof_mid]

    # Check displacement field
    # Note that this beam is asymmetrical, so the maximum displacement is not in the middle
    assert np.isclose(u_mid, -0.009092623429924918)
    assert not np.isclose(min(u_y), u_mid)
    assert np.isclose(max(u_y), 0)

    f = K @ u

    dofs_y.remove(337)
    dofs_y.remove(1002)

    bodyforces_y = f[dofs_y]
    reactions_y = f[[337, 1002]]

    # Check force equilibrium
    assert np.isclose(sum(reactions_y), 1)
    assert np.isclose(-sum(bodyforces_y), sum(reactions_y))
