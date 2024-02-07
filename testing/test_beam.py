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
    monkeypatch.chdir("examples/beam")


@pytest.fixture
def props():
    return pu.parse_file("beam.pro")


@pytest.mark.rank2
@pytest.mark.beam
@pytest.mark.core
def test_point_load(props):
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

    u_mid = u[globdat["dofSpace"].get_dof(3, "dy")]
    u_y = u[len(u) // 2 :]

    # Check displacement field
    assert np.isclose(u_mid, -0.01069531888112345)
    assert np.isclose(min(u_y), u_mid)
    assert np.isclose(max(u_y), 0)

    f = K @ u
    bodyforces_y = f[715:]
    reactions_y = f[713:715]

    # Check force equilibrium
    assert np.isclose(sum(reactions_y), 1)
    assert np.isclose(-sum(bodyforces_y), sum(reactions_y))


@pytest.mark.rank2
@pytest.mark.beam
@pytest.mark.core
def test_point_load_roll(props):
    props["model"]["diri"]["groups"] = "[lb,lb,rb]"
    props["model"]["diri"]["dofs"] = "[dx,dy,dy]"
    props["model"]["diri"]["values"] = "[0,0,0]"

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

    u_mid = u[globdat["dofSpace"].get_dof(3, "dy")]
    u_y = u[len(u) // 2 :]

    # Check displacement field
    assert np.isclose(u_mid, -0.01876697733149987)
    assert np.isclose(min(u_y), u_mid)
    assert np.isclose(max(u_y), 0)

    f = K @ u
    bodyforces_y = f[715:]
    reactions_y = f[713:715]

    # Check force equilibrium
    assert np.isclose(sum(reactions_y), 1)
    assert np.isclose(-sum(bodyforces_y), sum(reactions_y))


@pytest.mark.rank2
@pytest.mark.beam
@pytest.mark.core
def test_body_load(props):
    props["model"]["diri"]["groups"] = "[lb,lb,rb]"
    props["model"]["diri"]["dofs"] = "[dx,dy,dy]"
    props["model"]["diri"]["values"] = "[0,0,0]"

    props["model"]["neum"]["values"] = "[0.0]"
    props["model"]["load"]["values"] = "[-0.2]"

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

    u_mid = u[globdat["dofSpace"].get_dof(3, "dy")]
    u_y = u[len(u) // 2 :]

    # Check displacement field
    assert np.isclose(u_mid, -0.046080900691179726)
    assert np.isclose(max(u_y), 0)

    f = K @ u
    bodyforces_y = f[715:]
    reactions_y = f[713:715]

    # Check force equilibrium
    assert np.isclose(sum(reactions_y), 3.9962615448382706)
    assert np.isclose(-sum(bodyforces_y), sum(reactions_y))
