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
    monkeypatch.chdir("examples/bar")


@pytest.fixture
def props():
    return pu.parse_file("bar.pro")


def mesher_lin(L, n):
    dx = L / n
    with open("bar.mesh", "w") as fmesh:
        fmesh.write("nodes (ID, x, [y], [z])\n")
        for i in range(n + 1):
            fmesh.write("%d %f\n" % (i, i * dx))
        fmesh.write("elements (node#1, node#2, [node#3, ...])\n")
        for i in range(n):
            fmesh.write("%d %d\n" % (i, i + 1))


def mesher_quad(L, n):
    dx = L / n / 2
    with open("bar.mesh", "w") as fmesh:
        fmesh.write("nodes (ID, x, [y], [z])\n")
        for i in range(2 * n + 1):
            fmesh.write("%d %f\n" % (i, i * dx))
        fmesh.write("elements (node#1, node#2, [node#3, ...])\n")
        for i in range(n):
            fmesh.write("%d %d %d\n" % (2 * i, 2 * i + 1, 2 * i + 2))


@pytest.mark.rank1
@pytest.mark.bar
@pytest.mark.core
def test_lin(props):
    props["model"]["bar"]["shape"]["type"] = "Line2"
    props["model"]["bar"]["shape"]["intScheme"] = "Gauss2"

    mesher_lin(10, 64)

    globdat = main.jive(props)

    elems = globdat["elemSet"]
    nodes = globdat["nodeSet"]

    assert len(nodes) == 65
    assert len(elems) == 64

    K = globdat["matrix0"]
    u = globdat["state0"]
    f = globdat["extForce"]
    c = globdat["constraints"]

    conman = Constrainer(c, K)
    Kc = conman.get_output_matrix()
    fc = conman.get_rhs(f)

    # Check solver solution
    assert np.isclose(Kc @ u, fc).all()

    u_left = u[0]

    # Check displacement field
    assert np.isclose(u_left, -0.9989842929382093)
    assert np.isclose(min(u), u_left)
    assert np.isclose(max(u), 0)


@pytest.mark.rank1
@pytest.mark.bar
@pytest.mark.core
def test_quad(props):
    props["model"]["bar"]["shape"]["type"] = "Line3"
    props["model"]["bar"]["shape"]["intScheme"] = "Gauss3"
    mesher_quad(10, 64)

    globdat = main.jive(props)

    elems = globdat["elemSet"]
    nodes = globdat["nodeSet"]

    assert len(nodes) == 129
    assert len(elems) == 64

    K = globdat["matrix0"]
    u = globdat["state0"]
    f = globdat["extForce"]
    c = globdat["constraints"]

    conman = Constrainer(c, K)
    Kc = conman.get_output_matrix()
    fc = conman.get_rhs(f)

    # Check solver solution
    assert np.isclose(Kc @ u, fc).all()

    u_left = u[0]

    assert np.isclose(u_left, -0.9999995829648631)
    assert np.isclose(min(u), u_left)
    assert np.isclose(max(u), 0)
