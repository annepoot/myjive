import pytest
import os

import numpy as np

import jive.util.proputils as pu
from jive.app import main
from jive.solver.constrainer import Constrainer

@pytest.fixture(autouse=True)
def change_test_dir(monkeypatch):
    cwd = os.getcwd()
    pyjivedir = cwd[:cwd.rfind('/bfem')] + '/bfem'
    monkeypatch.chdir(pyjivedir)
    monkeypatch.chdir('gpexamples/gpbeam')

@pytest.fixture
def props():
    return pu.parse_file('beam.pro')

@pytest.mark.rank2
@pytest.mark.beam
@pytest.mark.gp
@pytest.mark.slow
def test_moments(props):
    globdat = main.jive(props)

    K = globdat['matrix0']
    u = globdat['state0']
    f = globdat['extForce']
    c = globdat['constraints']
    cdofs = c.get_constraints()[0]

    conman = Constrainer(c, K)
    Kc = conman.get_output_matrix()
    fc = conman.get_rhs(f)

    # Check solver solution
    assert np.isclose(Kc @ u, fc).all()

    mid = len(u)//2

    mean = globdat['gp']['mean']
    u_prior = mean['prior']['state0']
    f_prior = mean['prior']['extForce']
    u_post = mean['posterior']['state0']
    f_post = mean['posterior']['extForce']

    assert np.isclose(f_prior, 0).all()
    assert np.isclose(u_prior, 0).all()
    assert np.isclose(f_post[:mid], 0).all()
    assert np.isclose(min(f_post[mid:]), -0.009175998993927223)
    assert np.isclose(max(f_post[mid:]), 0)
    assert np.isclose(f_post[:mid], 0).all()
    assert np.isclose(min(f_post[mid:]), -0.009175998993927223)
    assert np.isclose(max(f_post[mid:]), 0)
    assert np.isclose(min(u_post[:mid]), 0)
    assert np.isclose(max(u_post[:mid]), 0.027936579812533398)
    assert np.isclose(min(u_post[mid:]), -0.045370434536314404)
    assert np.isclose(max(u_post[mid:]), 0)

    std = globdat['gp']['std']
    std_u_prior = std['prior']['state0']
    std_f_prior = std['prior']['extForce']
    std_u_post = std['posterior']['state0']
    std_f_post = std['posterior']['extForce']

    pdnoise = 1e-8

    assert np.isclose(std_f_prior[cdofs], pdnoise).all()
    assert np.isclose(min(np.delete(std_f_prior, cdofs)), 0.010944043765354355)
    assert np.isclose(max(np.delete(std_f_prior, cdofs)), 0.01895563984188281)
    assert np.isclose(std_u_prior[cdofs], pdnoise).all()
    assert np.isclose(min(np.delete(std_u_prior, cdofs)), 0.0003913478691337096)
    assert np.isclose(max(np.delete(std_u_prior, cdofs)), 0.010131784571359227)

    assert np.isclose(std_f_post[cdofs], pdnoise).all()
    assert np.isclose(min(np.delete(std_f_post, cdofs)), 0.010422140625598358)
    assert np.isclose(max(np.delete(std_f_post, cdofs)), 0.01894185273985648)
    assert np.isclose(std_u_post[cdofs], pdnoise).all()
    assert np.isclose(min(np.delete(std_u_post, cdofs)), 0.0000832662709979621)
    assert np.isclose(max(np.delete(std_u_post, cdofs)), 0.00024539965935899366)
