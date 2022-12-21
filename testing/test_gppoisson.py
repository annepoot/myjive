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
    monkeypatch.chdir('gpexamples/gppoisson')

@pytest.fixture
def props():
    props = pu.parse_file('poisson.pro')

    props['gpsolver'] = {}
    props['gpsolver']['type'] = 'GPSolver'
    props['gpsolver']['solver'] = {'type':'cholmod'}

    props['model']['gp']['explicitInverse'] = 'True'

    return props

@pytest.mark.rank2
@pytest.mark.poisson
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

    f_prior = globdat['f_prior']
    u_prior = globdat['u_prior']
    f_post = globdat['f_post']
    u_post = globdat['u_post']

    assert np.isclose(f_prior, 0).all()
    assert np.isclose(u_prior, 0).all()
    assert np.isclose(min(f_post), -0.011230547987848664)
    assert np.isclose(max(f_post), 0.01197235074186509)
    assert np.isclose(min(u_post), -0.2199448680398467)
    assert np.isclose(max(u_post), 0.36923905503709703)

    std_f_prior = globdat['std_f_prior']
    std_u_prior = globdat['std_u_prior']
    std_f_post = globdat['std_f_post']
    std_u_post = globdat['std_u_post']

    pdnoise = 1e-8

    assert np.isclose(std_f_prior[cdofs], pdnoise).all()
    assert np.isclose(min(np.delete(std_f_prior, cdofs)), 0.0024066963042805633)
    assert np.isclose(max(np.delete(std_f_prior, cdofs)), 0.008306578554038014)
    assert np.isclose(std_u_prior[cdofs], pdnoise).all()
    assert np.isclose(min(np.delete(std_u_prior, cdofs)), 0.008823523324660434)
    assert np.isclose(max(np.delete(std_u_prior, cdofs)), 0.09038467965264092)

    assert np.isclose(std_f_post[cdofs], pdnoise).all()
    assert np.isclose(min(np.delete(std_f_post, cdofs)), 0.0017171551839928435)
    assert np.isclose(max(np.delete(std_f_post, cdofs)), 0.00547682965009234)
    assert np.isclose(std_u_post[cdofs], pdnoise).all()
    assert np.isclose(min(np.delete(std_u_post, cdofs)), 0.0004613891118080122)
    assert np.isclose(max(np.delete(std_u_post, cdofs)), 0.0027349786517457527)
