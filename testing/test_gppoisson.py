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

    mean = globdat['gp']['mean']
    u_prior = mean['prior']['state0']
    f_prior = mean['prior']['extForce']
    u_post = mean['posterior']['state0']
    f_post = mean['posterior']['extForce']

    assert np.isclose(f_prior, 0).all()
    assert np.isclose(u_prior, 0).all()
    assert np.isclose(min(f_post), -0.009424503380738927)
    assert np.isclose(max(f_post), 0.014600669208726491)
    assert np.isclose(min(u_post), -0.22654577520208105)
    assert np.isclose(max(u_post), 0.3634450816888768)

    std = globdat['gp']['std']
    std_u_prior = std['prior']['state0']
    std_f_prior = std['prior']['extForce']
    std_u_post = std['posterior']['state0']
    std_f_post = std['posterior']['extForce']

    pdnoise = 1e-8

    assert np.isclose(std_f_prior[cdofs], pdnoise).all()
    assert np.isclose(min(np.delete(std_f_prior, cdofs)), 0.007070936127273382)
    assert np.isclose(max(np.delete(std_f_prior, cdofs)), 0.023273413258675314)
    assert np.isclose(std_u_prior[cdofs], pdnoise).all()
    assert np.isclose(min(np.delete(std_u_prior, cdofs)), 0.020402828323300988)
    assert np.isclose(max(np.delete(std_u_prior, cdofs)), 0.23877100696987844)

    assert np.isclose(std_f_post[cdofs], pdnoise).all()
    assert np.isclose(min(np.delete(std_f_post, cdofs)), 0.004781054493305154)
    assert np.isclose(max(np.delete(std_f_post, cdofs)), 0.015630087563907064)
    assert np.isclose(std_u_post[cdofs], pdnoise).all()
    assert np.isclose(min(np.delete(std_u_post, cdofs)), 0.0012844150074634022)
    assert np.isclose(max(np.delete(std_u_post, cdofs)), 0.0074937822133110295)
