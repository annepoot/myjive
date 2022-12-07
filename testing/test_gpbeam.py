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

    f_prior = globdat['f_prior']
    u_prior = globdat['u_prior']
    f_post = globdat['f_post']
    u_post = globdat['u_post']

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

    std_f_prior = np.sqrt(globdat['var_f_prior'].diagonal())
    std_f_post = np.sqrt(globdat['var_f_post'].diagonal())
    std_u_prior = np.sqrt(globdat['var_u_prior'].diagonal())
    std_u_post = np.sqrt(globdat['var_u_post'].diagonal())

    pdnoise = 1e-8

    assert np.isclose(std_f_prior[cdofs], pdnoise).all()
    assert np.isclose(min(np.delete(std_f_prior, cdofs)), 0.016049220130768718)
    assert np.isclose(max(np.delete(std_f_prior, cdofs)), 0.0277980650673971)
    assert np.isclose(std_u_prior[cdofs], pdnoise).all()
    assert np.isclose(min(np.delete(std_u_prior, cdofs)), 0.0005739037812811659)
    assert np.isclose(max(np.delete(std_u_prior, cdofs)), 0.014858058354835457)

    assert np.isclose(std_f_post[cdofs], pdnoise).all()
    assert np.isclose(min(np.delete(std_f_post, cdofs)), 0.015283859669512692)
    assert np.isclose(max(np.delete(std_f_post, cdofs)), 0.027777846184013986)
    assert np.isclose(std_u_post[cdofs], pdnoise).all()
    assert np.isclose(min(np.delete(std_u_post, cdofs)), 0.00012210831218971617)
    assert np.isclose(max(np.delete(std_u_post, cdofs)), 0.0003598736662344794)

@pytest.mark.rank2
@pytest.mark.beam
@pytest.mark.gp
@pytest.mark.sampling
@pytest.mark.slow
def test_samples(props):
    props['gpsolver']['nsample'] = '1000'
    props['gpsolver']['seed'] = '0'

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

    samples_f_prior = globdat['samples_f_prior']
    samples_u_prior = globdat['samples_u_prior']
    samples_f_post = globdat['samples_f_post']
    samples_u_post = globdat['samples_u_post']

    f_prior = globdat['f_prior']
    u_prior = globdat['u_prior']
    f_post = globdat['f_post']
    u_post = globdat['u_post']

    std_f_prior = np.sqrt(globdat['var_f_prior'].diagonal())
    std_f_post = np.sqrt(globdat['var_f_post'].diagonal())
    std_u_prior = np.sqrt(globdat['var_u_prior'].diagonal())
    std_u_post = np.sqrt(globdat['var_u_post'].diagonal())

    # Check if the sample mean deviates less than 0.15 std from the true mean
    tol = 0.15
    assert np.all(abs(np.delete((np.mean(samples_f_prior, axis=1)-f_prior) / std_f_prior, cdofs)) < tol)
    assert np.all(abs(np.delete((np.mean(samples_u_prior, axis=1)-u_prior) / std_u_prior, cdofs)) < tol)
    assert np.all(abs(np.delete((np.mean(samples_f_post, axis=1)-f_post) / std_f_post, cdofs)) < tol)
    assert np.all(abs(np.delete((np.mean(samples_u_post, axis=1)-u_post) / std_u_post, cdofs)) < tol)

    pdnoise = 1e-8

    # Check if the sample std deviates less than 0.1 std from the true std
    tol = 0.1
    assert np.all(abs((np.std(samples_f_prior, axis=1)-std_f_prior) / std_f_prior) < tol)
    assert np.isclose(np.std(samples_u_prior[cdofs], axis=1), pdnoise).all()
    assert np.all(abs(np.delete((np.std(samples_u_prior, axis=1)-std_u_prior) / std_u_prior, cdofs)) < tol)
    assert np.all(abs((np.std(samples_f_post, axis=1)-std_f_post) / std_f_post) < tol)
    assert np.isclose(np.std(samples_u_post[cdofs], axis=1), pdnoise).all()
    assert np.all(abs(np.delete((np.std(samples_u_post, axis=1)-std_u_post) / std_u_post, cdofs)) < tol)
