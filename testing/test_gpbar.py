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
    monkeypatch.chdir('gpexamples/gpbar')

@pytest.fixture
def props():
    return pu.parse_file('2nodebar.pro')

@pytest.mark.rank1
@pytest.mark.bar
@pytest.mark.gp
def test_moments(props):
    globdat = main.jive(props)

    K = globdat['matrix0']
    u = globdat['state0']
    f = globdat['extForce']
    c = globdat['constraints']

    conman = Constrainer(c, K)
    Kc = conman.get_output_matrix()
    fc = conman.get_rhs(f)

    # Check solver solution
    assert np.isclose(Kc @ u, fc).all()

    p25 = (len(u)-1)//4
    mid = (len(u)-1)//2
    p75 = (len(u)-1)//4*3

    f_prior = globdat['f_prior']
    u_prior = globdat['u_prior']
    f_post = globdat['f_post']
    u_post = globdat['u_post']

    assert np.isclose(f_prior[mid], 5)
    assert np.isclose(np.delete(f_prior, mid), 0).all()
    assert np.isclose(u_prior[mid], 2.4972362848453950)
    assert np.isclose(max(u_prior), u_prior[mid])
    assert np.isclose(min(u_prior), 0)

    assert np.isclose(f_post[mid], 5.1479018481190035)
    assert np.isclose(max(f_post), f_post[mid])
    assert np.isclose(min(f_post), 0)
    assert np.isclose(u_post[mid], 3.4637802350491290)
    assert np.isclose(max(u_post), u_post[mid])
    assert np.isclose(min(u_post), 0)

    std_f_prior = np.sqrt(globdat['var_f_prior'].diagonal())
    std_f_post = np.sqrt(globdat['var_f_post'].diagonal())
    std_u_prior = np.sqrt(globdat['var_u_prior'].diagonal())
    std_u_post = np.sqrt(globdat['var_u_post'].diagonal())

    pdnoise = 1e-8

    assert np.isclose(std_f_prior[[0,-1]], pdnoise).all()
    assert np.isclose(std_f_prior[1:-1], 0.4692721376883415).all()
    assert np.isclose(std_u_prior[[0,-1]], pdnoise).all()
    assert np.isclose(std_u_prior[mid], 0.7270013033665275)
    assert np.isclose(max(std_u_prior), std_u_prior[mid])
    assert np.isclose(std_u_prior[[p25,p75]], 0.7125685455762724).all()
    assert np.isclose(min(std_u_prior), pdnoise)

    assert np.isclose(std_f_post[[0,-1]], pdnoise).all()
    assert np.isclose(std_f_post[mid], 0.45152306702950545)
    assert np.isclose(std_f_post[[p25,p75]], 0.44521407891141995).all()
    assert np.isclose(std_u_post[[0,-1]], pdnoise).all()
    assert np.isclose(std_u_post[mid], 0.2808655798174858)
    assert np.isclose(std_u_post[[p25,p75]], 0.2243322103286353).all()
    assert np.isclose(min(std_u_post), pdnoise)

@pytest.mark.rank1
@pytest.mark.bar
@pytest.mark.gp
@pytest.mark.sampling
def test_samples(props):
    props['gpsolver']['nsample'] = '1000'
    props['gpsolver']['seed'] = '0'

    globdat = main.jive(props)

    K = globdat['matrix0']
    u = globdat['state0']
    f = globdat['extForce']
    c = globdat['constraints']

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

    # Check if the sample mean deviates less than 0.1 std from the true mean
    tol = 0.1
    assert np.all(abs((np.mean(samples_f_prior, axis=1)-f_prior) / std_f_prior)[1:-1] < tol)
    assert np.all(abs((np.mean(samples_u_prior, axis=1)-u_prior) / std_u_prior)[1:-1] < tol)
    assert np.all(abs((np.mean(samples_f_post, axis=1)-f_post) / std_f_post)[1:-1] < tol)
    assert np.all(abs((np.mean(samples_u_post, axis=1)-u_post) / std_u_post)[1:-1] < tol)

    # Check if the sample std deviates less than 0.1 std from the true std
    tol = 0.1
    assert np.all(abs((np.std(samples_f_prior, axis=1)-std_f_prior) / std_f_prior)[1:-1] < tol)
    assert np.isclose(np.std(samples_u_prior[[0,-1]], axis=1), 0).all()
    assert np.all(abs((np.std(samples_u_prior[1:-1], axis=1)-std_u_prior[1:-1]) / std_u_prior[1:-1]) < tol)
    assert np.all(abs((np.std(samples_f_post, axis=1)-std_f_post) / std_f_post)[1:-1] < tol)
    assert np.isclose(np.std(samples_u_post[[0,-1]], axis=1), 0).all()
    assert np.all(abs((np.std(samples_u_post[1:-1], axis=1)-std_u_post[1:-1]) / std_u_post[1:-1]) < tol)
