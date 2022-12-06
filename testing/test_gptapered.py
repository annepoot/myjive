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
    monkeypatch.chdir('gpexamples/gptapered')

@pytest.fixture
def props():
    return pu.parse_file('tapered.pro')

@pytest.mark.rank1
@pytest.mark.tapered
@pytest.mark.gp
@pytest.mark.sampling
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

    assert np.isclose(f_prior, 0).all()
    assert np.isclose(u_prior[[0,p25,mid,p75,-1]],
                      [0,0.08152985674807625,0.19461968909080446,0.38045804531309960,1]).all()
    assert np.isclose(max(u_prior), 1)
    assert np.isclose(min(u_prior), 0)

    assert np.isclose(f_post[[0,p25,mid,p75,-1]],
                      [0,0.019810267809282217,0.013532366060860505,0.019810267809282276,0]).all()
    assert np.isclose(u_post[[0,p25,mid,p75,-1]],
                      [0,0.653934988731335,1.2239871504443058,1.6050644073146396,1]).all()

    std_f_prior = np.sqrt(globdat['var_f_prior'].diagonal())
    std_f_post = np.sqrt(globdat['var_f_post'].diagonal())
    std_u_prior = np.sqrt(globdat['var_u_prior'].diagonal())
    std_u_post = np.sqrt(globdat['var_u_post'].diagonal())

    pdnoise = 1e-8

    assert np.isclose(std_f_prior[[0,-1]], pdnoise).all()
    assert np.isclose(std_f_prior[1:-1], 0.049514007157186166).all()
    assert np.isclose(std_u_prior[[0,-1]], pdnoise).all()
    assert np.isclose(std_u_prior[[p25,mid,p75]],
                      [0.3044290651558541,0.574543364169389,0.7347500354133456]).all()

    assert np.isclose(std_f_post[[0,-1]], pdnoise).all()
    assert np.isclose(std_f_post[[p25,mid,p75]],
                      [0.04581972397793565,0.045582927956165295,0.04581972397764517]).all()
    assert np.isclose(std_u_post[[0,-1]], pdnoise).all()
    assert np.isclose(std_u_post[[p25,mid,p75]],
                      [0.022150508336711568,0.05134911035792135,0.10832989885724671]).all()

@pytest.mark.rank1
@pytest.mark.tapered
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
    assert np.all(abs((np.mean(samples_f_prior, axis=1)-f_prior) / std_f_prior) < tol)
    assert np.all(abs((np.mean(samples_u_prior, axis=1)-u_prior) / std_u_prior) < tol)
    assert np.all(abs((np.mean(samples_f_post, axis=1)-f_post) / std_f_post) < tol)
    assert np.all(abs((np.mean(samples_u_post, axis=1)-u_post) / std_u_post) < tol)

    # Check if the sample std deviates less than 0.1 std from the true std
    tol = 0.1
    assert np.all(abs((np.std(samples_f_prior, axis=1)-std_f_prior) / std_f_prior)[1:-1] < tol)
    assert np.isclose(np.std(samples_u_prior[[0,-1]], axis=1), 0).all()
    assert np.all(abs((np.std(samples_u_prior, axis=1)-std_u_prior) / std_u_prior)[1:-1] < tol)
    assert np.all(abs((np.std(samples_f_post, axis=1)-std_f_post) / std_f_post)[1:-1] < tol)
    assert np.isclose(np.std(samples_u_post[[0,-1]], axis=1), 0).all()
    assert np.all(abs((np.std(samples_u_post, axis=1)-std_u_post) / std_u_post)[1:-1] < tol)
