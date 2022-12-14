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

    std_f_prior = globdat['std_f_prior']
    std_u_prior = globdat['std_u_prior']
    std_f_post = globdat['std_f_post']
    std_u_post = globdat['std_u_post']

    pdnoise = 1e-8

    assert np.isclose(std_f_prior[[0,-1]], pdnoise).all()
    assert np.isclose(std_f_prior[1:-1], 0.3918362293782328).all()
    assert np.isclose(std_u_prior[[0,-1]], pdnoise).all()
    assert np.isclose(std_u_prior[mid], 0.6070367843858733)
    assert np.isclose(max(std_u_prior), std_u_prior[mid])
    assert np.isclose(std_u_prior[[p25,p75]], 0.5949856163367298).all()
    assert np.isclose(min(std_u_prior), pdnoise)

    assert np.isclose(std_f_post[[0,-1]], pdnoise).all()
    assert np.isclose(std_f_post[mid], 0.37701598252510565)
    assert np.isclose(std_f_post[[p25,p75]], 0.3717480581869951).all()
    assert np.isclose(std_u_post[[0,-1]], pdnoise).all()
    assert np.isclose(std_u_post[mid], 0.234519164730104)
    assert np.isclose(std_u_post[[p25,p75]], 0.18731452470532461).all()
    assert np.isclose(min(std_u_post), pdnoise)
