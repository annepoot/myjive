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

    mean = globdat['gp']['mean']
    u_prior = mean['prior']['state0']
    f_prior = mean['prior']['extForce']
    u_post = mean['posterior']['state0']
    f_post = mean['posterior']['extForce']

    assert np.isclose(f_prior[:-2], 0).all()
    assert np.isclose(u_prior[[0,p25,mid,p75,-1]],
                      [0,0.06279688471198644,0.1508750175699451,0.2995989479395923,1]).all()
    assert np.isclose(max(u_prior), 1)
    assert np.isclose(min(u_prior), 0)

    assert np.isclose(f_post[[0,p25,mid,p75,-1]],
                      [0,0.05943080354952812,0.0405970982094388,0.059430803549521016,1]).all()
    assert np.isclose(u_post[[0,p25,mid,p75,-1]],
                      [0,0.6160305321061196,1.171961930023314,1.5878743681112304,1]).all()

    std = globdat['gp']['std']
    std_u_prior = std['prior']['state0']
    std_f_prior = std['prior']['extForce']
    std_u_post = std['posterior']['state0']
    std_f_post = std['posterior']['extForce']

    pdnoise = 1e-8

    assert np.isclose(std_f_prior[[0,-1]], pdnoise).all()
    assert np.isclose(std_f_prior[1:-1], 0.12677313820948438).all()
    assert np.isclose(std_u_prior[[0,-1]], pdnoise).all()
    assert np.isclose(std_u_prior[[p25,mid,p75]],
                      [0.25076103626305596,0.4889263637601685,0.6710847606302828]).all()

    assert np.isclose(std_f_post[[0,-1]], pdnoise).all()
    assert np.isclose(std_f_post[[p25,mid,p75]],
                      [0.11731448397869852,0.11670820351838872,0.11731448397872997]).all()
    assert np.isclose(std_u_post[[0,-1]], pdnoise).all()
    assert np.isclose(std_u_post[[p25,mid,p75]],
                      [0.027364163749027375,0.06473001483700104,0.13518029401857337]).all()
