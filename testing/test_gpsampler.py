import pytest
import os

import numpy as np

import jive.util.proputils as pu
from jive.app import main

def compare_gpexact_gpsampler(props):
    # First run the exact solve
    props['gpsolver']['type'] = 'GPSolver'
    props['gpsolver']['nsample'] = '0'
    props['gpsolver']['seed'] = '0'
    props['gpsolver']['solver'] = {'type':'cholmod'}
    props['model']['gp']['explicitInverse'] = 'True'

    globdat = main.jive(props)

    c = globdat['constraints']
    cdofs = c.get_constraints()[0]

    # Get the exact prior and posterior means
    mean = globdat['gp']['mean']
    u_prior = mean['prior']['state0']
    f_prior = mean['prior']['extForce']
    u_post = mean['posterior']['state0']
    f_post = mean['posterior']['extForce']

    # Get the exact prior and posterior standard deviations
    std = globdat['gp']['std']
    std_u_prior = std['prior']['state0']
    std_f_prior = std['prior']['extForce']
    std_u_post = std['posterior']['state0']
    std_f_post = std['posterior']['extForce']

    # Now, run the sampled solve
    props['gpsolver']['type'] = 'GPSampler'
    props['gpsolver']['nsample'] = '1000'
    props['gpsolver']['seed'] = '0'
    props['gpsolver']['solver'] = {'type':'cholmod'}
    props['model']['gp']['explicitInverse'] = 'False'
    props['model']['gp']['solver'] = {'type':'cholmod'}

    globdat = main.jive(props)

    # Get the sampled prior and posterior mean
    mean_s = globdat['gp']['mean']
    u_prior_s = mean_s['prior']['state0']
    f_prior_s = mean_s['prior']['extForce']
    u_post_s = mean_s['posterior']['state0']
    f_post_s = mean_s['posterior']['extForce']

    # Get the sampled prior and posterior standard deviations
    std_s = globdat['gp']['std']
    std_u_prior_s = std_s['prior']['state0']
    std_f_prior_s = std_s['prior']['extForce']
    std_u_post_s = std_s['posterior']['state0']
    std_f_post_s = std_s['posterior']['extForce']

    # Check if the sample mean deviates less than 0.15 std from the true mean
    tol = 0.15
    assert np.all(abs(np.delete((f_prior_s-f_prior) / std_f_prior, cdofs)) < tol)
    assert np.all(abs(np.delete((u_prior_s-u_prior) / std_u_prior, cdofs)) < tol)
    assert np.all(abs(np.delete((f_post_s-f_post) / std_f_post, cdofs)) < tol)
    assert np.all(abs(np.delete((u_post_s-u_post) / std_u_post, cdofs)) < tol)

    pdnoise = 1e-8

    # Check if the sample std deviates less than 0.1 std from the true std
    tol = 0.1
    assert np.all(abs((std_f_prior_s-std_f_prior) / std_f_prior) < tol)
    assert np.isclose(std_u_prior_s[cdofs], pdnoise).all()
    assert np.all(abs(np.delete((std_u_prior_s-std_u_prior) / std_u_prior, cdofs)) < tol)
    assert np.all(abs((std_f_post_s-std_f_post) / std_f_post) < tol)
    assert np.isclose(std_u_post_s[cdofs], pdnoise).all()
    assert np.all(abs(np.delete((std_u_post_s-std_u_post) / std_u_post, cdofs)) < tol)

def change_test_dir(target, monkeypatch):
    cwd = os.getcwd()
    pyjivedir = cwd[:cwd.rfind('/bfem')] + '/bfem'
    monkeypatch.chdir(pyjivedir)
    monkeypatch.chdir('gpexamples/' + target)

@pytest.fixture
def props_bar(monkeypatch):
    change_test_dir('gpbar', monkeypatch)
    return pu.parse_file('2nodebar.pro')

@pytest.fixture
def props_beam(monkeypatch):
    change_test_dir('gpbeam', monkeypatch)
    return pu.parse_file('beam.pro')

@pytest.fixture
def props_cantilever(monkeypatch):
    change_test_dir('gpcantilever', monkeypatch)
    return pu.parse_file('cantilever.pro')

@pytest.fixture
def props_poisson(monkeypatch):
    change_test_dir('gppoisson', monkeypatch)
    return pu.parse_file('poisson.pro')

@pytest.fixture
def props_tapered(monkeypatch):
    change_test_dir('gptapered', monkeypatch)
    return pu.parse_file('tapered.pro')

@pytest.mark.rank1
@pytest.mark.bar
@pytest.mark.gp
@pytest.mark.sampling
@pytest.mark.slow
def test_sampler_bar(props_bar):
    compare_gpexact_gpsampler(props_bar)

@pytest.mark.rank2
@pytest.mark.beam
@pytest.mark.gp
@pytest.mark.sampling
@pytest.mark.slow
def test_sampler_beam(props_beam):
    compare_gpexact_gpsampler(props_beam)

@pytest.mark.rank2
@pytest.mark.cantilever
@pytest.mark.gp
@pytest.mark.sampling
@pytest.mark.slow
def test_sampler_cantilever(props_cantilever):
    compare_gpexact_gpsampler(props_cantilever)

@pytest.mark.rank2
@pytest.mark.poisson
@pytest.mark.gp
@pytest.mark.sampling
@pytest.mark.slow
def test_sampler_poisson(props_poisson):
    compare_gpexact_gpsampler(props_poisson)

@pytest.mark.rank1
@pytest.mark.tapered
@pytest.mark.gp
@pytest.mark.sampling
@pytest.mark.slow
def test_sampler_tapered(props_tapered):
    compare_gpexact_gpsampler(props_tapered)
