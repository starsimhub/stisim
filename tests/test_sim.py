"""
Sim constructor tests: verify that sims can be built in various ways
and that parameters get routed to the correct modules.

Not for scientific validation — that lives in test_hiv.py, test_stis.py,
test_networks.py, and test_hiv_interventions.py.
"""
import os
import numpy as np
import starsim as ss
import stisim as sti
import hivsim

TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_data')


# %% Minimal sim tests

def test_minimal_hiv():
    """ Simplest possible HIV sim — no manually specified modules """
    sim = sti.Sim(diseases='hiv', n_agents=500, dur=10)
    sim.run()
    assert sim.results.hiv.cum_infections[-1] > 0
    return sim


def test_minimal_sti():
    """ Simplest possible non-HIV STI sim """
    sim = sti.Sim(diseases='ng', n_agents=500, dur=10)
    sim.run()
    assert sim.results.ng.cum_infections[-1] > 0
    return sim


def test_multi_disease():
    """ Multiple diseases in one sim """
    sim = sti.Sim(diseases=['hiv', 'syph', 'ng'], n_agents=1000, dur=10)
    sim.run()
    assert len(sim.diseases) == 3
    for name in ['hiv', 'syph', 'ng']:
        assert name in sim.diseases
    return sim


def test_time():
    """ Check that simulation time vector is initialized correctly """
    sim = sti.Sim(start=2010, diseases='hiv')
    sim.run()
    assert sim.pars.start == sim.t.yearvec[0]
    return sim


# %% Constructor flexibility tests

def test_sim_creation():
    """ Test various constructor patterns for sti.Sim """

    start = 2010
    stop = 2020

    pars = dict(
        start=start,
        stop=stop,
    )

    nw_pars = dict(debut=ss.lognorm_ex(20, 5))
    sti_pars = dict(ng=dict(eff_condom=0.6))

    # Test 1: default networks with custom pars, demographics from location string, and diseases from disease names with custom pars
    sim1 = sti.Sim(
        pars=pars,
        nw_pars=nw_pars,
        demographics='zimbabwe',
        datafolder=TEST_DATA_DIR,
        diseases=['ng', 'ct', 'tv', 'bv', 'hiv'],
        sti_pars=sti_pars,
    )
    sim1.init()

    assert sim1.diseases.ng.pars.eff_condom == 0.6, "Disease parameter not set correctly"
    assert len(sim1.diseases) == 5, "Incorrect number of diseases initialized"

    # Test 2: mix of strings and modules
    demographics = [ss.Pregnancy(), ss.Deaths()]
    networks = sti.StructuredSexual()
    diseases = [sti.Gonorrhea(), 'hiv']

    sim2 = sti.Sim(
        pars=pars,
        networks=networks,
        demographics=demographics,
        diseases=diseases,
    )

    sim2.init()

    assert isinstance(sim2.networks.structuredsexual, sti.StructuredSexual), "Network not initialized correctly"
    assert len(sim2.diseases) == 2, "Incorrect number of diseases initialized"
    assert len(sim2.demographics) == 2, "Incorrect number of demographics initialized"

    # Test 3: flat pars dict
    pars = dict(
        start=2010,
        beta_m2f=0.05,
        prop_f0=0.45,
        location='zimbabwe',
        datafolder=TEST_DATA_DIR,
        diseases=['ng', 'ct', 'tv'],
        ng=dict(eff_condom=0.6),
    )

    sim3 = sti.Sim(**pars)
    sim3.init()

    assert sim3.diseases.ng.pars.beta_m2f == pars['beta_m2f'], "Disease parameter not set correctly"
    assert sim3.diseases.ct.pars.beta_m2f == pars['beta_m2f'], "Disease parameter not set correctly"
    assert sim3.diseases.ng.pars.eff_condom == pars['ng']['eff_condom'], "Disease parameter not set correctly"
    assert sim3.networks.structuredsexual.pars.prop_f0 == pars['prop_f0'], "Network parameter not set correctly"
    assert len(sim3.networks) == 2, "Default networks not added"
    assert len(sim3.diseases) == 3, "Incorrect number of diseases initialized"

    return


# %% hivsim.Sim tests (#349)

def test_hivsim_defaults():
    """ hivsim.Sim() with no args should create a complete HIV sim with sensible defaults """
    sim = hivsim.Sim(n_agents=500, dur=10)
    sim.init()

    # Check default modules are present
    assert 'hiv' in sim.diseases, "HIV disease not added"
    assert len(sim.diseases) == 1, "Should have exactly 1 disease"
    assert len(sim.networks) == 2, "Should have 2 networks (sexual + maternal)"
    assert len(sim.demographics) == 2, "Should have 2 demographics (pregnancy + deaths)"
    assert len(sim.interventions) == 4, "Should have 4 interventions (test, ART, VMMC, PrEP)"

    sim.run()
    assert sim.results.hiv.cum_infections[-1] > 0
    return sim


def test_hivsim_pars():
    """ Parameters passed to hivsim.Sim should route to the HIV module """
    beta = 0.05
    sim = hivsim.Sim(beta_m2f=beta, n_agents=500, dur=5)
    sim.init()
    assert sim.diseases.hiv.pars.beta_m2f == beta, f"beta_m2f not routed to HIV: {sim.diseases.hiv.pars.beta_m2f}"
    return sim


def test_hivsim_custom_modules():
    """ User-supplied modules should replace defaults """
    custom_art = sti.ART()
    sim = hivsim.Sim(interventions=[custom_art], n_agents=500, dur=5)
    sim.init()

    # User passed 1 intervention, so defaults should be replaced
    assert len(sim.interventions) == 1, f"Expected 1 intervention, got {len(sim.interventions)}"

    # But other defaults should still be present
    assert len(sim.networks) == 2, "Default networks should still be present"
    assert len(sim.demographics) == 2, "Default demographics should still be present"
    return sim


# %% Demo tests (merged from test_examples.py)

def test_demo_simple():
    """ hivsim.demo('simple') and make_sim() should produce identical results """
    from hivsim_examples.simple.sim import make_sim

    seed = 42
    kw = dict(n_agents=200, dur=5)

    sim1 = make_sim(rand_seed=seed, **kw)
    sim1.run()

    sim2 = hivsim.demo('simple', run=False, rand_seed=seed, **kw)
    sim2.run()

    prev1 = sim1.results.hiv.prevalence[:]
    prev2 = sim2.results.hiv.prevalence[:]
    assert np.allclose(prev1, prev2), 'make_sim and hivsim.demo should produce identical results'


def test_demo_zimbabwe():
    """ hivsim.demo('zimbabwe') and make_sim() should produce identical results """
    from hivsim_examples.zimbabwe.sim import make_sim

    seed = 42
    sim1 = make_sim(rand_seed=seed, n_agents=500, stop=1995)
    sim1.run()

    sim2 = hivsim.demo('zimbabwe', run=False, rand_seed=seed, n_agents=500, stop=1995)
    sim2.run()

    prev1 = sim1.results.hiv.prevalence[:]
    prev2 = sim2.results.hiv.prevalence[:]
    assert np.allclose(prev1, prev2), 'make_sim and hivsim.demo should produce identical results'

    # Check population is scaled to ~10M (Zimbabwe 1990)
    n_alive_start = sim1.results.n_alive[0]
    assert np.isclose(n_alive_start, 9_980_999, rtol=0.05), \
        f'n_alive at t=0 ({n_alive_start}) should be close to 9,980,999'

    # Run a full version (no plot when running via pytest)
    sim3 = hivsim.demo('zimbabwe', plot=False)
    return


# %% Run as a script
if __name__ == '__main__':

    test_minimal_hiv()
    test_minimal_sti()
    test_multi_disease()
    test_time()
    test_sim_creation()
    test_hivsim_defaults()
    test_hivsim_pars()
    test_hivsim_custom_modules()
    test_demo_simple()
    test_demo_zimbabwe()

    print('All tests passed!')
