"""
Tests for stisim_examples package
"""
import pytest
import stisim as sti
import hivsim as hs
import stisim_examples as stx
from stisim_examples.loaders import load_location_data

# Shared defaults for fast tests
kw = dict(n_agents=100, dur=5)


def test_data_loading():
    """Test location listing and data loading."""
    locs = stx.list_locations()
    assert 'demo' in locs
    assert 'zimbabwe' in locs
    assert 'kenya' in locs

    data = load_location_data('demo', diseases='hiv')
    assert data.location == 'demo'
    assert 'hiv' in data.diseases
    assert 'condom_use' in data
    assert 'art_coverage' in data
    assert 'vmmc_coverage' in data

    # Zimbabwe should also have hiv_data for plotting
    zim_data = load_location_data('zimbabwe', diseases='hiv')
    assert 'hiv_data' in zim_data
    assert 'hiv.prevalence' in zim_data.hiv_data.columns

    with pytest.raises(ValueError, match="not found"):
        load_location_data('invalid_location', diseases='hiv')


def test_hiv_demo():
    """Test demo HIV sim with full data integration and convenience wrappers."""
    sim = stx.Sim(demographics='demo', diseases='hiv', **kw)

    # Check disease, network, and intervention instances were created
    assert any(isinstance(d, sti.HIV) for d in sim.pars.diseases)
    nw = [n for n in sim.pars.networks if hasattr(n, 'pars') and hasattr(n.pars, 'condom_data')]
    assert len(nw) > 0 and nw[0].pars.condom_data is not None
    intv_names = [intv.name for intv in sim.pars.interventions]
    for name in ['art', 'vmmc', 'fsw_testing', 'other_testing', 'low_cd4_testing']:
        assert name in intv_names, f'{name} not found in {intv_names}'

    sim.run()
    assert len(sim.results.hiv.prevalence) > 0

    # Test convenience wrappers produce equivalent results
    import hivsim_examples as hx
    sim2 = stx.HIVSim(location='demo', **kw)
    sim2.run()
    assert 'hiv' in sim2.diseases

    sim3 = hx.Sim(location='demo', **kw)
    sim3.run()
    assert 'hiv' in sim3.diseases


def test_hiv_zim():
    """Test Zimbabwe HIV sim with full data integration and plotting."""
    sim = stx.Sim(demographics='zimbabwe', diseases='hiv')

    assert any(isinstance(d, sti.HIV) for d in sim.pars.diseases)
    nw = [n for n in sim.pars.networks if hasattr(n, 'pars') and hasattr(n.pars, 'condom_data')]
    assert len(nw) > 0 and nw[0].pars.condom_data is not None
    intv_names = [intv.name for intv in sim.pars.interventions]
    assert 'art' in intv_names
    assert 'vmmc' in intv_names
    assert 'fsw_testing' in intv_names

    sim.run()
    assert len(sim.results.hiv.prevalence) > 0

    # Test data attached for plotting
    assert sim.data is not None

    # Test plot_hiv
    import matplotlib
    matplotlib.use('Agg')
    fig = sim.plot(annualize=True)  


def test_add_interventions():
    """Test that user-provided interventions are additive to packaged defaults."""
    custom_prep = sti.Prep()
    sim = stx.Sim(demographics='demo', diseases='hiv', interventions=[custom_prep], **kw)

    intv_names = [intv.name for intv in sim.pars.interventions]
    # Should have both the packaged interventions AND the user's
    assert 'art' in intv_names
    assert 'vmmc' in intv_names
    assert 'prep' in intv_names
    sim.run()


def test_demographics():
    """Test base stisim and hivsim with demographics parameter."""
    try:
        sim1 = sti.Sim(demographics='demo', diseases='hiv', **kw)
        sim1.run()
    except (FileNotFoundError, ValueError) as e:
        pytest.skip(f"Demographic data files not available: {e}")

    try:
        sim2 = hs.Sim(location='demo', **kw)
        sim2.run()
        assert 'hiv' in sim2.diseases
    except (FileNotFoundError, ValueError) as e:
        pytest.skip(f"Demographic data files not available: {e}")

    # Invalid location
    with pytest.raises(ValueError, match="not found"):
        stx.Sim(demographics='invalid_location', diseases='hiv')


def test_errors():
    """Test error handling for missing required parameters."""
    with pytest.raises(ValueError, match="demographics parameter is required"):
        stx.Sim(diseases='hiv')
    with pytest.raises(ValueError, match="diseases parameter is required"):
        stx.Sim(demographics='demo')


def test_sim_creation():
    """Test equivalent ways of constructing the same sim via stisim_examples."""
    import hivsim_examples as hx
    seed = 0

    # These three all go through the stisim_examples factory
    sim1 = stx.Sim(demographics='zimbabwe', diseases='hiv', rand_seed=seed, **kw)
    sim2 = stx.HIVSim(location='zimbabwe', rand_seed=seed, **kw)
    sim3 = hx.Sim(location='zimbabwe', rand_seed=seed, **kw)
    for sim in [sim1, sim2, sim3]:
        sim.run()

    prev1 = sim1.results.hiv.prevalence[:]
    prev2 = sim2.results.hiv.prevalence[:]
    prev3 = sim3.results.hiv.prevalence[:]
    assert all(prev1 == prev2), 'stx.Sim and stx.HIVSim differ'
    assert all(prev1 == prev3), 'stx.Sim and hx.Sim differ'


if __name__ == '__main__':
    test_data_loading()
    test_hiv_demo()
    test_hiv_zim()
    test_add_interventions()
    test_demographics()
    test_errors()
    test_sim_creation()
    print('All tests passed!')


