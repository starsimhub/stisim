"""
Tests for stisim_examples package
"""
import pytest
import stisim as sti
import hivsim as hs
import stisim_examples as stx
from stisim_examples.loaders import load_location_data
from stisim_examples.demo.sim import make_sim as make_demo_sim
from stisim_examples.zimbabwe.sim import make_sim as make_zim_sim

# Shared defaults for fast tests
kw = dict(n_agents=100, dur=5)


def test_data_loading():
    """Test location listing and data loading."""
    locs = stx.list_locations()
    assert 'demo' in locs
    assert 'zimbabwe' in locs

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
    """Test demo HIV sim via make_sim."""
    sim = make_demo_sim(**kw)

    # Check disease, network, and intervention instances were created
    assert any(isinstance(d, sti.HIV) for d in sim.pars.diseases)
    nw = [n for n in sim.pars.networks if hasattr(n, 'pars') and hasattr(n.pars, 'condom_data')]
    assert len(nw) > 0 and nw[0].pars.condom_data is not None
    intv_names = [intv.name for intv in sim.pars.interventions]
    for name in ['art', 'vmmc', 'fsw_testing', 'other_testing', 'low_cd4_testing']:
        assert name in intv_names, f'{name} not found in {intv_names}'

    sim.run()
    assert len(sim.results.hiv.prevalence) > 0

    # Test hivsim_examples convenience wrapper
    import hivsim_examples as hx
    sim2 = hx.Sim(location='demo', **kw)
    sim2.run()
    assert 'hiv' in sim2.diseases


def test_hiv_zim():
    """Test Zimbabwe HIV sim via make_sim with full data integration."""
    sim = make_zim_sim()

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
    import sciris
    sciris.savefig(filename='my_plot.png')


def test_add_interventions():
    """Test that user-provided interventions are additive to data-driven defaults."""
    custom_prep = sti.Prep()
    sim = make_demo_sim(interventions=[custom_prep], **kw)

    intv_names = [intv.name for intv in sim.pars.interventions]
    # Should have both the data-driven interventions AND the user's
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


def test_dataloader():
    """Test DataLoader class: loading data and making modules."""
    import sciris as sc

    data_path = sc.thispath(stx.__file__) / 'zimbabwe'
    dl = sti.DataLoader(data_path=data_path, location='zimbabwe', diseases='hiv')
    dl.load()

    # Check data was loaded
    assert 'hiv' in dl.data.diseases
    assert 'condom_use' in dl.data
    assert 'art_coverage' in dl.data
    assert 'vmmc_coverage' in dl.data
    assert 'hiv_data' in dl.data

    # Check module creation
    modules = dl.make_modules()
    assert any(isinstance(d, sti.HIV) for d in modules.diseases)
    assert len(modules.networks) > 0
    assert len(modules.interventions) > 0
    assert modules.data is not None

    # Test sti.Sim with data_path directly
    sim = sti.Sim(demographics='zimbabwe', diseases='hiv', data_path=data_path, **kw)
    sim.run()
    assert sim.data is not None
    assert len(sim.results.hiv.prevalence) > 0


def test_sim_equivalence():
    """Test that make_sim and sti.Sim(data_path=...) give equivalent results."""
    import sciris as sc

    seed = 42
    sim1 = make_zim_sim(rand_seed=seed, **kw)
    sim1.run()

    # Same thing via hivsim_examples
    import hivsim_examples as hx
    sim2 = hx.Sim(location='zimbabwe', rand_seed=seed, **kw)
    sim2.run()

    prev1 = sim1.results.hiv.prevalence[:]
    prev2 = sim2.results.hiv.prevalence[:]
    assert all(prev1 == prev2), 'make_zim_sim and hx.Sim differ'


if __name__ == '__main__':
    test_data_loading()
    test_hiv_demo()
    test_hiv_zim()
    test_add_interventions()
    test_demographics()
    test_dataloader()
    test_sim_equivalence()
    print('All tests passed!')
