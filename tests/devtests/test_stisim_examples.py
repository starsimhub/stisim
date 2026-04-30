"""
Tests for stisim_examples package
"""
import pytest
import sciris as sc
import starsim as ss
import stisim as sti
import hivsim


def test_list_locations():
    """Test listing available locations."""
    import stisim_examples as stx
    locs = stx.list_locations()
    assert 'demo' in locs
    assert 'zimbabwe' in locs
    assert 'kenya' in locs
    print(f'✓ Available locations: {list(locs.keys())}')


def test_load_location_data():
    """Test loading location data."""
    import stisim_examples as stx
    from stisim_examples.loaders import load_location_data

    # Load demo HIV data
    data = load_location_data('demo', diseases='hiv')
    assert data.location == 'demo'
    assert 'hiv' in data.diseases
    assert 'condom_use' in data
    assert 'art_coverage' in data
    print(f'✓ Successfully loaded demo location data')


def test_stx_sim_demo_hiv():
    """Test demo HIV simulation via stisim_examples."""
    import stisim_examples as stx

    sim = stx.Sim(demographics='demo', diseases='hiv', n_agents=100, dur=5)
    sim.run()
    assert 'hiv' in sim.diseases
    assert len(sim.results.hiv.prevalence) > 0
    print(f'✓ Demo HIV simulation ran successfully')


def test_demo_integration():
    """Test that demo sim has fully integrated data (networks, interventions)."""
    import stisim_examples as stx

    sim = stx.Sim(demographics='demo', diseases='hiv', n_agents=200, dur=10)

    # Check disease instances were created (before init, modules are in sim.pars)
    diseases = sim.pars.diseases
    assert isinstance(diseases, list), f'Expected list of diseases, got {type(diseases)}'
    assert any(isinstance(d, sti.HIV) for d in diseases), 'HIV instance not found in diseases'

    # Check network has condom data
    networks = sim.pars.networks
    nw = [n for n in networks if hasattr(n, 'pars') and hasattr(n.pars, 'condom_data')]
    assert len(nw) > 0, 'No network with condom_data found'
    assert nw[0].pars.condom_data is not None, 'Condom data is None'

    # Check interventions exist
    intvs = sim.pars.interventions
    intv_names = [intv.name for intv in intvs]
    assert 'art' in intv_names, f'ART not found in interventions: {intv_names}'
    assert 'vmmc' in intv_names, f'VMMC not found in interventions: {intv_names}'
    assert 'fsw_testing' in intv_names, f'FSW testing not found: {intv_names}'
    assert 'other_testing' in intv_names, f'Other testing not found: {intv_names}'
    assert 'low_cd4_testing' in intv_names, f'Low CD4 testing not found: {intv_names}'

    # Run and verify results
    sim.run()
    assert len(sim.results.hiv.prevalence) > 0
    print(f'✓ Demo integration test passed')


def test_zimbabwe_integration():
    """Test that Zimbabwe sim has fully integrated data."""
    import stisim_examples as stx

    sim = stx.Sim(demographics='zimbabwe', diseases='hiv', n_agents=200, dur=10)

    # Check disease instances
    diseases = sim.pars.diseases
    assert any(isinstance(d, sti.HIV) for d in diseases)

    # Check network has condom data
    networks = sim.pars.networks
    nw = [n for n in networks if hasattr(n, 'pars') and hasattr(n.pars, 'condom_data')]
    assert len(nw) > 0 and nw[0].pars.condom_data is not None

    # Check interventions exist
    intv_names = [intv.name for intv in sim.pars.interventions]
    assert 'art' in intv_names
    assert 'vmmc' in intv_names
    assert 'fsw_testing' in intv_names

    # Run and verify
    sim.run()
    assert len(sim.results.hiv.prevalence) > 0
    print(f'✓ Zimbabwe integration test passed')


def test_user_override_interventions():
    """Test that user-provided interventions override configured defaults."""
    import stisim_examples as stx

    # User provides their own interventions - should use those, not the defaults
    custom_art = sti.ART()
    sim = stx.Sim(demographics='demo', diseases='hiv', n_agents=100, dur=5,
                   interventions=[custom_art])

    # Should only have the user's intervention, not the configured defaults
    assert len(sim.pars.interventions) == 1
    sim.run()
    print(f'✓ User override of interventions works')


def test_stx_hivsim():
    """Test HIVSim convenience function."""
    import stisim_examples as stx

    sim = stx.HIVSim(location='demo', n_agents=100, dur=5)
    sim.run()
    assert 'hiv' in sim.diseases
    print(f'✓ HIVSim convenience function works')


def test_hx_sim():
    """Test hivsim_examples wrapper."""
    import hivsim_examples as hx

    sim = hx.Sim(location='demo', n_agents=100, dur=5)
    sim.run()
    assert 'hiv' in sim.diseases
    print(f'✓ hivsim_examples wrapper works')


def test_base_stisim_demographics():
    """Test base stisim with demographics parameter."""
    # Note: This requires demographic data files to exist
    # For now, we skip this test if data isn't available
    try:
        sim = sti.Sim(demographics='demo', diseases='hiv', n_agents=100, dur=5)
        sim.run()
        print(f'✓ Base stisim with demographics works')
    except (FileNotFoundError, ValueError) as e:
        pytest.skip(f"Demographic data files not available: {e}")


def test_base_hivsim_location():
    """Test base hivsim with location parameter."""
    try:
        sim = hivsim.Sim(location='demo', n_agents=100, dur=5)
        sim.run()
        assert 'hiv' in sim.diseases
        print(f'✓ Base hivsim with location parameter works')
    except (FileNotFoundError, ValueError) as e:
        pytest.skip(f"Demographic data files not available: {e}")


def test_stx_sim_missing_demographics():
    """Test error handling for missing demographics parameter."""
    import stisim_examples as stx

    with pytest.raises(ValueError, match="demographics parameter is required"):
        sim = stx.Sim(diseases='hiv')


def test_stx_sim_missing_diseases():
    """Test error handling for missing diseases parameter."""
    import stisim_examples as stx

    with pytest.raises(ValueError, match="diseases parameter is required"):
        sim = stx.Sim(demographics='demo')


def test_stx_sim_invalid_location():
    """Test error handling for invalid location."""
    import stisim_examples as stx

    with pytest.raises(ValueError, match="not found"):
        sim = stx.Sim(demographics='invalid_location', diseases='hiv')


def test_stx_sim_zimbabwe_hiv():
    """Test Zimbabwe HIV simulation via stisim_examples."""
    import stisim_examples as stx

    sim = stx.Sim(demographics='zimbabwe', diseases='hiv', n_agents=100, dur=5)
    sim.run()
    assert 'hiv' in sim.diseases
    assert len(sim.results.hiv.prevalence) > 0
    print(f'✓ Zimbabwe HIV simulation ran successfully')


def test_hivsim_zimbabwe():
    """Test Zimbabwe via hivsim."""
    sim = hivsim.Sim(location='zimbabwe', n_agents=100, dur=5)
    sim.run()
    assert 'hiv' in sim.diseases
    print(f'✓ hivsim Zimbabwe simulation works')


def test_stx_hivsim_zimbabwe():
    """Test Zimbabwe via stisim_examples HIVSim."""
    import stisim_examples as stx

    sim = stx.HIVSim(location='zimbabwe', n_agents=100, dur=5)
    sim.run()
    assert 'hiv' in sim.diseases
    print(f'✓ stisim_examples HIVSim Zimbabwe works')


if __name__ == '__main__':
    # Run tests
    test_list_locations()
    test_load_location_data()
    test_stx_sim_demo_hiv()
    test_demo_integration()
    test_zimbabwe_integration()
    test_user_override_interventions()
    test_stx_hivsim()
    test_hx_sim()
    test_base_stisim_demographics()
    test_base_hivsim_location()
    test_stx_sim_missing_demographics()
    test_stx_sim_missing_diseases()
    test_stx_sim_invalid_location()
    test_stx_sim_zimbabwe_hiv()
    test_hivsim_zimbabwe()
    test_stx_hivsim_zimbabwe()

    print('\n✅ All tests passed!')
