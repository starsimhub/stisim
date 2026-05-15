"""
Test calibration API: dot-notation parameter routing, get_pars, make_calib_sims
"""
import sciris as sc
import numpy as np
import starsim as ss
import stisim as sti
import hivsim

n_agents = 500
do_plot = False
sc.options(interactive=False)


def make_sim():
    """Create a minimal Zimbabwe HIV sim for calibration tests."""
    return hivsim.demo('zimbabwe', run=False, n_agents=n_agents, dt=ss.year, dur=10)


@sc.timer()
def test_set_sim_pars(do_plot=do_plot):
    """Test that set_sim_pars routes dot-notation keys to the correct modules."""
    sc.heading('Testing set_sim_pars')
    sim = make_sim()

    sti.set_sim_pars(sim, {
        'hiv.beta_m2f': 0.123,
        'hiv.eff_condom': 0.88,
        'fsw_testing.rel_test': 2.5,
        'rand_seed': 42,          # metadata — should be skipped
        'mismatch': 0.5,          # metadata — should be skipped
    })

    assert sim.get_module('hiv').pars['beta_m2f'] == 0.123, 'Expected beta_m2f to be set'
    assert sim.get_module('hiv').pars['eff_condom'] == 0.88, 'Expected eff_condom to be set'
    assert sim.get_module('fsw_testing').pars['rel_test'] == 2.5, 'Expected intervention par to be set'
    return sim


@sc.timer()
def test_default_build_fn(do_plot=do_plot):
    """Test that default_build_fn sets pars before init, then initializes."""
    sc.heading('Testing default_build_fn')
    sim = make_sim()

    calib_pars = {
        'hiv.beta_m2f':   dict(low=0.01, high=0.10, value=0.077),
        'hiv.eff_condom':  dict(low=0.5, high=0.99, value=0.88),
    }
    sim = sti.default_build_fn(sim, calib_pars)

    assert sim.initialized, 'Expected sim to be initialized after default_build_fn'
    assert sim.diseases.hiv.pars['beta_m2f'] == 0.077, 'Expected beta_m2f to match calibration value'
    assert sim.diseases.hiv.pars['eff_condom'] == 0.88, 'Expected eff_condom to match calibration value'
    return sim


@sc.timer()
def test_calibration(do_plot=do_plot):
    """Test end-to-end calibration with dot notation and get_pars."""
    sc.heading('Testing calibration')
    sim = make_sim()
    data = sim.data[['year', 'hiv.prevalence', 'hiv.new_infections']].rename(columns={'year': 'time'})

    # Test nested format (flattened internally by Calibration)
    calib_pars = dict(
        hiv=dict(
            beta_m2f=dict(low=0.01, high=0.10, guess=0.035),
            eff_condom=dict(low=0.5, high=0.99, guess=0.95),
        ),
    )

    calib = sti.Calibration(
        sim=sim, calib_pars=calib_pars, data=data,
        total_trials=2, n_workers=2, die=True, reseed=False,
    )
    calib.calibrate()

    # get_pars should return flat dicts without metadata
    par_sets = calib.get_pars()
    assert len(par_sets) == 2, f'Expected 2 par sets, got {len(par_sets)}'
    assert 'hiv.beta_m2f' in par_sets[0], 'Expected dot-notation key in get_pars output'
    assert 'mismatch' not in par_sets[0], 'Expected metadata stripped from get_pars output'

    par_sets_1 = calib.get_pars(n=1)
    assert len(par_sets_1) == 1, f'Expected 1 par set with n=1, got {len(par_sets_1)}'

    return calib


if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    T = sc.timer()

    test_set_sim_pars(do_plot=do_plot)
    test_default_build_fn(do_plot=do_plot)
    calib = test_calibration(do_plot=do_plot)

    T.toc()
    print('Done.')
