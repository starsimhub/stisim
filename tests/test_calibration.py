"""
Test calibration of an STIsim sim with Starsim's sim.fit(): parameters given as ss.par
with bounds, by name through sti.Sim, and the data attached to the sim
"""
import sciris as sc
import starsim as ss
import hivsim

n_agents = 500
do_plot = False
sc.options(interactive=False)


def make_sim(**kwargs):
    """Create a minimal Zimbabwe HIV sim for calibration tests."""
    sim = hivsim.demo('zimbabwe', run=False, n_agents=n_agents, dt=ss.year, stop=2000, verbose=0, **kwargs)
    sim.data = sim.data[['year', 'hiv.prevalence', 'hiv.new_infections']]
    return sim


@sc.timer()
def test_bounded_pars(do_plot=do_plot):
    """Parameters with bounds given to sti.Sim by name reach the disease and the network with their bounds."""
    sc.heading('Testing parameters with bounds')
    sim = make_sim(beta_m2f=ss.par(0.035, bounds=[0.01, 0.10]), prop_f0=ss.par(0.79, bounds=[0.6, 0.9]))
    sim.init()
    fitted, _ = ss.fit.bounded_pars(sim)
    assert set(fitted) == {'beta_m2f', 'prop_f0'}, f'Expected beta_m2f and prop_f0 to be found with their bounds, got {list(fitted)}'
    assert sim.diseases.hiv.value('beta_m2f') == 0.035, 'Expected the module to see the value'
    return sim


@sc.timer()
def test_calibration(do_plot=do_plot):
    """Test end-to-end calibration, and running the best parameter sets."""
    sc.heading('Testing calibration')
    sim = make_sim(beta_m2f=ss.par(0.035, bounds=[0.01, 0.10]), eff_condom=ss.par(0.95, bounds=[0.5, 0.99]))
    fit = sim.fit(n_trials=3, verbose=False)
    assert set(fit.pars) == {'beta_m2f', 'eff_condom'}, f'Expected both parameters to be fitted, got {list(fit.pars)}'
    trials = fit.to_df()
    assert len(trials) == 3, f'Expected 3 trials, got {len(trials)}'

    # Run the best two parameter sets, as STIsim's make_calib_sims() did
    best = trials.nsmallest(2, 'objective').drop(columns='objective').to_dict('records')
    msim = ss.parallel(fit.sim, pars=best, verbose=0)
    assert len(msim.sims) == 2, f'Expected 2 runs, got {len(msim.sims)}'
    if do_plot:
        fit.plot()
    return fit


if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    T = sc.timer()

    test_bounded_pars(do_plot=do_plot)
    fit = test_calibration(do_plot=do_plot)

    T.toc()
    print('Done.')
