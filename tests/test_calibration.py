"""
Test calibration using the Zimbabwe HIV example
"""
import sciris as sc
import stisim as sti
import hivsim as hs


def test_calibration():
    """Test calibration with hs.demo('zimbabwe') and default build_fn."""
    sc.heading('Beginning calibration')

    sim = hs.demo('zimbabwe', run=False, n_agents=500)

    # Use the UNAIDS data that's already loaded in the zimbabwe example
    data = sim.data[['year', 'hiv.prevalence', 'hiv.new_infections']].rename(columns={'year': 'time'})

    calib_pars = dict(
        hiv_beta_m2f=dict(low=0.01, high=0.10, guess=0.035),
        hiv_eff_condom=dict(low=0.5, high=0.99, guess=0.95),
        nw_f1_conc=dict(low=0.005, high=0.3, guess=0.16),
    )

    calib = sti.Calibration(
        sim=sim,
        calib_pars=calib_pars,
        total_trials=2,
        n_workers=1,
        die=True,
        reseed=False,
        debug=True,
        data=data,
    )

    calib.calibrate()
    calib.check_fit()

    return calib


if __name__ == '__main__':
    T = sc.tic()
    calib = test_calibration()
    sc.toc(T)
    print('Done.')
