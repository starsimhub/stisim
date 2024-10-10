"""
Test calibration
"""

#%% Imports and settings
import starsim as ss
import sciris as sc
import stisim as sti
import numpy as np
import pylab as pl
import pandas as pd

do_plot = 1
do_save = 0
n_agents = 2e3


#%% Define the tests
def make_sim():

    syph = sti.Syphilis(
        beta={'structuredsexual': [0.5, 0.25], 'maternal': [0.99, 0]},
        init_prev=0.05,
    )
    hiv = sti.HIV(
        beta={'structuredsexual': [1, 1], 'maternal': [1, 0]},
        beta_m2f=0.05,
        beta_f2m=0.025,
        beta_m2c=0.025,
        init_prev=0.15,
    )
    connector = sti.hiv_syph(hiv, syph, rel_sus_hiv_syph=2, rel_trans_hiv_syph=2)
    pregnancy = ss.Pregnancy(fertility_rate=20)
    death = ss.Deaths(death_rate=10)
    sexual = sti.StructuredSexual(prop_f1=0.2)
    maternal = ss.MaternalNet()

    sim = ss.Sim(
        dt=1,
        n_agents=n_agents,
        total_pop=9980999,
        start=1990,
        dur=40,
        diseases=[syph, hiv],
        networks=[sexual, maternal],
        connectors=connector,
        demographics=[pregnancy, death],
    )

    return sim


def test_calibration(do_plot=True):

    sc.heading('Testing calibration')

    # Define the calibration parameters
    calib_pars = dict(
        beta_m2f = dict(low=0.01, high=0.10, guess=0.05, path=('diseases', 'hiv', 'beta_m2f')),
        beta_f2m = dict(low=0.005, high=0.05, guess=0.025, path=('diseases', 'hiv', 'beta_f2m')),
        rel_trans_latent = dict(low=0.00, high=0.2, guess=0.1, path=('diseases', 'syphilis', 'rel_trans_latent')),
        prop_f1 = dict(low=0.1, high=0.45, guess=0.15, path=('networks', 'structuredsexual', 'prop_f1')),
        prop_m1 = dict(low=0.15, high=0.5, guess=0.21, path=('networks', 'structuredsexual', 'prop_m1')),
        f1_conc = dict(low=0.005, high=0.1, guess=0.01, path=('networks', 'structuredsexual', 'f1_conc')),
        m1_conc = dict(low=0.005, high=0.1, guess=0.01, path=('networks', 'structuredsexual', 'f1_conc')),
        p_pair_form = dict(low=0.4, high=0.9, guess=0.5, path=('networks', 'structuredsexual', 'p_pair_form')),
    )

    # Make the sim
    sim = make_sim()

    # Weight the data
    weights = {
        'n_alive': 1,
        'hiv.prevalence': 1,
        'hiv.n_infected': 1,
        'hiv.new_infections': 1,
        'hiv.new_deaths': 1,
        'syphilis.prevalence': 1
        }

    data = pd.read_csv(sti.root/'tests'/'test_data'/'zimbabwe_calib.csv')

    # Make the calibration
    calib = ss.Calibration(
        calib_pars=calib_pars,
        sim=sim,
        data=data,
        weights=weights,
        total_trials=4, n_workers=2, die=True
    )

    calib.calibrate(confirm_fit=True)

    print(f'Fit with original pars: {calib.before_fit}')
    print(f'Fit with best-fit pars: {calib.after_fit}')
    if calib.after_fit <= calib.before_fit:
        print(f'✓ Calibration improved fit ({calib.after_fit} <= {calib.before_fit})')
    else:
        print(f"✗ Calibration did not improve fit, but this isn't guaranteed ({calib.after_fit} > {calib.before_fit})")

    return sim, calib


#%% Run as a script
if __name__ == '__main__':

    T = sc.tic()

    sim, calib = test_calibration()

    sc.toc(T)
    print('Done.')
