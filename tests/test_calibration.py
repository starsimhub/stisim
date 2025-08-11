"""
Test calibration
"""

#%% Imports and settings
import starsim as ss
import sciris as sc
import stisim as sti
import pandas as pd

do_plot = 1
do_save = 0
n_agents = 2e3

# Settings
debug = True  # If True, this will do smaller runs that can be run locally for debugging
do_save = True


def make_sim():

    syph = sti.Syphilis(beta_m2f=0.06, init_prev=0.05)
    hiv = sti.HIV(beta_m2f=0.05, beta_m2c=0.025, init_prev=0.15)
    connector = sti.hiv_syph(hiv, syph, rel_sus_hiv_syph=2, rel_trans_hiv_syph=2)
    pregnancy = ss.Pregnancy(fertility_rate=20)
    death = ss.Deaths(death_rate=10)
    sexual = sti.StructuredSexual(prop_f2=0.2)
    maternal = ss.MaternalNet()

    sim = ss.Sim(
        dt=1,
        n_agents=n_agents,
        total_pop=9980999,
        start=1990,
        dur=33,
        diseases=[syph, hiv],
        networks=[sexual, maternal],
        connectors=connector,
        demographics=[pregnancy, death],
    )

    sim.init()

    return sim


def build_sim(sim, calib_pars):

    hiv = sim.diseases.hiv
    syph = sim.diseases.syphilis
    nw = sim.networks.structuredsexual

    # Apply the calibration parameters
    for k, pars in calib_pars.items():  # Loop over the calibration parameters
        if k == 'rand_seed':
            sim.pars.rand_seed = pars
            continue

        v = pars['value']
        if 'hiv_' in k:  # HIV parameters
            k = k.replace('hiv_', '')  # Strip off identifying part of parameter name
            hiv.pars[k] = v
        elif 'syph_' in k:  # Syphilis parameters
            k = k.replace('syph_', '')  # As above
            syph.pars[k] = v
        elif 'nw_' in k:  # Network parameters
            k = k.replace('nw_', '')  # As above
            if 'pair_form' in k:
                nw.pars[k].set(v)
            else:
                nw.pars[k] = v
        else:
            raise NotImplementedError(f'Parameter {k} not recognized')

    return sim


def run_calib(calib_pars=None):
    sc.heading('Beginning calibration')

    # Make the sim and data
    sim = make_sim()
    data = pd.read_csv('test_data/zimbabwe_calib.csv')

    # Make the calibration
    calib = sti.Calibration(
        sim=sim,
        calib_pars=calib_pars,
        build_fn=build_sim,
        total_trials=2,
        n_workers=1,
        die=True,
        reseed=False,
        debug=debug,
        data=data,
        save_results=True,
    )

    # Perform the calibration
    sc.printcyan('\nPeforming calibration...')
    calib.calibrate()
    calib.check_fit()
    calib.plot_optuna('plot_param_importances')

    sc.printcyan('\nShrinking calibration...')
    cal = calib.shrink()

    return sim, calib, cal


#%% Run as a script
if __name__ == '__main__':

    T = sc.tic()

    # Define the calibration parameters
    calib_pars = dict(
        hiv_beta_m2f = dict(low=0.01, high=0.10, guess=0.05),
        syph_beta_m2f = dict(low=0.01, high=0.10, guess=0.05),
        syph_rel_trans_latent = dict(low=0.00, high=0.2, guess=0.1),
        nw_prop_f1 = dict(low=0.1, high=0.45, guess=0.15),
        nw_prop_m1 = dict(low=0.15, high=0.5, guess=0.21),
        nw_f1_conc = dict(low=0.005, high=0.1, guess=0.01),
        nw_m1_conc = dict(low=0.005, high=0.1, guess=0.01),
        nw_p_pair_form = dict(low=0.4, high=0.9, guess=0.5),
    )

    sim, calib, cal = run_calib(calib_pars=calib_pars)

    sc.toc(T)
    print('Done.')




