"""
HIV natural history verification tests

Tests to ensure appropriate behavior of HIV as a disease absent any treatment.
"""

import matplotlib.pyplot as plt
import pytest
import sciris as sc
import sys

from itertools import chain
from pathlib import Path
from statistics import mean

import starsim as ss
import stisim as sti

from stisim import ART, HIVTest

tests_directory = Path(__file__).resolve().parent
sys.path.append(str(tests_directory))

from hiv_natural_history_analyzers import CD4ByUIDTracker, RelativeInfectivityTracker, TimeToAIDSTracker, \
    SexualTransmissionCountTracker, MTCTransmissionCountTracker

from testlib import build_testing_sim


verbose = False
do_plot = False
sc.options(interactive=False)


@sc.timer()
def test_cd4_counts_decline_untreated():
    sc.heading("Ensuring CD4 counts decline without treatment")

    analyzer = CD4ByUIDTracker(subpop=CD4ByUIDTracker.INFECTED)
    sim = build_testing_sim(analyzers=[analyzer],
                            maternal_network=None, prior_network=None, sexual_network=None,
                            pregnancy=None, death=None,
                            n_agents=5, duration=5)
    sim.run()
    results = sim.results[analyzer.name][analyzer.result_name]

    assert len(results) > 0, f"Test requires at least one agent to have HIV, 0 found."

    found_decrease = False
    for uid, cd4_timeseries in results.items():
        # compute change in cd4 at all timesteps for the agent
        cd4_deltas = [cd4_timeseries[i+1] - cd4_timeseries[i] for i in range(len(cd4_timeseries)-1)]
        minimum_delta = min(cd4_deltas)

        # assert CD4 can never go back up (untreated)
        assert minimum_delta <= 0, f"CD4 erroneously went up for at least one untreated agent"
        if not found_decrease and minimum_delta < 0:
            found_decrease = True
    # ensure we have at least ONE agent with an actual decreasing CD4 (not just stable)
    assert found_decrease is True, f"CD4 is stable for all untreated agents, no decrease observed"
    return sim


@sc.timer()
def test_time_from_infection_to_aids_untreated():
    """Note: if a different definition of AIDS is represented in HIV module, update to use that."""
    sc.heading("Regression: Ensuring mean time from infection to AIDS (to falling state) is reasonable.")

    result_tolerance = 0.03  # fraction of the expected value
    sim = build_testing_sim(analyzers=[TimeToAIDSTracker()], n_agents=2000, duration=5)
    sim.run()
    results = sim.results
    times_to_aids = list(chain(*results.timetoaidstracker['hiv.ti_to_aids']))

    # ensure we have at least ONE agent that progressed to AIDS before computing and checking mean
    assert len(times_to_aids) > 0, "Failed to generate at least one HIV infection for testing"

    mean_time = mean(times_to_aids)
    expected_mean = 121.02083587646484  # months, from 10k agents, 25 years 5% prevalence
    if verbose:
        print(f"{len(times_to_aids)} agents progressed to AIDS.")
        print(f"mean time to AIDS stage: {mean_time} timesteps.")
        print(f"expected_mean: {expected_mean} timesteps.")
        print(f"min time to AIDS stage: {min(times_to_aids)} timesteps.")
        print(f"max time to AIDS stage: {max(times_to_aids)} timesteps.")
    delta = result_tolerance * expected_mean
    msg = f"mean time to AIDS: {mean_time} not within delta: {delta} of expected: {expected_mean}"
    assert mean_time == pytest.approx(expected=expected_mean, rel=result_tolerance), msg
    return sim


@sc.timer()
def test_latent_transmission_ratio_is_1():
    sc.heading("Ensuring latent HIV transmission ratio is always 1.")

    sim = build_testing_sim(analyzers=[RelativeInfectivityTracker(states=['latent'])], n_agents=100, duration=1)
    sim.run()
    latent_ratios = sim.results['relativeinfectivitytracker']['hiv.latent_rel_trans']
    latent_ratios = list(set(chain(*latent_ratios)))

    assert len(latent_ratios) > 0, "Cannot test latent HIV transmission ratios, no latent transmission detected."
    assert len(latent_ratios) == 1, f"Multiple latent HIV transmission ratios detected (only 1 is valid): {latent_ratios}"
    assert latent_ratios[0] == 1, f"Only one latent HIV transmission ratio detected, but it is not 1: {latent_ratios}"
    return sim


@sc.timer()
def test_acute_transmission_higher_than_latent():
    sc.heading("Checking HIV transmission ratio acute > latent.")

    sim = build_testing_sim(analyzers=[RelativeInfectivityTracker(states=['acute'])], n_agents=500, duration=1)
    sim.run()
    acute_ratios = sim.results['relativeinfectivitytracker']['hiv.acute_rel_trans']
    acute_ratios = list(set(chain(*acute_ratios)))

    assert len(acute_ratios) > 0, "Cannot test acute HIV transmission ratios, no acute transmission detected."
    minimum_ratio = min(acute_ratios)
    assert minimum_ratio > 1, f"Acute HIV transmission ratios must be greater than 1. Lowest detected value: {minimum_ratio}"
    return sim


@sc.timer()
def test_aids_transmission_is_higher_than_latent():
    sc.heading("Checking HIV transmission ratio AIDS > latent.")

    sim = build_testing_sim(analyzers=[RelativeInfectivityTracker(states=['aids'])], n_agents=25, duration=10)
    sim.run()
    aids_ratios = sim.results['relativeinfectivitytracker']['hiv.aids_rel_trans']
    aids_ratios = list(set(chain(*aids_ratios)))

    assert len(aids_ratios) > 0, "Cannot test AIDS HIV transmission ratios, no AIDS transmission detected."
    minimum_ratio = min(aids_ratios)
    assert minimum_ratio > 1, f"AIDS HIV transmission ratios must always be > 1. Minimum detected value: {minimum_ratio}"
    return sim


@sc.timer()
def test_no_sexual_transmission_without_network():
    sc.heading("Ensuring no sexual transmission if there is no sexual network.")

    analyzer = SexualTransmissionCountTracker()
    sim = build_testing_sim(analyzers=[analyzer], sexual_network=None, n_agents=100, duration=1)
    sim.run()
    n_hiv_transmissions = sum(sim.results[analyzer.name][analyzer.result_name])

    assert n_hiv_transmissions == 0, f"Expected no sexual transmissions, but {n_hiv_transmissions} were detected."
    return sim


def _run_beta_test(baseline_m2f, baseline_m2c, mode: str, multiplier=2, result_tolerance=0.16,
                   n_agents=50000, duration=1, init_prev=0.2, fertility=10, rand_seed=1):
    if mode == 'sexual':
        analyzer = SexualTransmissionCountTracker()
    elif mode == 'mtc':
        analyzer = MTCTransmissionCountTracker()
    else:
        raise ValueError(f'Unknown beta test mode: {mode}')

    pregnancy = ss.Pregnancy(fertility_rate=fertility)

    # run baseline sim
    baseline_hiv = sti.HIV(beta_m2f=baseline_m2f, beta_m2c=baseline_m2c, init_prev=init_prev)
    sim = build_testing_sim(diseases=[baseline_hiv], pregnancy=pregnancy,
                            analyzers=[analyzer],
                            n_agents=n_agents, duration=duration)
    sim.pars['rand_seed'] = rand_seed
    sim.run()

    hiv_transmissions_baseline = sum(sim.results[analyzer.name][analyzer.result_name])

    # ensure at least one such transmission occurs
    assert hiv_transmissions_baseline > 0, f"Cannot assess effect of beta, no {mode} HIV transmissions occurred in baseline."

    # Now multiply betas as selected and run test
    test_hiv = sti.HIV(beta_m2f=multiplier * baseline_m2f, beta_m2c=multiplier * baseline_m2c, init_prev=init_prev)
    sim = build_testing_sim(diseases=[test_hiv], pregnancy=pregnancy,
                            analyzers=[analyzer],
                            n_agents=n_agents, duration=duration)
    sim.pars['rand_seed'] = rand_seed
    sim.run()

    hiv_transmissions_test = sum(sim.results[analyzer.name][analyzer.result_name])

    # ensure at least one such transmission occurs
    assert hiv_transmissions_test > 0, f"Cannot assess effect of beta, no {mode} HIV transmissions occurred in test."

    # check transmission ratio, expected to be approximately increased by "multiplier" if base beta is low enough
    expected_ratio = multiplier
    delta = result_tolerance * expected_ratio
    test_ratio = hiv_transmissions_test / hiv_transmissions_baseline

    if verbose:
        msg = f"{mode} transmission ratio: {test_ratio} ({hiv_transmissions_test}/{hiv_transmissions_baseline}) expected: {expected_ratio}"
        print(msg)

    msg = f"{mode} transmission ratio: {test_ratio} not within delta: {delta} of expected: {expected_ratio}"
    assert test_ratio == pytest.approx(expected=expected_ratio, rel=result_tolerance), msg



@sc.timer()
def test_doubling_hiv_maternal_beta_doubles_transmissions():
    sc.heading("Checking that doubling mtc beta roughly doubles mtc transmissions.")

    # setting baseline beta low, fertility HIGH, prevalence HIGH to generate enough births/transmissions quickly
    _run_beta_test(baseline_m2f=0, baseline_m2c=0.0025, mode='mtc', multiplier=2, fertility=1000,
                   duration=5, n_agents=40000, init_prev=1.0)


@sc.timer()
def test_doubling_hiv_sexual_beta_doubles_transmissions():
    sc.heading("Checking that doubling sexual beta roughly doubles sexual transmissions at low infectivity.")

    # doubling both beta for m2f (implicitly f2m) (sexual transmission only, no mother-to-child transmission)
    # This happens to be a realistic baseline m2f beta
    _run_beta_test(baseline_m2f=0.001, baseline_m2c=0, mode='sexual', multiplier=2, duration=1, n_agents=100000)


@sc.timer()
def test_mtc_transmission_occurs():
    sc.heading("Ensuring that pre-term mother-to-child transmission occurs.")

    # setting fertility rate super high to enable shrinking the test agent/timewise
    analyzer = MTCTransmissionCountTracker()
    sim = build_testing_sim(analyzers=[analyzer], n_agents=500, duration=1, pregnancy=ss.Pregnancy(fertility_rate=1000))
    sim.run()

    mtc_tranmissions = sim.results[analyzer.name][analyzer.result_name]
    total_mtc_tranmissions = sum(mtc_tranmissions)

    if verbose:
        print(f"{total_mtc_tranmissions} mother-to-child transmissions were recorded")
    assert total_mtc_tranmissions > 0, f"Expected MTC transmissions to occur, but none were recorded."
    return sim


@sc.timer()
def test_cd4_rises_on_ART():
    sc.heading("Ensuring that agent CD4 levels rise when on ART (monotonically in a short test)")

    # To keep test small, infect everyone with HIV and put everyone on ART immediately
    analyzer = CD4ByUIDTracker(subpop=CD4ByUIDTracker.ONART)
    test_intervention = HIVTest(test_prob_data=1.0, dt_scale=False)  # everyone tests, first timestep
    initial_art_intervention = ART(art_initiation=1.0)  # everyone diagnosed starts ART.

    duration = 1  # years
    hiv = sti.HIV(beta_m2f=0.05, beta_m2c=0.1, init_prev=1.0, dur_on_art=ss.constant(v=ss.years(duration)))
    sim = build_testing_sim(analyzers=[analyzer], diseases=[hiv],
                            death=None, maternal_network=None, prior_network=None, sexual_network=None,
                            interventions=[test_intervention, initial_art_intervention],
                            n_agents=5, duration=duration)
    sim.run()
    cd4_ts_by_uid = sim.results[analyzer.name][analyzer.result_name]

    # Ensure at least one agent timeseries was recorded
    assert len(cd4_ts_by_uid.keys()) > 0

    # assert that cd4 cound must go up monotonically while on ART (it eventually maxes out (delta=0) after a long time)
    for uid, cd4_ts in cd4_ts_by_uid.items():
        delta_cd4 = [cd4_ts[i+1] - cd4_ts[i] for i in range(len(cd4_ts)-1)]
        for delta in delta_cd4:
            assert delta > 0, f"Expected CD4 count on ART to go up, but instead changed by: {delta}"
    return sim


@sc.timer()
def test_art_increases_longevity():
    sc.heading("Ensuring that untreated HIV always kills and ART can prevent this (to some degree).")

    test_intervention = HIVTest(test_prob_data=1, dt_scale=False)  # everyone tests first timestep
    initial_art_intervention = ART(art_initiation=0.5)  # 50% of diagnosed agents uptake ART
    interventions = [test_intervention, initial_art_intervention]

    # infect everyone immediately
    # agents are on art for the full sim length
    #   ... and more (due to bug: https://github.com/starsimhub/stisim/issues/336)
    duration = 20  # years by default
    disease = sti.HIV(init_prev=1.0, dur_on_art=ss.constant(v=ss.years(duration*100)))

    sim = build_testing_sim(n_agents=10, duration=duration, pregnancy=None, death=None,
                            interventions=interventions, diseases=[disease])
    sim.run()

    # at the end of the sim, assert that all agents are on ART (no non-ART agents alive)
    # at the end of the sim, assert that at least one person is alive
    hiv = sim.diseases.hiv
    n_alive = len(sim.people.alive.uids)
    n_alive_not_on_art = len((sim.people.alive & ~hiv.on_art).uids)

    if verbose:
        print(f"{n_alive} agents are alive, "
              f"{n_alive-n_alive_not_on_art} agents are alive and on ART, "
              f"{n_alive_not_on_art} are alive and NOT on ART.")

    assert n_alive > 0, f"Expected at least one agent to still be living, found none."
    assert n_alive_not_on_art == 0, f"Expected no agents to be still living without ART, there are {n_alive_not_on_art} ."

    return sim


@sc.timer()
def test_no_hiv_with_no_outbreaks():
    sc.heading("Ensuring that HIV infections remains zero without any seeding infections/events.")

    disease = sti.HIV(beta_m2f=0.05, beta_m2c=0.1, init_prev=0)
    sim = build_testing_sim(n_agents=1000, duration=3, diseases=[disease])  # , analyzers=[analyzer])
    sim.run()

    # HIV infections should be 0 at all timesteps
    n_infected = sim.results.hiv.n_infected
    unique_values = list(set(n_infected))

    assert len(unique_values) == 1, f"Found {len(unique_values)} unique HIV infection counts, but there should only be one."
    assert unique_values[0] == 0, f"Found a single unique infection count: {unique_values[0]}, but it is not 0 as expected."

    return sim


@sc.timer()
def test_increased_testing_speeds_diagnosis():
    sc.heading("Ensuring that increased testing diagnoses agents more quickly (short test to prevent conflating HIV deaths)")

    base_prob = 0.05
    base_testing = HIVTest(test_prob_data=base_prob, dt_scale=False)
    higher_prob = base_prob * 2
    higher_testing = HIVTest(test_prob_data=higher_prob, dt_scale=False)

    n_agents = 100
    duration = 1  # years by default
    disease = sti.HIV(init_prev=1.0)

    sim1 = build_testing_sim(n_agents=n_agents, duration=duration,
                            diseases=[disease],
                            pregnancy=None, death=None,
                            prior_network=None, sexual_network=None, maternal_network=None,
                            interventions=[base_testing])
    sim2 = build_testing_sim(n_agents=n_agents, duration=duration,
                            diseases=[disease],
                            pregnancy=None, death=None,
                            prior_network=None, sexual_network=None, maternal_network=None,
                            interventions=[higher_testing])
    sims = [sim1, sim2]
    msim = ss.parallel(*sims)

    hiv = sim1.diseases.hiv
    n_diagnosed_base = hiv.results.n_diagnosed[-1]  # total diagnoses by end of simulation

    hiv = sim2.diseases.hiv
    n_diagnosed_higher = hiv.results.n_diagnosed[-1]  # total diagnoses by end of simulation

    if True:
        print(f"Testing rate: {base_prob}/dt n_diagnosed: {n_diagnosed_base} after {duration} year(s)\n"
              f"Testing rate: {higher_prob}/dt n_diagnosed: {n_diagnosed_higher} after {duration} year(s)\n")

    assert n_diagnosed_base > 0, "Expected at least one agent to be diagnosed, found none."
    assert n_diagnosed_higher > 0, "Expected at least one agent to be diagnosed, found none."
    assert n_diagnosed_higher > n_diagnosed_base, (f"Expected higher testing rates to diagnose agents faster, "
                                                   f"but went from {n_diagnosed_base} to {n_diagnosed_higher} diagnoses in {duration} year(s) instead.")

    return sims


# Not currently implemented in hivsim, so leaving this partially-completed test commented out for future work
# def test_perinatally_infected_progress_faster(self):
#     sim = build_testing_sim(diseases=self.diseases, demographics=self.demographics,
#                             interventions=self.interventions, networks=self.networks,
#                             analyzers=[BirthTracker()],
#                             n_agents=500, duration=30)
#     sim.run()
#
#     peri_infected_uids = sim.results['birthtracker']['hiv.perinatally_infected_uids']
#     tis_acute = sim.results['birthtracker']['hiv.tis_acute']
#     tis_latent = sim.results['birthtracker']['hiv.tis_latent']
#     tis_falling = sim.results['birthtracker']['hiv.tis_falling']


if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    timer = sc.timer()

    test_cd4_counts_decline_untreated()
    test_time_from_infection_to_aids_untreated()
    test_latent_transmission_ratio_is_1()
    test_acute_transmission_higher_than_latent()
    test_aids_transmission_is_higher_than_latent()
    test_no_sexual_transmission_without_network()
    test_doubling_hiv_maternal_beta_doubles_transmissions()
    test_doubling_hiv_sexual_beta_doubles_transmissions()
    test_mtc_transmission_occurs()
    test_cd4_rises_on_ART()
    test_art_increases_longevity()
    test_no_hiv_with_no_outbreaks()
    test_increased_testing_speeds_diagnosis()

    sc.heading("Total:")
    timer.toc()

    if do_plot:
        plt.show()
