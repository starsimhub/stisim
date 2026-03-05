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

tests_directory = Path(__file__).resolve().parent
sys.path.append(str(tests_directory))

from hiv_natural_history_analyzers import CD4ByUIDTracker, RelativeInfectivityTracker, TimeToAIDSTracker
from testlib import build_testing_sim


verbose = False
do_plot = False
sc.options(interactive=False)


@sc.timer()
def test_cd4_counts_decline_untreated():
    sc.heading("Ensuring CD4 counts decline without treatment")

    sim = build_testing_sim(analyzers=[CD4ByUIDTracker()],
                            maternal_network=None, prior_network=None, sexual_network=None,
                            pregnancy=None, death=None,
                            n_agents=5, duration=5)
    sim.run()
    results = sim.results['cd4byuidtracker']['hiv.ts_cd4']

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


def test_latent_transmission_ratio_is_always_1(self):
    sim = build_testing_sim(diseases=self.diseases, demographics=self.demographics,
                            interventions=self.interventions, networks=self.networks,
                            analyzers=[RelativeInfectivityTracker(states=['latent'])],
                            n_agents=5, duration=5)
    sim.run()
    latent_ratios = sim.results['relativeinfectivitytracker']['hiv.latent_infectivity_ratios']
    latent_ratios = list(set(chain(*latent_ratios)))

    self.assertEqual(len(latent_ratios), 1)
    self.assertEqual(latent_ratios[0], 1)


def test_acute_transmission_is_higher_than_latent(self):
    # since test_latent_transmission_ratio_is_always_1 ensures latent ratio is 1, we check > 1 here
    sim = build_testing_sim(diseases=self.diseases, demographics=self.demographics,
                            interventions=self.interventions, networks=self.networks,
                            analyzers=[RelativeInfectivityTracker(states=['acute'])],
                            n_agents=50, duration=5)
    sim.run()
    acute_ratios = sim.results['relativeinfectivitytracker']['hiv.acute_infectivity_ratios']
    acute_ratios = list(set(chain(*acute_ratios)))

    self.assertGreater(len(acute_ratios), 0)  # make sure we compare at least one item
    minimum_ratio = min(acute_ratios)
    self.assertGreater(minimum_ratio, 1)

if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    timer = sc.timer()

    test_cd4_counts_decline_untreated()
    test_time_from_infection_to_aids_untreated()

    sc.heading("Total:")
    timer.toc()

    if do_plot:
        plt.show()
