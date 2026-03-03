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
from statistics import median, mean

tests_directory = Path(__file__).resolve().parent
sys.path.append(str(tests_directory))

from hiv_natural_history_analyzers import CD4ByUIDTracker, TimeToAIDSTracker
from testlib import build_testing_sim


verbose = False
do_plot = False
sc.options(interactive=False)


@sc.timer()
def test_cd4_counts_decline_over_time_without_treatment():
    sc.heading("Ensuring CD4 counts decline without treatment")

    sim = build_testing_sim(analyzers=[CD4ByUIDTracker()], n_agents=5, duration=5)
    sim.run()
    results = sim.results['cd4byuidtracker']['hiv.ts_cd4']

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
def test_median_time_from_infection_to_aids_without_treatment():
    sc.heading("Ensuring median time from infection to AIDS (to falling state) is median(acute_dur + latent_dur)")

    result_tolerance = 0.03  # fraction of the median
    sim = build_testing_sim(analyzers=[TimeToAIDSTracker()], n_agents=500, duration=5)
    sim.run()
    results = sim.results
    times_to_aids = list(chain(*results.timetoaidstracker['hiv.ti_to_aids']))

    # ensure we have at least ONE agent that progressed to AIDS before computing and checking median
    assert len(times_to_aids) > 0, "Failed to generate at least one HIV infection for testing"

    """
    import starsim as ss
    from statistics import median
    dur_acute = ss.lognorm_ex(3, 1)  # months
    dur_latent = ss.lognorm_ex(120, 36)  # months
    dur_acute.init(force=True)
    dur_latent.init(force=True)
    ti_falling = dur_acute.rvs(1000000) + dur_latent.rvs(1000000)
    expected_median = median(ti_falling)
    """
    median_time = median(times_to_aids)
    expected_median = 117.73883435759095  # months, from 1m random generations as above
    if verbose:
        print(f"{len(times_to_aids)} agents progressed to AIDS.")
        print(f"median time to AIDS stage: {median_time} timesteps.")
        print(f"expected_median: {expected_median} timesteps.")
        print(f"min time to AIDS stage: {min(times_to_aids)} timesteps.")
        print(f"mean time to AIDS: {mean(times_to_aids)} timesteps.")
        print(f"max time to AIDS stage: {max(times_to_aids)} timesteps.")
    delta = result_tolerance * expected_median
    msg = f"median time to AIDS: {median_time} not within delta: {delta} of expected: {expected_median}"
    assert median_time == pytest.approx(expected=expected_median, rel=result_tolerance), msg
    return sim


if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    timer = sc.timer()

    test_cd4_counts_decline_over_time_without_treatment()
    test_median_time_from_infection_to_aids_without_treatment()

    sc.heading("Total:")
    timer.toc()

    if do_plot:
        plt.show()