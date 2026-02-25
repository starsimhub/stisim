import starsim as ss
import stisim as sti
import unittest

from itertools import chain
from statistics import median, mean

from tests.hiv_natural_history_analyzers import CD4ByUIDTracker, RelativeInfectivityTracker, TimeToAIDSTracker
from tests.testlib import build_testing_sim

verbose = False


class TestHIVNaturalHistoryVerification(unittest.TestCase):

    def setUp(self):
        # default test setup; individual tests can replace/add to them.
        self.diseases = [sti.HIV(beta_m2f=0.05, beta_m2c=0.1, init_prev=0.05)]

        pregnancy = ss.Pregnancy(fertility_rate=10)
        death = ss.Deaths(death_rate=10)
        self.demographics = [pregnancy, death]

        sexual = sti.StructuredSexual(recall_prior=True)
        prior = sti.PriorPartners()
        maternal = ss.MaternalNet()
        self.networks = [sexual, prior, maternal]

        self.interventions = []

    def test_cd4_counts_decline_over_time_without_treatment(self):
        sim = build_testing_sim(diseases=self.diseases, demographics=self.demographics,
                                interventions=self.interventions, networks=self.networks,
                                analyzers=[CD4ByUIDTracker()],
                                n_agents=5, duration=5)
        sim.run()
        results = sim.results['cd4byuidtracker']['hiv.ts_cd4']

        # ensure we have at least ONE agent with an actual decreasing CD4 (not just stable with the LessEqual assert)
        found_decrease = False
        for uid, cd4_timeseries in results.items():
            # compute change in cd4 at all timesteps for the agent
            cd4_deltas = [cd4_timeseries[i+1] - cd4_timeseries[i] for i in range(len(cd4_timeseries)-1)]
            minimum_delta = min(cd4_deltas)
            # assert CD4 can never go back up (untreated)
            self.assertLessEqual(minimum_delta, 0, msg=f"{minimum_delta} is not less than 0")
            if not found_decrease and minimum_delta < 0:
                found_decrease = True
        self.assertTrue(found_decrease)

    def test_median_time_from_infection_to_aids_without_treatment_using_run(self):
        result_tolerance = 0.03  # fraction of the median
        sim = build_testing_sim(diseases=self.diseases, demographics=self.demographics,
                                interventions=self.interventions, networks=self.networks,
                                analyzers=[TimeToAIDSTracker()],
                                n_agents=200, duration=25)
        sim.run()
        results = sim.results
        times_to_aids = list(chain(*results.timetoaidstracker['hiv.ti_to_aids']))

        # ensure we have at least ONE agent that progressed to AIDS before computing and checking median
        self.assertGreater(len(times_to_aids), 0)
        median_time = median(times_to_aids)
        if verbose:
            print(f"{len(times_to_aids)} agents progressed to AIDS.")

        expected_median = 117.89986811530153  # months, from 500k random generations from time-to-falling code.
        delta = result_tolerance * expected_median
        if verbose:
            print(f"median time to AIDS stage: {median_time} timesteps.")
            print(f"expected_median: {expected_median} timesteps.")
            print(f"min time to AIDS stage: {min(times_to_aids)} timesteps.")
            print(f"mean time to AIDS: {mean(times_to_aids)} timesteps.")
            print(f"max time to AIDS stage: {max(times_to_aids)} timesteps.")
        self.assertAlmostEqual(median_time, expected_median, delta=delta)

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
