import starsim as ss
import stisim as sti
import unittest

from itertools import chain
from statistics import median, mean

from tests.hiv_natural_history_analyzers import CD4ByUIDTracker, RelativeInfectivityTracker, TimeToAIDSTracker, \
    TransmissionTracker, BreastfeedingTransmissionTracker
from tests.testlib import build_testing_sim

verbose = False


class TestHIVNaturalHistoryVerification(unittest.TestCase):

    def setUp(self):
        # default test setup; individual tests can replace/add to them.
        self.init_prev = 0.05
        self.beta_m2f = 0.05
        self.beta_m2c = 0.1
        self.diseases = [sti.HIV(beta_m2f=self.beta_m2f, beta_m2c=self.beta_m2c, init_prev=self.init_prev)]

        self.pregnancy = ss.Pregnancy(fertility_rate=10)
        self.death = ss.Deaths(death_rate=10)
        self.demographics = [self.pregnancy, self.death]

        self.sexual = sti.StructuredSexual(recall_prior=True)
        self.prior = sti.PriorPartners()
        self.maternal = ss.MaternalNet()
        self.networks = [self.sexual, self.prior, self.maternal]

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

    def test_falling_transmission_is_higher_than_latent(self):
        # since test_latent_transmission_ratio_is_always_1 ensures latent ratio is 1, we check > 1 here
        sim = build_testing_sim(diseases=self.diseases, demographics=self.demographics,
                                interventions=self.interventions, networks=self.networks,
                                analyzers=[RelativeInfectivityTracker(states=['falling'])],
                                n_agents=50, duration=5)
        sim.run()
        acute_ratios = sim.results['relativeinfectivitytracker']['hiv.falling_infectivity_ratios']
        acute_ratios = list(set(chain(*acute_ratios)))

        self.assertGreater(len(acute_ratios), 0)  # make sure we compare at least one item
        minimum_ratio = min(acute_ratios)
        self.assertGreater(minimum_ratio, 1)

    def test_no_transmission_if_no_sexual_network(self):
        networks = [self.prior, self.maternal]
        sim = build_testing_sim(diseases=self.diseases, demographics=self.demographics,
                                interventions=self.interventions, networks=networks,
                                analyzers=[TransmissionTracker()],
                                n_agents=50, duration=5)
        sim.run()
        hiv_transmisions = sum(sim.results['transmissiontracker']['hiv.n_transmissions'])
        self.assertEqual(0, hiv_transmisions)

    def test_breastfeeding_hiv_transmission_occurs(self):
        sim = build_testing_sim(diseases=self.diseases, demographics=self.demographics,
                                interventions=self.interventions, networks=self.networks,
                                analyzers=[BreastfeedingTransmissionTracker()],
                                n_agents=50, duration=10)
        sim.run()
        hiv_transmissions = sum(sim.results['breastfeedingtransmissiontracker']['hiv.n_bf_transmissions'])

        # ensure at least one such transmission occurs
        self.assertGreater(hiv_transmissions, 0)

    def _run_beta_test(self, baseline_m2f, baseline_m2c, multiplier=2, result_tolerance=0.1,
                       n_agents=50000, duration=1, init_prev=0.05, fertility=None):
        if fertility is None:
            demographics = self.demographics
        else:
            pregnancy = ss.Pregnancy(fertility_rate=fertility)
            demographics = [pregnancy, self.death]

        baseline_hiv = sti.HIV(beta_m2f=baseline_m2f, beta_m2c=baseline_m2c, init_prev=init_prev)
        sim = build_testing_sim(diseases=[baseline_hiv], demographics=demographics,
                                interventions=self.interventions, networks=self.networks,
                                analyzers=[TransmissionTracker()],
                                n_agents=n_agents, duration=duration)
        sim.run()
        hiv_transmissions_baseline = sum(sim.results['transmissiontracker']['hiv.n_transmissions'])

        # ensure at least one such transmission occurs
        self.assertGreater(hiv_transmissions_baseline, 0)

        # Now multiply betas as selected
        test_hiv = sti.HIV(beta_m2f=multiplier * baseline_m2f, beta_m2c=multiplier * baseline_m2c, init_prev=init_prev)
        sim = build_testing_sim(diseases=[test_hiv], demographics=demographics,
                                interventions=self.interventions, networks=self.networks,
                                analyzers=[TransmissionTracker()],
                                n_agents=n_agents, duration=duration)
        sim.run()
        hiv_transmissions_test = sum(sim.results['transmissiontracker']['hiv.n_transmissions'])

        # ensure at least one such transmission occurs
        self.assertGreater(hiv_transmissions_test, 0)

        # check transmission ratio, expected to be approximately increased by "multiplier" if base beta is low enough
        expected_ratio = multiplier
        delta = result_tolerance * expected_ratio
        test_ratio = hiv_transmissions_test / hiv_transmissions_baseline
        msg = f"baseline: {hiv_transmissions_baseline} test: {hiv_transmissions_test} ratio: {test_ratio} delta: {delta}"
        self.assertAlmostEqual(test_ratio, expected_ratio, delta=delta, msg=msg)
        if verbose:
            print(msg)

    def test_doubling_hiv_maternal_beta_doubles_relevant_transmissions(self):
        # setting baseline beta low, fertility HIGH, prevalence HIGH to generate enough births/transmissions quickly
        self._run_beta_test(baseline_m2f=0, baseline_m2c=0.0025, multiplier=2, fertility=1000,
                            duration=5, n_agents=3000, init_prev=1.0)

    def test_doubling_hiv_sexual_beta_doubles_relevant_transmissions(self):
        # doubling both beta for m2f (implicitly f2m) (sexual transmission only, no mother-to-child transmission)
        # This happens to be a realistic baseline m2f beta
        self._run_beta_test(baseline_m2f=0.001, baseline_m2c=0, multiplier=2)

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
