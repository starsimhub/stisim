import unittest

import starsim as ss
import stisim as sti

from tests.testlib import build_testing_sim


class TestHIVNaturalHistoryVerification(unittest.TestCase):

    def setUp(self):
        self.steps = 25

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
                                interventions=self.interventions, networks=self.networks)

        # sim.run_one_step() does not work for the first step due to an initialization bug:
        # https://github.com/starsimhub/starsim/issues/1136
        sim.run(until='1900-01-01')

        # now iterate for a while, recording CD4 counts of living people along the way (timeseries)
        cd4 = {}
        for step in range(self.steps):
            sim.run_one_step()
            for person in sim.people:
                if bool(person['alive']) is True and bool(person['hiv.infected']) is True:
                    uid = person['uid']
                    cd4_count = person['hiv.cd4']
                    cd4[uid] = cd4[uid] + [cd4_count] if uid in cd4 else [cd4_count]

        # Verify that CD4 counts cannot increase (stable and/or decrease) at all steps for all HIV+ individuals who
        # are not on treatment. No one is on treatment in this test.
        checks_performed = 0
        found_decrease = False
        for id, cd4_counts in cd4.items():
            for index in range(len(cd4_counts)-1):
                self.assertTrue(cd4_counts[index] >= cd4_counts[index+1])
            checks_performed += len(cd4_counts)-1

            # Ensure that at least ONE agent had a decreasing cd4 count (not just stable)
            if not found_decrease:
                if cd4_counts[0] > cd4_counts[-1]:
                    found_decrease = True

        self.assertGreater(checks_performed, 0)  # just in case, make sure we actually compared some CD4 counts
        self.assertTrue(found_decrease)
