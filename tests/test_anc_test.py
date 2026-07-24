"""
Tests for ANCTest (antenatal/infant PMTCT screening, issue #324).

Prior to this test file, ANCTest had no test coverage at all, and was
broken end-to-end: `init_pre` used the `'month'` alias with `ss.dur`
(crashes under current starsim) and computed `active_diseases` *after*
calling `super().init_pre(sim)`, which internally calls `init_results()` --
so the per-disease results were never created.
"""

import starsim as ss
import stisim as sti


def _make_sim(visit_prob=1.0, n_agents=2000, dur=8, rand_seed=1, art_coverage=0.6):
    infant_hiv = sti.InfantHIVTest(name='infant_hiv')
    anc = sti.ANCTest(visit_prob=visit_prob, newborn_tests={'hiv': infant_hiv})
    art = sti.ART(coverage=art_coverage)
    sim = sti.Sim(
        diseases=[sti.HIV(init_prev=0.3, beta_m2f=0.05, beta_m2c=0.1)],
        demographics=[ss.Pregnancy(fertility_rate=50), ss.Deaths()],
        networks=[sti.StructuredSexual(), ss.MaternalNet(), ss.BreastfeedingNet()],
        interventions=[anc, infant_hiv, art],
        n_agents=n_agents, dur=dur, start=2000, verbose=0, rand_seed=rand_seed,
    )
    sim.run()
    return sim


def test_anctest_runs_and_detects_hiv():
    """ ANCTest should run without crashing and record HIV-positive ANC attendees """
    sim = _make_sim()
    anc_res = sim.results.anctest
    assert anc_res.n_attended.sum() > 0
    assert anc_res.n_hiv_positive.sum() > 0


def test_anctest_schedules_immediate_art():
    """An HIV+ agent attending ANC should be diagnosed and have ART scheduled
    for the same timestep (no dx-to-tx delay). Uses a single manually-forced
    ANC visit rather than a full stochastic run, since ART flips `on_art` to
    True within the same step it's scheduled -- a full-sim analyzer can't
    distinguish "just scheduled" from "already treated" after the fact.
    """
    infant_hiv = sti.InfantHIVTest(name='infant_hiv')
    anc = sti.ANCTest(visit_prob=1.0, newborn_tests={'hiv': infant_hiv})
    art = sti.ART(coverage=0.6)
    sim = sti.Sim(
        diseases=[sti.HIV(init_prev=0)],
        demographics=[ss.Pregnancy(fertility_rate=50), ss.Deaths()],
        networks=[sti.StructuredSexual(), ss.MaternalNet(), ss.BreastfeedingNet()],
        interventions=[anc, infant_hiv, art],
        n_agents=500, dur=1, start=2000, verbose=0, rand_seed=1,
    )
    sim.init()

    hiv = sim.diseases.hiv
    uid = sim.people.female.uids[0]
    hiv.susceptible[uid] = False
    hiv.infected[uid] = True
    hiv.on_art[uid] = False
    sim.interventions.anctest.ti_visit[uid] = sim.ti  # force an ANC visit this step

    sim.interventions.anctest.step()

    assert hiv.diagnosed[uid]
    assert hiv.ti_diagnosed[uid] == sim.ti
    assert hiv.ti_art[uid] == sim.ti


def test_anctest_schedules_newborn_test():
    """ Newborns of HIV+ mothers diagnosed at ANC should get an infant HIV test """
    sim = _make_sim()
    assert sim.results.infant_hiv.new_tests.sum() > 0


if __name__ == '__main__':
    test_anctest_runs_and_detects_hiv()
    test_anctest_schedules_immediate_art()
    test_anctest_schedules_newborn_test()
