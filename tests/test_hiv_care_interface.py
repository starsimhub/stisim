"""
What HIV gives care outside the HIV interventions (e.g. a cascade of care that decides itself who is tested,
treated and protected): WHO stage, first symptoms, PMTCT for mothers not on ART, ART in an arm the caller chooses,
and PrEP courses started without a Prep intervention.
"""
import numpy as np
import sciris as sc
import starsim as ss
import stisim as sti
import hivsim

sc.options(interactive=False)


class Every(ss.Intervention):
    """ Calls `fn(sim)` every step """
    def __init__(self, fn, name='every'):
        super().__init__(name=name)
        self.fn = fn

    def step(self):
        self.fn(self.sim)


def test_who_stage_and_symptoms():
    """ First symptoms come once, as the CD4 count falls below cd4_symptoms, when WHO stage is 4 """
    seen = []
    def record(sim):
        hiv = sim.diseases.hiv
        new = hiv.new_symptoms()
        seen.append(new)
        assert np.all(hiv.cd4[new] < hiv.pars.cd4_symptoms)
        assert np.all(hiv.who_stage(new) == 4)
        uninfected = (~hiv.infected & sim.people.alive).uids
        assert np.all(hiv.who_stage(uninfected) == 0)
    hiv = sti.HIV(init_prev=1.0, beta_m2f=0, beta_m2c=0, beta_breastfeed=0)
    sim = ss.Sim(dt=1/12, n_agents=300, start=2000, dur=20, diseases=hiv, interventions=Every(record), rand_seed=1)
    sim.run()
    uids = np.concatenate(seen)
    assert len(uids) > 50, 'Expected people to become symptomatic over 20 years untreated'
    assert len(np.unique(uids)) == len(uids), 'First symptoms should come once per person'
    stages = sim.diseases.hiv.who_stage(sim.people.alive.uids)
    assert set(np.unique(stages)) <= {1, 2, 3, 4}
    return sim


def test_start_pmtct():
    """ PMTCT given to mothers not on ART cuts transmission in pregnancy """
    def make(eff):
        def give(sim):
            hiv, preg = sim.diseases.hiv, sim.demographics.pregnancy
            mothers = (preg.P & hiv.infected).uids
            if eff:
                hiv.start_pmtct(mothers, eff=eff, dur=ss.years(1))
        return hivsim.Sim(n_agents=1000, dur=10, init_prev=0.3, beta_m2c=0.1, rand_seed=2, interventions=Every(give), label=str(eff))
    sims = ss.parallel(make(0), make(1.0)).sims
    none, full = [s.results.hiv.new_infections_prenatal.sum() for s in sims]
    assert none > 10, f'Expected prenatal infections without PMTCT, got {none}'
    assert full == 0, f'Expected no prenatal infections with full PMTCT, got {full}'
    return sims


def test_start_art_in_chosen_arm():
    """ ART started in the arm the caller chooses, without changing ART's own probability of suppression """
    chosen = {}
    def start(sim):
        hiv, art = sim.diseases.hiv, sim.interventions.art
        if sim.ti == 12:
            uids = (hiv.infected & ~hiv.on_art).uids
            arm = np.arange(len(uids)) % 2 == 0
            art.start_art(uids, effective=arm)
            chosen.update(uids=uids, arm=arm)
            assert np.array_equal(art.on_effective_art[uids], arm)
            assert np.array_equal(art.on_nonsuppressive_art[uids], ~arm)
    art = sti.ART(p_effective_art=ss.bernoulli(p=0.3), dur_on_art=ss.constant(ss.years(100)))
    sim = ss.Sim(dt=1/12, n_agents=200, start=2000, dur=3, diseases=sti.HIV(init_prev=0.5, beta_m2f=0, beta_m2c=0, beta_breastfeed=0), interventions=[Every(start), art], rand_seed=3)
    sim.run()
    assert len(chosen['uids']) > 20
    assert art.pars.p_effective_art.pars.p == 0.3
    return sim


def test_prep_without_intervention():
    """ A PrEP course started by other care (source 0) ends when due, by HIV """
    started = {}
    def start(sim):
        hiv = sim.diseases.hiv
        if sim.ti == 6:
            uids = (~hiv.infected & sim.people.alive).uids[:50]
            hiv.start_prep(uids, eff=0.75, dur=ss.months(3))
            started.update(uids=uids, ti=sim.ti)
        if 'uids' in started:
            on = hiv.on_prep[started['uids']] | ~sim.people.alive[started['uids']]
            assert on.all() if sim.ti < started['ti'] + 3 else not hiv.on_prep[started['uids']].any()
    sim = ss.Sim(dt=1/12, n_agents=200, start=2000, dur=2, diseases=sti.HIV(init_prev=0.1, beta_m2f=0, beta_m2c=0, beta_breastfeed=0), interventions=Every(start), rand_seed=4)
    sim.run()
    return sim


if __name__ == '__main__':
    test_who_stage_and_symptoms()
    test_start_pmtct()
    test_start_art_in_chosen_arm()
    test_prep_without_intervention()
