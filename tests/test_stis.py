"""
STI scientific validation: verify that disease parameters affect epi dynamics
in the expected directions and that disease-specific features work correctly.

Covers all non-HIV STIs (syphilis, gonorrhea, chlamydia, trichomoniasis, BV, GUD)
and HIV epi parameter sensitivity. HIV-specific scientific validation (CD4, ART
effects, natural history) lives in test_hiv.py.
"""
import numpy as np
import pytest
import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti

n_agents = 500
debug = False


def test_syph_epi():
    """ Test that syphilis epi parameters affect dynamics in the expected direction """
    sc.heading('Test epi dynamics of syphilis')

    base_pars = dict(n_agents=n_agents)

    # Define the parameters to vary
    par_effects = dict(
        rel_trans_primary=[0.5, 1.5],  # Relative transmissibility during primary infection
        init_prev=[0.1, 0.9],
        beta_m2f=[0.2, 0.5]  # Beta for male to female transmission; opposite direction uses half this value
    )

    # Build all sims and run in one parallel call
    sims = []
    for par, (lo, hi) in par_effects.items():
        simpardict_lo = sc.mergedicts(dict(beta_m2f=0.3, init_prev=0.05), {par: lo})
        simpardict_hi = sc.mergedicts(dict(beta_m2f=0.3, init_prev=0.05), {par: hi})
        sims.append(sti.Sim(sc.mergedicts(base_pars, dict(diseases=sti.Syphilis(**simpardict_lo), dur=10)), label=f'{par}={lo}'))
        sims.append(sti.Sim(sc.mergedicts(base_pars, dict(diseases=sti.Syphilis(**simpardict_hi), dur=10)), label=f'{par}={hi}'))

    for sim in sims:
        sim.run()

    # Check results pairwise
    for i, (par, (lo, hi)) in enumerate(par_effects.items()):
        s0, s1 = sims[2*i], sims[2*i+1]
        ind = 1 if par == 'init_prev' else -1
        v0 = s0.results.syph.cum_infections[ind]
        v1 = s1.results.syph.cum_infections[ind]
        print(f'Checking with varying {par:10s} ... ', end='')
        assert v0 <= v1, f'Expected infections to be lower with {par}={lo} than with {par}={hi}, but {v0} > {v1})'
        print(f'✓ ({v0} <= {v1})')

    return sims


def test_syph_congenital_death_timing():
    """
    Regression test for #588: syphilis fetal-loss outcomes (miscarriage,
    stillborn) must be scheduled *before* ti_delivery so the baby's age is
    still <0 when the death fires. Pregnancy.process_prenatal_deaths only
    counts a death as a fetal loss when age<0 — age flips to 0 exactly at
    ti_delivery, so a death scheduled at ti_delivery is misclassified as a
    live birth. Also covers ti_miscarriage, which was assigned but never
    fired (dead code — no request_death, no result counter).
    """
    sc.heading('Test syphilis congenital death timing (#588)')

    def make_sim(outcome_idx):
        p = np.zeros(5)
        p[outcome_idx] = 1.0
        birth_outcomes = sc.objdict(mat_active=ss.choice(a=5, p=p), early=ss.choice(a=5, p=p), late=ss.choice(a=5, p=p))
        return sti.Sim(n_agents=n_agents, dur=12, start=2000,
                        diseases=sti.Syphilis(init_prev=0.3, beta_m2f=0.5, beta_m2c=1.0, birth_outcomes=birth_outcomes),
                        demographics=[ss.Pregnancy(), ss.Deaths()])

    # Force every congenital outcome to be 'stillborn' (index 2)
    sim_still = make_sim(2)
    sim_still.run()
    n_stillborn = sim_still.results.syph.new_stillborns.sum()
    n_pregnancy_losses = sim_still.results.pregnancy.stillbirths.sum() + sim_still.results.pregnancy.miscarriages.sum()
    assert n_stillborn > 0, 'Expected some congenital stillbirths with birth_outcomes forced to stillborn'
    # Every syph-flagged stillbirth should register as a Pregnancy fetal loss too;
    # a large gap means deaths are firing at/after ti_delivery and being
    # misclassified as live births instead.
    assert n_pregnancy_losses >= 0.8 * n_stillborn, (
        f'syph recorded {n_stillborn} stillbirths but Pregnancy recorded only {n_pregnancy_losses} fetal losses — '
        f'stillbirths are being misclassified as live births'
    )

    # Force every congenital outcome to be 'miscarriage' (index 0)
    sim_misc = make_sim(0)
    sim_misc.run()
    n_miscarriages = sim_misc.results.syph.new_miscarriages.sum()
    assert n_miscarriages > 0, 'ti_miscarriage is scheduled but never fires a death'

    return sim_still, sim_misc


def test_stis(which='discharging', n_agents=1e3, start=2010, stop=2020):
    """ Test running grouped STIs: discharging (NG, CT, TV, BV) and ulcerative (syphilis, GUD) """
    sc.heading('Test STI sim')

    if which == 'discharging':
        sti_pars = dict(
            ng=dict(beta_m2f=0.06, init_prev=0.02),
            ct=dict(beta_m2f=0.06, init_prev=0.05),
            tv=dict(beta_m2f=0.1, init_prev=0.1),
            bv=dict(),
        )
    elif which == 'ulcerative':
        sti_pars = dict(
            syph=dict(beta_m2f=0.1, init_prev=0.01),
            gud=dict(prevalence=0.05),
        )

    pregnancy = ss.Pregnancy(fertility_rate=10)
    death = ss.Deaths(death_rate=10)

    sim = sti.Sim(
        n_agents=n_agents,
        start=start,
        stop=stop,
        diseases=list(sti_pars.keys()),
        sti_pars=sti_pars,
        demographics=[pregnancy, death],
    )

    sim.run(verbose=1/12)

    return sim


def test_optional_stratification(n_agents=500):
    """ Test that age_bins=None and sex_keys=None drop stratified results without changing the rest """
    sc.heading('Test optional age/sex stratification of results')

    # Starsim 3 diseases: flat results by sex and age
    def make_sim(**kw):
        dis = [sti.HIV(init_prev=0.1, **kw)]
        return ss.Sim(n_agents=n_agents, start=2000, stop=2005, rand_seed=1, diseases=dis,
                      networks=[sti.StructuredSexual(), ss.MaternalNet()],
                      demographics=[ss.Pregnancy(fertility_rate=50), ss.Deaths(death_rate=10)])

    s_full = make_sim()
    s_none = make_sim(age_bins=None, sex_keys=None)
    s_full.run(); s_none.run()

    for name in ['hiv']:
        r_full = s_full.diseases[name].results
        r_none = s_none.diseases[name].results
        assert len(r_none) < len(r_full), f'{name}: expected fewer results without stratification'
        assert not [k for k in r_none.keys() if k.endswith(('_f', '_m')) or '_15_20' in k], \
            f'{name}: stratified results should not be defined'
        assert 'prevalence' in r_none and 'new_infections' in r_none, f'{name}: main results should remain'
        for key in r_none.keys():
            if key in r_full:
                assert np.allclose(np.asarray(r_full[key], float), np.asarray(r_none[key], float), equal_nan=True), \
                    f'{name}: result {key} changed when stratification was disabled'

    # Declared diseases: sti.Sim stratifies by sex and age, so results can be selected; sex_keys=None drops the flat ones
    kw = dict(n_agents=n_agents, start=2000, stop=2005, rand_seed=1)
    sim = sti.Sim(diseases=sti.Gonorrhea(init_prev=0.1), **kw).run()
    r = sim.results.ng
    assert np.allclose(r.n_infected_f + r.n_infected_m, r.n_infected)
    assert np.allclose(r.n_infected.by('age').values.sum(axis=1), r.n_infected)
    sim_none = sti.Sim(diseases=sti.Gonorrhea(init_prev=0.1, sex_keys=None), **kw).run()
    assert 'prevalence_f' not in sim_none.results.ng and np.array_equal(sim_none.results.ng.prevalence, r.prevalence)
    sim_noage = sti.Sim(diseases=sti.Gonorrhea(init_prev=0.1, age_bins=None), **kw).run()
    assert [s.name for s in sim_noage.strats] == ['sex']
    return s_full, s_none


def test_seis_natural_history(n_agents=5000):
    """ Test the bacterial STIs' natural history: symptoms by sex, care seeking, TV persistence, and PID """
    sc.heading('Test bacterial STI natural history')
    kw = dict(start=2000, stop=2004, dt=ss.months(1), n_agents=n_agents, stratify=dict(sex=['female', 'male']), verbose=0)
    no_infection = dict(beta=ss.permonth(0), init_prev=0.5)

    # Share symptomatic, by sex
    for cls in [sti.Gonorrhea, sti.Chlamydia, sti.Trichomoniasis]:
        sim = ss.Sim(cls(**no_infection), **kw).run()
        d = sim.diseases[0]
        r = d.results
        p_symp = d.pars.p_symp.value.values
        ppl = sim.people
        seeded = ~np.isnan(d.ti_E.raw[ppl.uid.raw]) # Everyone infected, all seeds since beta is 0
        symp = ~np.isnan(d.ti_P.raw[ppl.uid.raw])
        for female,p in zip([True, False], p_symp):
            sex = ppl.female.raw[ppl.uid.raw] == female
            share = (symp & sex).sum()/(seeded & sex).sum()
            assert np.isclose(share, p, atol=0.05), f'{d.name} female={female}: {share:.2f} symptomatic, expected {p}'
        assert r.new_care_seekers.sum() > 0

    # TV: most asymptomatic infections in women persist; men's clear
    sim = ss.Sim(sti.Trichomoniasis(**no_infection), **kw).run()
    r = sim.results.tv
    left = r.n_infected.select(sex='female')[-1]/r.n_infected.select(sex='female')[0], r.n_infected.select(sex='male')[-1]/r.n_infected.select(sex='male')[0]
    assert left[0] > 0.4 and left[1] < 0.05, f'Expected persistence in women only, not {left}'

    # PID: women only
    sim = ss.Sim(sti.Chlamydia(p_pid=0.3, **no_infection), **kw).run()
    r = sim.results.ct
    assert r.cum_pid_a.select(sex='female')[-1] > 0 and r.cum_pid_a.select(sex='male')[-1] == 0
    return sim


def test_abm_ode():
    """ Agents and an ODE agree for each bacterial STI under homogeneous mixing, at weekly steps (they differ by order dt) """
    sc.heading('Test ABM vs ODE agreement for NG, CT, and TV')
    pop = 20_000
    kw = dict(start=2000, stop=2012, dt=ss.weeks(1), stratify=dict(sex=['female', 'male']), verbose=0)
    out = sc.objdict()
    for cls,beta in [(sti.Gonorrhea, 0.16), (sti.Chlamydia, 0.09), (sti.Trichomoniasis, 0.06)]:
        res = sc.objdict()
        for solver,size,seeds in [('abm', dict(n_agents=pop), [0, 1]), ('ode', dict(total_pop=pop), [0])]:
            vals = []
            for seed in seeds:
                sim = ss.Sim(cls(beta=ss.permonth(beta), init_prev=0.05), solver=solver, rand_seed=seed, **size, **kw).run()
                r = sim.results[sim.diseases[0].name]
                vals.append([r.prevalence_f[-1], r.prevalence_m[-1], r.cum_infections[-1], r.new_care_seekers.sum()])
            res[solver] = np.mean(vals, axis=0)
        rel = np.abs(res.abm/res.ode - 1)
        out[cls.__name__] = sc.objdict(abm=res.abm, ode=res.ode, rel=rel)
        print(f'{cls.__name__}: ABM {np.round(res.abm, 3)}, ODE {np.round(res.ode, 3)}, differences {np.round(rel, 3)}')
        assert np.all(rel < 0.1), f'{cls.__name__}: ABM and ODE differ by {np.round(rel, 3)} (prevalence F, M; infections; care seekers)'
    return out


def test_bv(n_agents=5000):
    """ BV: women enter at 15 at their stable CST, risk factors raise the rate of BV, and treatment cures """
    sc.heading('Test BV')
    kw = dict(n_agents=n_agents, start=2000, stop=2010, verbose=0, rand_seed=1, demographics=[ss.Pregnancy(fertility_rate=ss.freqperyear(80)), ss.Deaths(death_rate=10)])
    sim = sti.Sim(diseases=sti.BV(), **kw).run()
    bv, ppl = sim.diseases.bv, sim.people
    r = sim.results.bv
    u = bv.U.uids
    assert not np.any(ppl.female[u] & (ppl.age[u] >= 15 + sim.t.dt.years)), 'Women 15+ are still without a CST' # Those who turned 15 in the last step move at the next
    assert not np.any(ppl.male[(~bv.U).uids]), 'Men have a CST'
    assert np.allclose(r.n_women.by('stable').values[-1]/r.n_women[-1], [0.1, 0.3, 0.6], atol=0.03)
    assert 0.3 < r.prevalence[-1] < 0.8 and np.isclose(r.prevalence[-1], r.cst4_prevalence[-1])
    assert np.allclose(r.cst1_prevalence + r.cst3_prevalence + r.cst4_prevalence + r.n_P/r.n_women, 1)
    assert r.clear_s.sum() + r.clear_a.sum() > 0.5*r.onset.sum() # BV clears (in STIsim 3, stable-CST4 girls never did)

    # Symptoms by stable CST
    cst4 = r.n_CST4p.by('stable').values.sum(axis=0) + r.n_CST4s.by('stable').values.sum(axis=0)
    share = cst4/(cst4 + r.n_CST4a.by('stable').values.sum(axis=0))
    assert share[0] > share[1] > share[2], f'Symptoms should be likeliest in stable CST1: {share}'

    # Risk factors: poor menstrual hygiene raises prevalence; without it, it is lower
    low = sti.Sim(diseases=sti.BV(p_poor_menstrual_hygiene=0.1), **kw).run().results.bv
    assert low.prevalence[-1] < r.prevalence[-1]
    hyg = r.n_infected.by('poor_menstrual_hygiene').values[-1]/r.n_women.by('poor_menstrual_hygiene').values[-1]
    assert hyg[1] > hyg[0], f'BV prevalence with and without poor hygiene: {hyg}'
    assert 'bv_uncirc_onset' in sim.what_modifies('bv.onset', output=True) # From the sexual network

    # Treatment: newly symptomatic women seek care, cured by stable CST; with a post-treatment window, recurrence
    tx = sti.Sim(diseases=sti.BV(), interventions=sti.treat_BV(), **kw).run()
    rt, t = tx.results.bv, tx.results.treat_bv
    assert np.isclose(t.new_treated.sum(), 0.9*rt.symptoms.sum(), rtol=0.05)
    assert np.isclose(t.cured.sum(), rt.tx.sum()) and np.isclose(t.new_doses.sum(), 14*t.new_treated.sum())
    assert rt.recur.sum() == 0 and np.isclose(rt.settle.sum(), rt.tx.sum(), rtol=0.01) # dur_post_tx=0: back at once, as in STIsim 3
    assert rt.prevalence[-1] < r.prevalence[-1]
    win = sti.Sim(diseases=sti.BV(dur_post_tx=ss.months(6)), interventions=sti.treat_BV(), **kw).run().results.bv
    assert win.recur.sum() > 0.2*win.tx.sum()
    with pytest.raises(ValueError):
        sti.BV(init_prev=0.1)
    with pytest.raises(NotImplementedError):
        sti.treat_BV(durable_cst1_change=True)
    return sim


def test_bv_abm_ode():
    """ BV as agents and as an ODE, with treatment: they agree (weekly steps; two risk factors, so the ODE is quick) """
    sc.heading('Test BV ABM vs ODE')
    pop = 20_000
    kw = dict(start=2000, stop=2008, dt=ss.weeks(1), stratify=[dict(sex=['female', 'male']), dict(age=[0, 15, 100])], verbose=0)
    res = sc.objdict()
    for solver,size,seeds in [('abm', dict(n_agents=pop), [0, 1]), ('ode', dict(total_pop=pop), [0])]:
        vals = []
        for seed in seeds:
            bv = sti.BV(rr_douching=1, rr_concurrency=1)
            sim = ss.Sim(diseases=bv, interventions=sti.treat_BV(), solver=solver, rand_seed=seed, **size, **kw).run()
            r = sim.results
            vals.append([r.bv.prevalence[-1], r.bv.cst1_prevalence[-1], r.bv.cst3_prevalence[-1], r.bv.onset.sum(), r.treat_bv.delivered.sum()])
        res[solver] = np.mean(vals, axis=0)
    rel = np.abs(res.abm/res.ode - 1)
    print(f'BV: ABM {np.round(res.abm, 3)}, ODE {np.round(res.ode, 3)}, differences {np.round(rel, 3)}')
    assert np.all(rel < 0.1), f'BV: ABM and ODE differ by {np.round(rel, 3)} (prevalence, CST1, CST3; onsets; treated)'
    return res


def test_gud():
    """ GUD: an SIS infection over the sexual network; none without seeding """
    sc.heading('Test GUD')
    kw = dict(n_agents=2000, start=2000, stop=2005, verbose=0)
    r = sti.Sim(diseases=sti.GUD(init_prev=0.1, beta_m2f=0.2), **kw).run().results.gud
    assert r.new_infections.sum() > 0 and r.recovery.sum() > 0 and 'prevalence_f' in r
    assert np.isclose(r.prevalence[0], 0.1, atol=0.03)
    r0 = sti.Sim(diseases=sti.GUD(), **kw).run().results.gud
    assert r0.n_infected.sum() == 0 # init_prev=0, as in STIsim 3
    return r


if __name__ == '__main__':
    sc.options(interactive=False)
    s1, s2 = test_syph_epi()
    test_syph_congenital_death_timing()
    test_stis(which='discharging')
    test_stis(which='ulcerative')
    test_optional_stratification()
    test_seis_natural_history()
    test_abm_ode()
    test_bv()
    test_bv_abm_ode()
    test_gud()
