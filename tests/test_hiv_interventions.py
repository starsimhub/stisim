"""
Test HIV interventions: ART, VMMC, PrEP

Tests cover:
    - ART coverage input format equivalence (scalar, dict, DataFrame, stratified)
    - ART functional effects (people go on ART, warning without testing)
    - Scientific validation skeletons (mortality, parameter sensitivity, duration)
"""

import warnings
import sciris as sc
import numpy as np
import pandas as pd
import starsim as ss
import stisim as sti
import hivsim
import matplotlib.pyplot as pl

n_agents = 1_000
do_plot  = False
sc.options(interactive=False)


# %% Coverage format tests
@sc.timer()
def test_art_specs(do_plot=do_plot):
    """ Check that scalar, dict, and DataFrame coverage all produce ART uptake """
    sc.heading('Testing ART coverage format equivalence...')

    years = np.arange(2000, 2021)
    p_vals = np.linspace(0, 0.9, len(years))

    specs = dict(
        scalar = sti.ART(coverage=0.5),
        dict_p = sti.ART(coverage={'year': [2000, 2020], 'value': [0, 0.9]}),
        df_p   = sti.ART(coverage=pd.DataFrame(index=years, data={'p_art': p_vals})),
        df_n   = sti.ART(coverage=pd.DataFrame(index=years, data={'n_art': p_vals * 1e6})),
    )

    sims = []
    for label, art in specs.items():
        sim = hivsim.demo('simple', run=False, plot=False, n_agents=n_agents, dur=5)
        sim.pars.interventions = [sti.HIVTest(name='hiv_test', test_prob_data=0.3), art]
        sim.label = label
        sims.append(sim)

    msim = ss.parallel(sims)

    results = dict()
    for sim in msim.sims:
        n_art = sim.results.hiv.n_on_art[-1]
        results[sim.label] = n_art
        assert n_art > 0, f'Expected people on ART with {sim.label} format, got {n_art}'

    if do_plot:
        fig, ax = pl.subplots()
        ax.bar(results.keys(), results.values())
        ax.set_ylabel('People on ART (final)')
        ax.set_title('ART coverage format comparison')

    return results


@sc.timer()
def test_art_mixed_coverage(do_plot=do_plot):
    """ Check that mixed n/p coverage format works (dual-column DF and dict with format list) """
    sc.heading('Testing mixed n/p coverage format...')

    # Dual-column DataFrame: n_art for early years, p_art for later
    years = np.arange(2000, 2021)
    df = pd.DataFrame(index=years)
    df['n_art'] = np.nan
    df['p_art'] = np.nan
    df.loc[2000:2010, 'n_art'] = np.linspace(0, 5000, 11)
    df.loc[2011:2020, 'p_art'] = np.linspace(0.5, 0.9, 10)

    art1 = sti.ART(coverage=df)
    sim1 = hivsim.demo('simple', run=False, plot=False, n_agents=n_agents, dur=5)
    sim1.pars.interventions = [sti.HIVTest(name='hiv_test', test_prob_data=0.3), art1]
    sim1.label = 'dual_column_df'

    # Dict with per-entry format list
    art2 = sti.ART(coverage={
        'year':   [2000, 2010, 2015, 2020],
        'value':  [0,    5000, 0.5,  0.9],
        'format': ['n',  'n',  'p',  'p'],
    })
    sim2 = hivsim.demo('simple', run=False, plot=False, n_agents=n_agents, dur=5)
    sim2.pars.interventions = [sti.HIVTest(name='hiv_test', test_prob_data=0.3), art2]
    sim2.label = 'mixed_format_dict'

    msim = ss.parallel([sim1, sim2])
    for sim in msim.sims:
        assert sim.results.hiv.n_on_art[-1] > 0, f'Expected people on ART with {sim.label}'

    return msim


@sc.timer()
def test_vmmc_specs(do_plot=do_plot):
    """ Check that scalar and dict VMMC coverage produce circumcisions """
    sc.heading('Testing VMMC coverage formats...')

    specs = dict(
        scalar = sti.VMMC(coverage=0.3),
        dict_p = sti.VMMC(coverage={'year': [2000, 2020], 'value': [0, 0.3]}),
        df_p   = sti.VMMC(coverage=pd.DataFrame(index=np.arange(2000, 2021), data={'p_vmmc': np.linspace(0, 0.3, 21)})),
    )

    sims = []
    for label, vmmc in specs.items():
        sim = hivsim.demo('simple', run=False, plot=False, n_agents=n_agents, dur=5)
        sim.pars.interventions = [sti.HIVTest(name='hiv_test', test_prob_data=0.1), vmmc]
        sim.label = label
        sims.append(sim)

    msim = ss.parallel(sims)

    results = dict()
    for sim in msim.sims:
        n_circ = sim.results.vmmc.n_circumcised[-1]
        results[sim.label] = n_circ
        assert n_circ > 0, f'Expected circumcisions with {sim.label} format, got {n_circ}'

    return results


@sc.timer()
def test_art_stratified_coverage(do_plot=do_plot):
    """ Check that age/sex stratified coverage data is parsed and applied """
    sc.heading('Testing stratified ART coverage...')

    rows = []
    for year in [2005, 2015]:
        for gender in [0, 1]:
            for age_bin in ['[15,25)', '[25,35)', '[35,45)']:
                p = 0.3 if year == 2005 else 0.7
                p *= (1.2 if gender == 0 else 1.0)  # Higher for women (0=female)
                rows.append(dict(Year=year, Gender=gender, AgeBin=age_bin, NationalARTPrevalence=min(p, 1.0)))
    strat_df = pd.DataFrame(rows)

    sim = hivsim.demo('simple', run=False, plot=False, n_agents=n_agents, dur=5)
    sim.pars.interventions = [sti.HIVTest(name='hiv_test', test_prob_data=0.3), sti.ART(coverage=strat_df)]
    sim.run()
    assert sim.results.hiv.n_on_art[-1] > 0, 'Expected people on ART with stratified coverage (canonical names)'

    # Test with lowercase column names + Sex alias + string gender values
    rows_lc = []
    for year in [2005, 2015]:
        for sex in ['f', 'm']:
            for age_bin in ['[15,25)', '[25,35)', '[35,45)']:
                p = 0.4 if year == 2005 else 0.8
                rows_lc.append(dict(year=year, sex=sex, agebin=age_bin, coverage=min(p, 1.0)))

    sim2 = hivsim.demo('simple', run=False, plot=False, n_agents=n_agents, dur=5)
    sim2.pars.interventions = [sti.HIVTest(name='hiv_test', test_prob_data=0.3), sti.ART(coverage=pd.DataFrame(rows_lc))]
    sim2.run()
    assert sim2.results.hiv.n_on_art[-1] > 0, 'Expected people on ART with stratified coverage (string sex values)'

    return sim, sim2


# %% Functional tests
@sc.timer()
def test_art_effects(do_plot=do_plot):
    """ Check that ART puts people on treatment with testing, and warns without """
    sc.heading('Testing ART functional effects...')

    # With testing → people on ART
    sim_with = hivsim.demo('simple', run=False, plot=False, n_agents=n_agents, dur=5)
    sim_with.run()
    n_art = sim_with.results.hiv.n_on_art[-1]
    assert n_art > 0, f'Expected people on ART with the default demo, got {n_art}'

    # Without testing → warning raised and no-one on ART
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        sim_without = sti.Sim(
            diseases=sti.HIV(beta_m2f=0.05, init_prev=0.05),
            networks=sti.StructuredSexual(),
            demographics=[ss.Pregnancy(fertility_rate=10), ss.Deaths(death_rate=10)],
            interventions=[sti.ART(coverage=0.9)],
            n_agents=n_agents, dur=5, verbose=-1,
        )
        sim_without.run()
        warn_texts = [str(x.message) for x in w]
        has_testing_warn = any('testing' in t.lower() for t in warn_texts)

    assert has_testing_warn, 'Expected a warning about missing HIV testing intervention'
    n_art_no_test = sim_without.results.hiv.n_on_art[-1]
    assert n_art_no_test == 0, f'Expected no-one on ART without testing, got {n_art_no_test}'

    if do_plot:
        fig, axes = pl.subplots(1, 2, figsize=(12, 5))
        axes[0].plot(sim_with.t.yearvec, sim_with.results.hiv.n_on_art)
        axes[0].set_title('With testing')
        axes[0].set_ylabel('People on ART')
        axes[1].plot(sim_without.t.yearvec, sim_without.results.hiv.n_on_art)
        axes[1].set_title('Without testing')

    return sim_with, sim_without


@sc.timer()
def test_art_no_pregnancy(do_plot=do_plot):
    """ Check that ART works when pregnancy is not in the sim (issue #313) """
    sc.heading('Testing ART without pregnancy module...')

    sim = sti.Sim(
        diseases=sti.HIV(beta_m2f=0.05, init_prev=0.05),
        networks=sti.StructuredSexual(),
        demographics=[ss.Deaths(death_rate=10)],
        interventions=[
            sti.HIVTest(name='hiv_test', test_prob_data=0.3),
            sti.ART(coverage=0.8),
        ],
        n_agents=n_agents, dur=20, verbose=-1,
    )
    sim.run()
    assert sim.results.hiv.n_on_art[-1] > 0, 'Expected people on ART without pregnancy module'

    return sim


# %% Coverage matching tests
@sc.timer()
def test_art_coverage_matching(do_plot=do_plot):
    """ Check that ART intervention matches supplied coverage data """
    sc.heading('Testing ART coverage matching...')

    # Simple constant proportion coverage
    target_p = 0.7
    sim = hivsim.demo('simple', run=False, plot=False, n_agents=1_000, dur=5)
    sim.pars.interventions = [
        sti.HIVTest(name='hiv_test', test_prob_data=0.5),
        sti.ART(coverage=target_p),
    ]
    sim.pars.analyzers = [sti.art_coverage()]
    sim.run()

    ac = sim.results.art_coverage
    # Check coverage in the last few years (after ramp-up)
    rtol = 0.3  # Generous tolerance for stochastic model with 1K agents
    final_p = np.mean(ac.p_art[-24:])  # Last 2 years of monthly data
    assert np.isclose(final_p, target_p, rtol=rtol), \
        f'Expected ART coverage ~{target_p}, got {final_p:.2f} (rtol={rtol})'

    if do_plot:
        fig, axes = pl.subplots(1, 2, figsize=(12, 5))
        axes[0].plot(sim.t.yearvec, ac.n_art, label='On ART')
        axes[0].set_ylabel('Count')
        axes[0].set_title('ART numbers')
        axes[1].plot(sim.t.yearvec, ac.p_art, label='Coverage')
        axes[1].axhline(target_p, color='k', ls='--', label=f'Target ({target_p})')
        axes[1].set_ylabel('Proportion')
        axes[1].set_title('ART coverage')
        axes[1].legend()

    return sim


@sc.timer()
def test_art_stratified_coverage_matching(do_plot=do_plot):
    """ Check that age-stratified ART coverage approximately matches supplied targets """
    sc.heading('Testing stratified ART coverage matching...')

    # Build stratified coverage data: higher coverage for older age groups and women.
    # 4 strata (2 age bins × 2 sexes) keeps enough infected agents per stratum at
    # n_agents=1k for the per-stratum assertion to be meaningful.
    age_bins = [15, 30, 100]
    targets = {}
    rows = []
    # Use constant coverage across the sim window (start=2000, dur=5) — this avoids
    # ramp-up artefacts and lets the assertion compare the final state directly to
    # the target. Both years carry the same per-stratum target.
    for year in [2000, 2005]:
        for gender in [0, 1]:
            for lo, hi in zip(age_bins[:-1], age_bins[1:]):
                ab = f'[{lo},{hi})'
                age_factor = 1 + (lo - 15) * 0.01  # Slightly higher for older
                sex_factor = 1.1 if gender == 0 else 1.0  # Higher for women (0=female)
                p = min(0.7 * age_factor * sex_factor, 0.95)
                rows.append(dict(Year=year, Gender=gender, AgeBin=ab, coverage=p))
                if year == 2005:
                    targets[(ab, gender)] = p

    strat_df = pd.DataFrame(rows)
    sim = hivsim.demo('simple', run=False, plot=False, n_agents=1_000, dur=5)
    sim.pars.interventions = [
        sti.HIVTest(name='hiv_test', test_prob_data=0.5),
        sti.ART(coverage=strat_df),
    ]
    sim.pars.analyzers = [sti.art_coverage(age_bins=age_bins)]
    sim.run()

    ac = sim.results.art_coverage

    # Check aggregate coverage is nonzero (stratified data also produces a meaningful total)
    final_p = np.mean(ac.p_art[-24:])
    assert final_p > 0.1, f'Expected meaningful aggregate ART coverage with stratified data, got {final_p:.2f}'

    # Per-stratum coverage should approximately match the supplied targets — this is the
    # whole point of stratified input. Tolerance is generous to allow for finite-population
    # noise (1k agents split across 4 strata) and for diagnosis lag (newly-eligible agents
    # must first be diagnosed by HIVTest before ART can initiate them).
    tol = 0.20
    for (ab, gender), target_p in targets.items():
        lo, hi = ss.parse_age_range(ab)
        sex = 'f' if gender == 0 else 'm'
        key = f'p_art_{sex}_{int(lo)}_{int(hi)}'
        measured_p = np.mean(ac[key][-24:])
        print(f'  {key}: target={target_p:.2f}, measured={measured_p:.2f}')
        assert abs(measured_p - target_p) < tol, (
            f'Stratum {key}: measured {measured_p:.2f} vs target {target_p:.2f} '
            f'(tol {tol}). Stratified ART correction is not preserving per-stratum targets.'
        )

    if do_plot:
        fig, axes = pl.subplots(1, 2, figsize=(14, 5))
        for sex, ax in [('f', axes[0]), ('m', axes[1])]:
            for i in range(len(age_bins) - 1):
                lo, hi = age_bins[i], age_bins[i + 1]
                key = f'p_art_{sex}_{lo}_{hi}'
                ax.plot(sim.t.yearvec, ac[key], label=f'{lo}-{hi}')
            ax.set_ylabel('ART coverage')
            ax.set_title(f'{"Women" if sex == "f" else "Men"}')
            ax.legend()
            ax.set_ylim(0, 1)

    return sim


# %% Scientific validation skeletons
@sc.timer()
def test_art_reduces_mortality(do_plot=do_plot):
    """ Check that HIV mortality is lower with ART than without """
    sc.heading('Testing ART reduces mortality...')

    # Without ART
    sim_no_art = hivsim.demo('simple', run=False, plot=False, n_agents=1_000, dur=5)
    sim_no_art.pars.interventions = [sti.HIVTest(name='hiv_test', test_prob_data=0.3)]
    sim_no_art.label = 'no_art'

    # With ART
    sim_art = hivsim.demo('simple', run=False, plot=False, n_agents=1_000, dur=5)
    sim_art.pars.interventions = [
        sti.HIVTest(name='hiv_test', test_prob_data=0.3),
        sti.ART(coverage=0.8),
    ]
    sim_art.label = 'with_art'

    msim = ss.parallel([sim_no_art, sim_art])
    sim_no_art, sim_art = msim.sims

    deaths_no_art = sim_no_art.results.hiv.new_deaths.sum()
    deaths_art    = sim_art.results.hiv.new_deaths.sum()

    if do_plot:
        fig, ax = pl.subplots()
        ax.plot(sim_no_art.t.yearvec, np.cumsum(sim_no_art.results.hiv.new_deaths), label='No ART')
        ax.plot(sim_art.t.yearvec, np.cumsum(sim_art.results.hiv.new_deaths), label='With ART')
        ax.legend()
        ax.set_ylabel('Cumulative HIV deaths')
        ax.set_title('ART mortality reduction')

    # TODO: uncomment once parameterized and validated with sufficient population
    # assert deaths_art < deaths_no_art, f'Expected fewer deaths with ART ({deaths_art}) than without ({deaths_no_art})'

    return sim_no_art, sim_art


@sc.timer()
def test_art_parameter_sensitivity(do_plot=do_plot):
    """ Check that ART coverage and initiation parameters affect outcomes as expected """
    sc.heading('Testing ART parameter sensitivity...')

    par_effects = dict(
        coverage       = [0.2, 0.9],
        art_initiation = [0.3, 0.95],
    )

    results = dict()
    sims = []
    for par, (lo, hi) in par_effects.items():
        for val in [lo, hi]:
            if par == 'coverage':
                art = sti.ART(coverage=val)
            else:
                # No coverage cap — art_initiation is the only constraint, giving
                # a clear signal: higher p → faster uptake → more on ART at end.
                art = sti.ART(art_initiation=ss.bernoulli(p=val))
            sim = hivsim.demo('simple', run=False, plot=False, n_agents=n_agents, dur=10, rand_seed=0)
            sim.pars.interventions = [sti.HIVTest(name='hiv_test', test_prob_data=0.3), art]
            sim.label = (par, val)
            sims.append(sim)
    
    msim = ss.parallel(sims)

    for sim in msim.sims:
        n_on_art = sim.results.hiv.n_on_art[-1]
        results[sim.label] = n_on_art

    for par, (lo, hi) in par_effects.items():
        v_lo = results[(par, lo)]
        v_hi = results[(par, hi)]
        assert v_lo <= v_hi, f'Expected higher {par} to increase ART uptake, but {par}={lo} gave {v_lo} vs {par}={hi} gave {v_hi}'

    return results


@sc.timer()
def test_art_duration(do_plot=do_plot):
    """ Check that agents remain on ART for expected durations """
    sc.heading('Testing ART duration tracking...')

    # Set a known constant ART duration and create a mini sim
    dur_years = 5
    run_years = dur_years - 1
    sim_kwargs = dict(n_agents=10, dur=run_years, start=1990, verbose=-1)

    # HIV parameters: everyone starts infected and with a fixed duration on ART
    hiv_pars = dict(init_prev=1.0, dur_on_art=ss.constant(v=ss.years(dur_years)))

    # Interventions: 100% testing and ART initiation to ensure everyone starts ART at ti~0
    interventions=[
        sti.HIVTest(test_prob_data=1.0, dt_scale=False),
        sti.ART(art_initiation=1.0),
    ]

    # Make the sim
    sim = hivsim.demo(
        'simple', plot=False,
        demographics=[ss.Deaths(death_rate=0)],  # No deaths to isolate ART duration effect
        interventions=interventions,
        **hiv_pars,
        **sim_kwargs,
    )

    # All agents should still be on ART — the sim is shorter than dur_on_art
    hiv = sim.diseases.hiv
    n_alive  = len(sim.people.alive.uids)
    n_on_art = len(hiv.on_art.uids)
    assert n_on_art == n_alive, f'Expected all {n_alive} living agents on ART, got {n_on_art}'

    # Everyone started ART at ti~0, so mean duration should approximate run_years
    on_art = hiv.on_art.uids
    ti_start = hiv.ti_art[on_art]
    durations_years = (sim.t.ti - ti_start) * float(sim.pars.dt)
    mean_dur = float(np.nanmean(durations_years))
    assert np.isclose(mean_dur, run_years, atol=0.5), f'Expected mean ART duration ~{run_years} years, got {mean_dur:.1f}'

    if do_plot:
        fig, ax = pl.subplots()
        ax.hist(durations_years[~np.isnan(durations_years)], bins=20)
        ax.set_xlabel('Years on ART')
        ax.set_ylabel('Count')
        ax.set_title(f'ART duration distribution (mean={mean_dur:.1f} years)')

    return sim


@sc.timer()
def test_art_dx2tx_delay(do_plot=do_plot):
    """ Check dur_dx2tx on HIVTest enforces a delay between diagnosis and ART start """
    sc.heading('Testing dur_dx2tx (diagnosis → treatment) delay...')

    # Build two sims that only differ in the dx→tx delay
    delay_years = 2
    sim_kwargs = dict(n_agents=200, dur=delay_years + 2, start=1990, verbose=-1)
    hiv_pars = dict(init_prev=1.0, dur_on_art=ss.constant(v=ss.years(10)))

    def make_sim(delay):
        return hivsim.demo(
            'simple', plot=False,
            demographics=[ss.Deaths(death_rate=0)],
            interventions=[
                sti.HIVTest(test_prob_data=1.0, dt_scale=False,
                            dur_dx2tx=ss.constant(ss.years(delay))),
                sti.ART(art_initiation=1.0),
            ],
            **hiv_pars, **sim_kwargs,
        )

    # With delay: agents diagnosed at ti=0 should start ART at ti = delay_ti
    sim = make_sim(delay_years)
    hiv = sim.diseases.hiv
    dx_uids = hiv.diagnosed.uids
    n_alive = len(sim.people.alive.uids)
    assert len(dx_uids) == n_alive, 'Expected all living agents diagnosed'

    delay_ti = int(round(delay_years / float(sim.pars.dt)))
    expected = hiv.ti_diagnosed[dx_uids] + delay_ti
    on_art_uids = hiv.on_art.uids
    assert np.allclose(hiv.ti_art[on_art_uids], expected[np.isin(dx_uids, on_art_uids)]), \
        f'Expected ti_art == ti_diagnosed + {delay_ti} for on-ART agents'

    # No-one should be on ART before the delay elapses
    # The sim runs (delay_years + 2) years; check n_on_art ramp
    n_on_art = sim.results.hiv.n_on_art
    early_ti = max(0, delay_ti - 1)
    assert n_on_art[early_ti] == 0, f'No agents should be on ART before ti={delay_ti}, got {n_on_art[early_ti]}'
    assert n_on_art[delay_ti] > 0, f'Agents should be on ART by ti={delay_ti}'

    # Default (no delay): agents diagnosed at ti=0 should be on ART at ti=0
    sim2 = make_sim(0)
    hiv2 = sim2.diseases.hiv
    dx2 = hiv2.diagnosed.uids
    assert np.allclose(hiv2.ti_art[dx2], hiv2.ti_diagnosed[dx2]), \
        'Default dur_dx2tx=0 should leave ti_art == ti_diagnosed'
    assert sim2.results.hiv.n_on_art[0] > 0, 'With no delay, agents should be on ART from ti=0'

    return sim, sim2


@sc.timer()
def test_pn():
    """ PartnerNotification: with PN, more HIV diagnoses occur than without """
    sc.heading('Testing PartnerNotification...')

    def newly_diagnosed(sim):
        hiv = sim.diseases.hiv
        return (hiv.ti_diagnosed == hiv.ti).uids

    def make_sim(use_pn, rand_seed=1):
        hiv_test = sti.HIVTest(name='hiv_test', test_prob_data=0.3)
        intvs = [hiv_test]
        if use_pn:
            intvs.append(sti.PartnerNotification(
                eligibility=newly_diagnosed,
                test=hiv_test,
                p_notify_current=ss.bernoulli(p=0.7),
                p_attends_current=ss.bernoulli(p=0.7),
                p_notify_previous=ss.bernoulli(p=0.3),
                p_attends_previous=ss.bernoulli(p=0.5),
                name='pn',
            ))
        return hivsim.demo(
            'simple', run=False, plot=False,
            n_agents=n_agents, dur=5, rand_seed=rand_seed, verbose=-1,
            networks=[sti.StructuredSexual(recall_prior=True),
                      sti.PriorPartners(dur_recall=ss.years(0.5))],
            interventions=intvs,
        )

    no_pn   = make_sim(use_pn=False); no_pn.run()
    with_pn = make_sim(use_pn=True);  with_pn.run()

    n_dx_no_pn   = int(no_pn.results.hiv_test.new_diagnoses.sum())
    n_dx_with_pn = int(with_pn.results.hiv_test.new_diagnoses.sum())
    assert n_dx_with_pn > n_dx_no_pn, \
        f'Expected more diagnoses with PN ({n_dx_with_pn}) than without ({n_dx_no_pn})'

    pn_res = with_pn.results.pn
    assert int(pn_res.new_attending.sum()) <= int(pn_res.new_notified.sum()), \
        'Attendance should not exceed notification'
    return n_dx_no_pn, n_dx_with_pn


def test_pn_rates():
    """ pn_rates: edge-type and partner-sex stratification of the probability callable """
    sc.heading('Testing pn_rates stratification...')

    def evaluate(rates, partner_edges, edge_types, uids, female=None):
        """ Run a pn_rates callable against stubbed module/sim state. """
        nw = sc.dictobj(edge_types=edge_types)
        module = sc.dictobj(current_partner_edges=partner_edges, _cur_nw=nw)
        sim = sc.dictobj(people=sc.dictobj(female=female))
        return sti.pn_rates(rates)(module, sim, ss.uids(uids))

    # Single edge type: matched partner gets its rate, non-partner gets 0
    one = {'stable': 0}
    assert np.allclose(evaluate({'stable': 0.2}, {1: [0]}, one, [1, 2]), [0.2, 0.0])

    # Two edge types, incl. a partner reachable via both edges (rates summed)
    two = {'stable': 0, 'casual': 1}
    assert np.allclose(evaluate({'stable': 0.2, 'casual': 0.1}, {1: [0], 2: [1], 3: [0, 1]}, two, [1, 2, 3]),
                       [0.2, 0.1, 0.3])

    # Multi-edge probabilities sum but cap at 1
    assert np.allclose(evaluate({'stable': 0.8}, {1: [0, 0]}, one, [1]), [1.0])

    # Edge type absent from `rates` contributes 0 (the "pn_rates 0" case)
    assert np.allclose(evaluate({'stable': 0.5}, {1: [2]}, {'stable': 0, 'onetime': 2}, [1]), [0.0])

    # Empty current_partner_edges (e.g. the prior channel) -> all zeros
    assert np.allclose(evaluate({'stable': 0.5}, {}, one, [1, 2]), [0.0, 0.0])

    # Per-sex rates resolve by partner sex (uid 1 female, uid 2 male)
    female = np.array([False, True, False])
    assert np.allclose(evaluate({'stable': {'f': 0.8, 'm': 0.5}}, {1: [0], 2: [0]}, one, [1, 2], female),
                       [0.8, 0.5])

    # Mixing scalar and per-sex specs raises
    try:
        sti.pn_rates({'stable': {'f': 0.8, 'm': 0.5}, 'casual': 0.1})
        raise AssertionError('Expected ValueError for mixed scalar/dict spec')
    except ValueError:
        pass
    return


def test_pn_rates_sim():
    """ pn_rates through PartnerNotification.step() across one and two networks """
    sc.heading('Testing pn_rates in a running sim...')

    def newly_diagnosed(sim):
        hiv = sim.diseases.hiv
        return (hiv.ti_diagnosed == hiv.ti).uids

    def make_sim(current_rates, use_prior=False):
        # Notify/attend probs are deterministic (matched edge -> 1, else 0),
        # so current notifications reflect exactly which partners matched.
        hiv_test = sti.HIVTest(name='hiv_test', test_prob_data=0.3)
        pn = sti.PartnerNotification(
            eligibility=newly_diagnosed, test=hiv_test,
            p_notify_current=ss.bernoulli(p=sti.pn_rates(current_rates)),
            p_attends_current=ss.bernoulli(p=1.0),
            p_notify_previous=ss.bernoulli(p=0.5 if use_prior else 0.0),
            p_attends_previous=ss.bernoulli(p=1.0),
            name='pn',
        )
        networks = [sti.StructuredSexual(recall_prior=use_prior)]
        if use_prior:
            networks.append(sti.PriorPartners(dur_recall=ss.years(0.5)))
        return hivsim.demo('simple', run=False, plot=False, n_agents=n_agents,
                           dur=5, rand_seed=1, verbose=-1,
                           networks=networks, interventions=[hiv_test, pn])

    edge_types = sti.StructuredSexual().edge_types

    # Single network: matching rates notify partners; an unmatched edge type notifies none
    match = make_sim({k: 1.0 for k in edge_types}); match.run()
    none  = make_sim({'no_such_edge_type': 1.0});   none.run()
    n_match = int(match.results.pn.new_notified_current.sum())
    assert n_match > 0, f'Expected current notifications with matching rates, got {n_match}'
    assert int(none.results.pn.new_notified_current.sum()) == 0, 'Unmatched edge type should notify nobody'

    # Two networks: with a prior network also configured, the current channel
    # still notifies via pn_rates. The prior channel is current-channel-blind:
    # pn_rates returns 0 there (covered deterministically in test_pn_rates).
    both = make_sim({k: 1.0 for k in edge_types}, use_prior=True); both.run()
    assert int(both.results.pn.new_notified_current.sum()) > 0, 'Expected current-channel notifications'
    return n_match


if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    T = sc.timer()

    r1  = test_art_specs(do_plot=do_plot)
    r2  = test_art_mixed_coverage(do_plot=do_plot)
    r3  = test_vmmc_specs(do_plot=do_plot)
    r4  = test_art_stratified_coverage(do_plot=do_plot)
    r5  = test_art_effects(do_plot=do_plot)
    r6  = test_art_no_pregnancy(do_plot=do_plot)
    r7  = test_art_coverage_matching(do_plot=do_plot)
    r8  = test_art_stratified_coverage_matching(do_plot=do_plot)
    r9  = test_art_reduces_mortality(do_plot=do_plot)
    r10 = test_art_parameter_sensitivity(do_plot=do_plot)
    r11 = test_art_duration(do_plot=do_plot)
    r12 = test_pn()
    r13 = test_pn_rates()
    r14 = test_pn_rates_sim()

    T.toc()

    if do_plot:
        pl.show()
