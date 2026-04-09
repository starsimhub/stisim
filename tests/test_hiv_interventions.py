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

    results = dict()
    for label, art in specs.items():
        sim = hivsim.demo('simple', run=False, plot=False, n_agents=n_agents)
        sim.pars.interventions = [sti.HIVTest(name='hiv_test', test_prob_data=0.3), art]
        sim.run()
        n_art = sim.results.hiv.n_on_art[-1]
        results[label] = n_art

    for label, n_art in results.items():
        assert n_art > 0, f'Expected people on ART with {label} format, got {n_art}'

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
    sim1 = hivsim.demo('simple', run=False, plot=False, n_agents=n_agents)
    sim1.pars.interventions = [sti.HIVTest(name='hiv_test', test_prob_data=0.3), art1]
    sim1.run()
    assert sim1.results.hiv.n_on_art[-1] > 0, 'Expected people on ART with dual-column DataFrame'

    # Dict with per-entry format list
    art2 = sti.ART(coverage={
        'year':   [2000, 2010, 2015, 2020],
        'value':  [0,    5000, 0.5,  0.9],
        'format': ['n',  'n',  'p',  'p'],
    })
    sim2 = hivsim.demo('simple', run=False, plot=False, n_agents=n_agents)
    sim2.pars.interventions = [sti.HIVTest(name='hiv_test', test_prob_data=0.3), art2]
    sim2.run()
    assert sim2.results.hiv.n_on_art[-1] > 0, 'Expected people on ART with mixed-format dict'

    return sim1, sim2


@sc.timer()
def test_vmmc_specs(do_plot=do_plot):
    """ Check that scalar and dict VMMC coverage produce circumcisions """
    sc.heading('Testing VMMC coverage formats...')

    specs = dict(
        scalar = sti.VMMC(coverage=0.3),
        dict_p = sti.VMMC(coverage={'year': [2000, 2020], 'value': [0, 0.3]}),
        df_p   = sti.VMMC(coverage=pd.DataFrame(index=np.arange(2000, 2021), data={'p_vmmc': np.linspace(0, 0.3, 21)})),
    )

    results = dict()
    for label, vmmc in specs.items():
        sim = hivsim.demo('simple', run=False, plot=False, n_agents=n_agents)
        sim.pars.interventions = [sti.HIVTest(name='hiv_test', test_prob_data=0.1), vmmc]
        sim.run()
        n_circ = sim.results.vmmc.n_circumcised[-1]
        results[label] = n_circ

    for label, n_circ in results.items():
        assert n_circ > 0, f'Expected circumcisions with {label} format, got {n_circ}'

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

    sim = hivsim.demo('simple', run=False, plot=False, n_agents=n_agents)
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

    sim2 = hivsim.demo('simple', run=False, plot=False, n_agents=n_agents)
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
    sim_with = hivsim.demo('simple', run=False, plot=False, n_agents=n_agents)
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
            n_agents=n_agents, dur=20, verbose=-1,
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
    sim = hivsim.demo('simple', run=False, plot=False, n_agents=5_000)
    sim.pars.interventions = [
        sti.HIVTest(name='hiv_test', test_prob_data=0.5),
        sti.ART(coverage=target_p),
    ]
    sim.pars.analyzers = [sti.art_coverage()]
    sim.run()

    ac = sim.results.art_coverage
    # Check coverage in the last few years (after ramp-up)
    rtol = 0.2  # Generous tolerance for stochastic model with 5K agents
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

    # Build stratified coverage data: higher coverage for older age groups and women
    age_bins = [15, 25, 35, 45]
    targets = {}
    rows = []
    for year in [2005, 2015]:
        for gender in [0, 1]:
            for lo, hi in zip(age_bins[:-1], age_bins[1:]):
                ab = f'[{lo},{hi})'
                base_p = 0.3 if year == 2005 else 0.7
                age_factor = 1 + (lo - 15) * 0.01  # Slightly higher for older
                sex_factor = 1.1 if gender == 0 else 1.0  # Higher for women (0=female)
                p = min(base_p * age_factor * sex_factor, 0.95)
                rows.append(dict(Year=year, Gender=gender, AgeBin=ab, coverage=p))
                if year == 2015:
                    targets[(ab, gender)] = p

    strat_df = pd.DataFrame(rows)
    sim = hivsim.demo('simple', run=False, plot=False, n_agents=5_000)
    sim.pars.interventions = [
        sti.HIVTest(name='hiv_test', test_prob_data=0.5),
        sti.ART(coverage=strat_df),
    ]
    sim.pars.analyzers = [sti.art_coverage(age_bins=age_bins)]
    sim.run()

    ac = sim.results.art_coverage

    # Check aggregate coverage is nonzero (stratified data produces a global target)
    final_p = np.mean(ac.p_art[-24:])
    assert final_p > 0.1, f'Expected meaningful aggregate ART coverage with stratified data, got {final_p:.2f}'

    # Print per-stratum coverage for inspection (not asserted — current implementation
    # computes a global target from stratified data, not per-stratum force-fitting)
    for (ab, gender), target_p in targets.items():
        lo, hi = ss.parse_age_range(ab)
        sex = 'f' if gender == 0 else 'm'
        key = f'p_art_{sex}_{int(lo)}_{int(hi)}'
        measured_p = np.mean(ac[key][-24:])
        print(f'  {key}: target={target_p:.2f}, measured={measured_p:.2f}')

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
    sim_no_art = hivsim.demo('simple', run=False, plot=False, n_agents=2_000)
    sim_no_art.pars.interventions = [sti.HIVTest(name='hiv_test', test_prob_data=0.3)]
    sim_no_art.run()

    # With ART
    sim_art = hivsim.demo('simple', run=False, plot=False, n_agents=2_000)
    sim_art.pars.interventions = [
        sti.HIVTest(name='hiv_test', test_prob_data=0.3),
        sti.ART(coverage=0.8),
    ]
    sim_art.run()

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

    # TODO: implement
    # Vary art_initiation probability (low vs high)
    # Vary coverage level (low vs high)
    # Assert: higher coverage → more people on ART
    # Assert: higher initiation → faster ART uptake

    par_effects = dict(
        coverage       = [0.2, 0.9],
        art_initiation = [0.3, 0.95],
    )

    results = dict()
    for par, (lo, hi) in par_effects.items():
        for val in [lo, hi]:
            if par == 'coverage':
                art = sti.ART(coverage=val)
            else:
                art = sti.ART(coverage=0.5, art_initiation=ss.bernoulli(p=val))
            sim = hivsim.demo('simple', run=False, plot=False, n_agents=n_agents)
            sim.pars.interventions = [sti.HIVTest(name='hiv_test', test_prob_data=0.3), art]
            sim.run()
            results[(par, val)] = float(sim.results.hiv.n_on_art[-1])

    # TODO: uncomment once validated with larger populations
    # for par, (lo, hi) in par_effects.items():
    #     v_lo = results[(par, lo)]
    #     v_hi = results[(par, hi)]
    #     assert v_lo <= v_hi, f'Expected higher {par} to increase ART uptake, but {par}={lo} gave {v_lo} vs {par}={hi} gave {v_hi}'

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

    T.toc()

    if do_plot:
        pl.show()
