"""
Benchmark and diagnostics for seven PFA variants.

Runs each variant at n=1k, 5k, 10k (LSA capped at 5k) for several reps,
then again on hivsim.demo('zimbabwe'). Saves results to
``pfa_comparison_results.obj`` for the companion notebook.

Usage::

    cd stisim/tests/devtests
    python devtest_pfa_comparison.py            # full run
    python devtest_pfa_comparison.py --quick    # n=1k, 1 rep, skip Zimbabwe
"""
import argparse
import sys
from pathlib import Path

import sciris as sc
import starsim as ss
import stisim as sti

sys.path.insert(0, str(Path(__file__).parent))
from pfa_diagnostics import (  # noqa: E402
    PartnersLastYearAnalyzer,
    PairFormationAgesAnalyzer,
    PairPrevalenceAnalyzer,
)


# (display name, match_method string key into matchers.MATCHERS)
VARIANTS = [
    ('SortBisect',       'sort_bisect'),
    ('SortPair',         'sort_pair'),
    ('LSA',              'lsa'),
    ('DesiredAgeBucket', 'desired_age_bucket'),
    ('GreedyOldEnough',  'greedy_old_enough'),
    ('KDTreeNN',         'kdtree_nn'),
    ('BandMatch',        'band_match'),
    ('ClosestAgeTaperedSeeking', 'closest_age_tapered_seeking')
]

LSA_MAX_N = 5_000


def run_generic(n_agents_list, n_reps, sim_years):
    """Generic benchmark on ss.Sim + MFNetwork variant + HIV."""
    results = sc.objdict()
    for n in n_agents_list:
        for name, method in VARIANTS:
            if name == 'LSA' and n > LSA_MAX_N:
                print(f'  skipping LSA at n={n} (capped at {LSA_MAX_N})')
                continue
            for rep in range(n_reps):
                key = (name, n, rep)
                print(f'Running {key} ...')
                sim = ss.Sim(
                    n_agents=n,
                    start='2000-01-01',
                    stop=f'{2000+sim_years}-01-01',
                    networks=sti.MFNetwork(match_method=method),
                    diseases=sti.HIV(),
                    analyzers=[
                        PartnersLastYearAnalyzer(),
                        PairFormationAgesAnalyzer(),
                        PairPrevalenceAnalyzer(),
                    ],
                    rand_seed=rep,
                    verbose=0,
                )
                with sc.timer() as t:
                    sim.run()
                mf_net = next(nw for nw in sim.networks() if isinstance(nw, sti.MFNetwork))
                results[key] = sc.objdict(
                    wall_time=t.elapsed,
                    partners_last_year=sim.analyzers.partnerslastyearanalyzer.records,
                    pair_formation_ages=sim.analyzers.pairformationagesanalyzer.records,
                    pair_prevalence=sim.analyzers.pairprevalenceanalyzer.records,
                    lifetime_partners=mf_net.lifetime_partners.values.copy(),
                    sex_male=sim.people.male.values.copy(),
                    age=sim.people.age.values.copy(),
                )
                print(f'    wall_time={t.elapsed:.2f}s')
    return results


def run_zimbabwe(n_reps, n_agents=10_000, checkpoint=None):
    """Calibrated benchmark on hivsim.demo('zimbabwe')."""
    import hivsim  # noqa: F401  (used to register hivsim_examples.zimbabwe with sc.importbyname)
    results = sc.objdict()
    for name, method in VARIANTS:
        if name == 'LSA':
            continue  # too slow for full hivsim runs
        for rep in range(n_reps):
            key = (name, 'zimbabwe', rep)
            print(f'Running {key} ...')
            sim = hivsim.demo('zimbabwe', run=False, plot=False, n_agents=n_agents)
            orig_net = None
            for net in sim.pars['networks']:
                if isinstance(net, sti.StructuredSexual):
                    orig_net = net
                    break
            if orig_net is None:
                raise RuntimeError('No StructuredSexual found in hivsim.demo(zimbabwe) networks')
            orig_net.pars.match_method = method
            sim.pars['rand_seed'] = rep
            sim.pars['analyzers'] = list(sim.pars.get('analyzers') or []) + [
                PartnersLastYearAnalyzer(),
                PairFormationAgesAnalyzer(),
                PairPrevalenceAnalyzer(),
            ]
            with sc.timer() as t:
                sim.run()
            mf_net = next(nw for nw in sim.networks() if isinstance(nw, sti.MFNetwork))
            results[key] = sc.objdict(
                wall_time=t.elapsed,
                hiv_prevalence=sim.results.hiv.prevalence.values.copy(),
                partners_last_year=sim.analyzers.partnerslastyearanalyzer.records,
                pair_formation_ages=sim.analyzers.pairformationagesanalyzer.records,
                pair_prevalence=sim.analyzers.pairprevalenceanalyzer.records,
                lifetime_partners=mf_net.lifetime_partners.values.copy(),
                sex_male=sim.people.male.values.copy(),
                age=sim.people.age.values.copy(),
            )
            print(f'    wall_time={t.elapsed:.2f}s')
            if checkpoint is not None:
                checkpoint(results)
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--quick', action='store_true', help='1k agents, 1 rep, skip Zimbabwe')
    parser.add_argument('--skip-zimbabwe', action='store_true')
    args = parser.parse_args()

    if args.quick:
        n_agents_list = [1_000]
        n_reps = 1
        sim_years = 5
    else:
        n_agents_list = [1_000, 5_000, 10_000]
        n_reps = 3
        sim_years = 20

    out = Path(__file__).parent / 'pfa_comparison_results.obj'
    all_results = sc.objdict()

    def save_zimbabwe(partial):
        all_results.zimbabwe = partial
        sc.save(out, all_results)

    all_results.generic = run_generic(n_agents_list, n_reps, sim_years)
    sc.save(out, all_results)  # save generic results before attempting Zimbabwe
    if not args.quick and not args.skip_zimbabwe:
        all_results.zimbabwe = run_zimbabwe(n_reps, checkpoint=save_zimbabwe)
        sc.save(out, all_results)
    print(f'Saved {out}')


if __name__ == '__main__':
    main()
