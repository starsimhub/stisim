"""
HIV natural history verification tests

Tests to ensure appropriate behavior of HIV as a disease absent any treatment.
"""
from itertools import chain

import matplotlib.pyplot as plt
import pytest
import sciris as sc
import sys

from pathlib import Path

import starsim as ss
import stisim as sti

from stisim.interventions.prep import SuppliedPrep
from stisim.interventions.product import Product
from stisim.interventions.supplies import Supplies
from stisim.interventions.supply import Supply

tests_directory = Path(__file__).resolve().parent
sys.path.append(str(tests_directory))

from hiv_natural_history_analyzers import PrepCoverageAnalyzer, PrepEfficacyAnalyzer, PrepDurationAnalyzer, \
    PrepCountsAnalyzer

from testlib import build_testing_sim


verbose = False
do_plot = False
sc.options(interactive=False)

"""
# These are done
default_eligibility
single_eligibility_group
multi_eligibility_groups
each_eligibility_group_has_coverage
single_group_coverage_targets_reached_but_not_exceeded_with_infinite_supply
multi_group_coverage_targets_reached_but_not_exceeded_with_infinite_supply
prep_eff_updates_correctly
prep_under_supply_limits

# Next
prep_eff_updates_correctly_with_reuptake  # later, not "promoting" reuptake for now
prep_eff_is_0_after_drop

PrepManager:
dropout_updates_agents_properly
update_eff_updates_agents_properly
uptake_updates_agents_and_product_properly

"""
pill = Product(name='test-pill', type='prep', delivery_mode='pill', cost=2, eff_by_ti=[0.95])
shot_1y_perfect = Product(name='test-shot', type='prep', delivery_mode='shot', cost=10, eff_by_ti=[1.0] * 12)
shot_1y_imperfect = Product(name='test-shot', type='prep', delivery_mode='shot', cost=10, eff_by_ti=[0.95, 0.9, 0.85, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3, 0.2])

def infinite_prep(product):
    return limited_prep(product=product, quantity=1000000)

def limited_prep(product, quantity):
    return Supplies(supplies=[Supply(product=product, quantity=quantity)])

@sc.timer()
def test_default_eligibility_and_coverage():
    duration = 1 # years
    supplies = infinite_prep(product=shot_1y_perfect)
    intervention = SuppliedPrep(name='this-is-a-test', supplies=supplies)
    target_coverage = intervention.coverages[0]  # should be 100%
    eligibilities = {'fsw': intervention.default_eligibilities[0]}  # there is only one
    analyzer = PrepCoverageAnalyzer(eligibilities=eligibilities, consider_new_infections=True)
    sim = build_testing_sim(analyzers=[analyzer], interventions=[intervention], n_agents=5000, duration=duration) # , death=None, diseases=diseases)
    sim.run()
    analyzer_results = sim.results[analyzer.name][analyzer.result_name]

    # Trimming off final reporting ti as starsim always (currently) runs +1 ti more than requested.
    analyzer_results = {group_name: coverage[:-1] for group_name, coverage in analyzer_results.items()}
    for group_name, coverage in analyzer_results.items():
        assert len(coverage) == (duration * 12)

    if verbose:
        for group_name, coverage in analyzer_results.items():
            print(f"group: {group_name} coverage: {coverage}")

    for group_name, coverage in analyzer_results.items():
        # There should be agents on prep with at all timesteps (given sim/product timescale)
        assert min(coverage) > 0, f"group: {group_name} No agents found on PrEP in 1+ timesteps. Agents should be on PrEP for all timesteps."

        # coverage should generally increase and approach the target (default: 100%). However, due to eligibility counts
        # changing over time, there MAY be temporary, small drops as new FSWs enter the scene.
        delta_coverage = [coverage[i+1] - coverage[i] for i in range(len(coverage)-1)]
        max_decrease = 0.05
        assert min(delta_coverage) >= -1*max_decrease, f"group: {group_name} Coverage decreased more than {max_decrease}% between timesteps at least once."
        # Final set of coverages should be really close (within 5%) to the target of 1 (100%)
        assert abs(target_coverage - min(coverage[-3:])) < 0.05, f"group {group_name} At least one late-sim coverage too far off of target coverage: {target_coverage}."
    return sim

@sc.timer()
def test_default_eligibility_and_modified_coverage():
    """HIV infection disabled to prevent over-covering the default population as the non-users contract HIV (denominator drops -> coverage goes up)"""
    duration = 1 # years
    supplies = infinite_prep(product=shot_1y_perfect)
    target_coverage = 0.8
    intervention = SuppliedPrep(name='this-is-a-test', supplies=supplies, coverages=[target_coverage])
    eligibilities = {'fsw': intervention.default_eligibilities[0]}  # there is only one
    analyzer = PrepCoverageAnalyzer(eligibilities=eligibilities, consider_new_infections=True)
    diseases = [sti.HIV(beta_m2f=0.05, beta_m2c=0.1, init_prev=0.0)]
    sim = build_testing_sim(analyzers=[analyzer], interventions=[intervention], n_agents=5000, duration=duration, diseases=diseases)
    sim.run()
    analyzer_results = sim.results[analyzer.name][analyzer.result_name]

    # Trimming off final reporting ti as starsim always (currently) runs +1 ti more than requested.
    analyzer_results = {group_name: coverage[:-1] for group_name, coverage in analyzer_results.items()}
    for group_name, coverage in analyzer_results.items():
        assert len(coverage) == (duration * 12)

    if verbose:
        for group_name, coverage in analyzer_results.items():
            print(f"group: {group_name} coverage: {coverage}")

    for group_name, coverage in analyzer_results.items():
        # There should be agents on prep with at all timesteps (given sim/product timescale)
        assert min(coverage) > 0, f"group: {group_name} No agents found on PrEP in 1+ timesteps. Agents should be on PrEP for all timesteps."

        # coverage should generally increase and approach the target (80%). However, due to eligibility counts
        # changing over time, there MAY be temporary, small drops as new FSWs enter/leave the scene or die.
        delta_coverage = [coverage[i+1] - coverage[i] for i in range(len(coverage)-1)]
        max_decrease = 0.05
        assert min(delta_coverage) >= -1*max_decrease, f"group: {group_name} Coverage decreased more than {max_decrease}% between timesteps at least once."
        # Final set of coverages should be really close (within 5%) to the target
        assert abs(target_coverage - min(coverage[-3:])) < 0.05, f"group {group_name} At least one late-sim coverage too far off of target coverage: {target_coverage}."
    return sim


@sc.timer()
def test_custom_eligibilities_and_coverages():
    """HIV infection disabled to prevent over-covering the default population as the non-users contract HIV (denominator drops -> coverage goes up)"""
    duration = 3 # years

    # create eligibility & coverage logic
    elig1 = SuppliedPrep.default_eligibilities[0]
    elig2 = lambda sim: (sim.people.female == False) & (~sim.diseases.hiv.infected)
    group_names = ['fsw', 'men']
    elig_funcs = [elig1, elig2]
    target_coverages = [0.8, 0.5]
    eligibilities = dict(zip(group_names, elig_funcs))
    coverage_by_group = dict(zip(group_names, target_coverages))

    # create the intervention & supplies
    supplies = infinite_prep(product=shot_1y_perfect)
    intervention = SuppliedPrep(name='this-is-a-test', supplies=supplies, eligibilities=elig_funcs, coverages=target_coverages)

    diseases = [sti.HIV(beta_m2f=0.05, beta_m2c=0.1, init_prev=0.0)]
    analyzer = PrepCoverageAnalyzer(eligibilities=eligibilities, consider_new_infections=True)
    sim = build_testing_sim(analyzers=[analyzer], interventions=[intervention], n_agents=5000, duration=duration, diseases=diseases)

    sim.run()
    analyzer_results = sim.results[analyzer.name][analyzer.result_name]

    # Trimming off final reporting ti as starsim always (currently) runs +1 ti more than requested.
    analyzer_results = {group_name: coverage[:-1] for group_name, coverage in analyzer_results.items()}
    for group_name, coverage in analyzer_results.items():
        assert len(coverage) == (duration * 12)

    if verbose:
        for group_name, coverage in analyzer_results.items():
            print(f"group: {group_name} coverage: {coverage}")

    for group_name, coverage in analyzer_results.items():
        # There should be agents on prep with at all timesteps (given sim/product timescale)
        assert min(coverage) > 0, f"group: {group_name} No agents found on PrEP in 1+ timesteps. Agents should be on PrEP for all timesteps."

        # coverage should generally increase and approach the target (80%). However, due to eligibility counts
        # changing over time, there MAY be temporary, small drops as new FSWs enter/leave the scene or die.
        delta_coverage = [coverage[i+1] - coverage[i] for i in range(len(coverage)-1)]
        max_decrease = 0.05
        assert min(delta_coverage) >= -1*max_decrease, f"group: {group_name} Coverage decreased more than {max_decrease}% between timesteps at least once."
        # Final set of coverages should be really close (within 5%) to the target
        target = coverage_by_group[group_name]
        assert abs(target - min(coverage[-3:])) < 0.05, f"group {group_name} At least one late-sim coverage too far off of target coverage: {target}."
    return sim

@sc.timer()
def test_prep_eff_updates_correctly():
    """HIV infection disabled to prevent over-covering the default population as the non-users contract HIV (denominator drops -> coverage goes up)"""
    duration = 2 # years

    # create eligibility & coverage logic
    elig1 = SuppliedPrep.default_eligibilities[0]
    elig2 = lambda sim: (sim.people.female == False) & (~sim.diseases.hiv.infected)
    group_names = ['fsw', 'men']
    elig_funcs = [elig1, elig2]
    target_coverages = [0.8, 0.5]

    # create the intervention & supplies
    product = shot_1y_imperfect
    supplies = infinite_prep(product=product)
    intervention = SuppliedPrep(name='this-is-a-test', supplies=supplies, eligibilities=elig_funcs, coverages=target_coverages)

    diseases = [sti.HIV(beta_m2f=0.05, beta_m2c=0.1, init_prev=0.0)]
    eff_analyzer = PrepEfficacyAnalyzer()
    duration_analyzer = PrepDurationAnalyzer()
    analyzers = [eff_analyzer, duration_analyzer]
    sim = build_testing_sim(analyzers=analyzers, interventions=[intervention], n_agents=5000, duration=duration, diseases=diseases)

    sim.run()
    efficacy_results = sim.results[eff_analyzer.name][eff_analyzer.result_name]
    efficacy_tis = sim.results[eff_analyzer.name][eff_analyzer.ti_name]
    duration_results = sim.results[duration_analyzer.name][duration_analyzer.result_name]
    duration_tis = sim.results[duration_analyzer.name][duration_analyzer.ti_name]

    # ensure prep was distributed
    n_timeseries = len(efficacy_results.keys())
    if verbose:
        print(f"{n_timeseries} agent timeseries recorded")
    assert n_timeseries > 0, "No agent PrEP efficacy timeseries were recorded"

    n_efficacies = sum([len(efficacies) for efficacies in efficacy_results.values()])
    if verbose:
        print(f"{n_efficacies} agent efficacies recorded")
    assert n_efficacies > 0, "Only empty PrEP efficacy timeseries were recorded, which is a bug"

    # ensure we recorded appropriate unique efficacies and durations (also for debugging as needed)
    unique_efficacies = set(chain(*efficacy_results.values()))
    n_unique_efficacies = len(unique_efficacies)
    unique_durations = set(chain(*duration_results.values()))
    n_unique_durations = len(unique_durations)

    if verbose:
        print(f"{n_unique_efficacies} unique agent efficacies recorded: {unique_efficacies}")
        print(f"{n_unique_durations} unique agent durations recorded: {unique_durations}")

    n_expected = len(product.eff_by_ti)
    assert n_unique_efficacies == n_expected, f"{n_unique_efficacies} recorded but {n_expected} expected"
    assert n_unique_durations == n_expected, f"{n_unique_durations} recorded but {n_expected} expected"

    # ensure that the analyzers are recording the exact same tis for every uid -- consistency check
    assert efficacy_tis.keys() == duration_tis.keys(), "There should be (but is not) a 1-1 correspondence in duration and efficacies uids recorded."
    for uid, tis in efficacy_tis.items():
        other_tis = duration_tis[uid]
        assert tis == other_tis, f"Corresponding time-index series between for uid: {uid} does not match between duration and efficacies recorded."

    # now that we know we are comparing the same uids at the same times, combine analyzer results to ensure proper prep efficacy was recorded
    for uid, efficacies in efficacy_results.items():
        durations = duration_results[uid]
        expected_eff = [product.efficacy_at_ti(ti=duration) for duration in durations]
        assert efficacies == expected_eff, f"At least one efficacy does not match expectations for uid: {uid}"

    return sim

@sc.timer()
def test_supply_limited_prep():
    duration = 3 # years
    n_supply = 100
    supplies = limited_prep(product=shot_1y_perfect, quantity=n_supply)
    intervention = SuppliedPrep(name='this-is-a-test', supplies=supplies)
    target_coverage = intervention.coverages[0]  # should be 100%
    eligibilities = {'fsw': intervention.default_eligibilities[0]}  # there is only one
    count_analyzer = PrepCountsAnalyzer(eligibilities=eligibilities, consider_new_infections=True)
    coverage_analyzer = PrepCoverageAnalyzer(eligibilities=eligibilities, consider_new_infections=True)
    sim = build_testing_sim(analyzers=[count_analyzer, coverage_analyzer], interventions=[intervention],
                            n_agents=5000, duration=duration)
    sim.run()

    analyzer_results = sim.results[count_analyzer.name][count_analyzer.result_name]

    # Trimming off final reporting ti as starsim always (currently) runs +1 ti more than requested.
    analyzer_results = {group_name: coverage[:-1] for group_name, coverage in analyzer_results.items()}
    for group_name, n_prep in analyzer_results.items():
        assert len(n_prep) == (duration * 12)

    if verbose:
        for group_name, n_prep in analyzer_results.items():
            print(f"group: {group_name} n_on_prep: {n_prep}")

    # Count checks
    # There should be agents on prep with at all timesteps (given sim/product timescale) (for when supplies are available, which is first 2 years in this setup)
    for group_name, n_prep in analyzer_results.items():
        for i in range(2*12):  # 2 years
            assert n_prep[i] > 0, f"group: {group_name} No agents found on PrEP at ti: {i}. Agents should be on PrEP for all timesteps for years 1-2."
            assert n_prep[i] < n_supply,  f"group: {group_name} too many agents on prep at ti: {i}. Exceeded supply limit."
        for i in range (2*12+1, duration*12):
            # The final year(s) of should have 0 agents on PrEP because it ran out
            assert n_prep[i] == 0, f"group: {group_name} no agents shoud be on PrEP for the final year (3), but {n_prep[i]} agents were at ti: {i}"

    # coverage checks, ensure that coverages are not exceeded and drop after a certain point (reaching 0 in year 3)
    analyzer_results = sim.results[coverage_analyzer.name][coverage_analyzer.result_name]
    if verbose:
        for group_name, coverage in analyzer_results.items():
            print(f"group: {group_name} coverage: {coverage}")

    # year 1: near full coverage (100%)
    # year 2: about 50 +/- 8% coverage
    # year 3: 0 coverage
    for group_name, coverage in analyzer_results.items():
        target = target_coverage
        tolerance = 0.05
        for i in range(9, 1*12):
            # coverage should approach the target by the late part of year1
            assert coverage[i] == pytest.approx(expected=target, abs=tolerance), f"Group: {group_name} Expected coverage at ti: {i} == {target} +/- {tolerance}, was {coverage[i]}"

        target = 0.6
        tolerance = 0.1
        for i in range(1*12+9, 2*12):
            # coverage should approach the target by the late part of year2  -- reduced due to supply constraint
            assert coverage[i] == pytest.approx(expected=target, abs=tolerance), f"Group: {group_name} Expected coverage at ti: {i} == {target} +/- {tolerance}, was {coverage[i]}"

        for i in range(2*12, len(coverage)):
            assert coverage[i] == 0, f"Year 3 coverages should be 0, but ti: {i} coverage is: {coverage[i]}"

    return sim




if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    timer = sc.timer()

    test_default_eligibility_and_coverage()
    test_default_eligibility_and_modified_coverage()
    test_custom_eligibilities_and_coverages()
    test_prep_eff_updates_correctly()
    test_supply_limited_prep()

    sc.heading("Total:")
    timer.toc()

    if do_plot:
        plt.show()
