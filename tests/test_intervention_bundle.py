"""
InterventionBundle tests.

Replicates tests/test_prep.py::test_default_eligibility_and_coverage, but wraps
the SuppliedPrep intervention in an InterventionBundle.
"""
import matplotlib.pyplot as plt
import sciris as sc
import sys

from pathlib import Path

import stisim as sti

from stisim.intervention_bundle import InterventionBundle
from stisim.logistics.product import Product
from stisim.logistics.supplies import Supplies
from stisim.logistics.supply import Supply

tests_directory = Path(__file__).resolve().parent
sys.path.append(str(tests_directory))

# from stisim.interventions.prep import SuppliedPrep
# from hiv_natural_history_analyzers import PrepCoverageAnalyzer

from testlib import build_testing_sim


verbose = False
do_plot = False
sc.options(interactive=False)

shot_1y_perfect = Product(name='test-shot', type='prep', delivery_mode='shot', cost=10, eff_by_ti=[1.0] * 12)


def infinite_prep(product):
    return limited_prep(product=product, quantity=1000000)


def limited_prep(product, quantity):
    return Supplies(supplies=[Supply(product=product, quantity=quantity)])


@sc.timer()
def test_default_eligibility_and_coverage_bundled():
    duration = 1  # years
    supplies = infinite_prep(product=shot_1y_perfect)
    intervention = SuppliedPrep(name='this-is-a-test', supplies=supplies)
    target_coverage = intervention.coverages[0]  # should be 100%
    eligibilities = {'fsw': intervention.default_eligibilities[0]}  # there is only one

    bundle = InterventionBundle(interventions=[intervention])

    analyzer = PrepCoverageAnalyzer(eligibilities=eligibilities, consider_new_infections=True)
    sim = build_testing_sim(analyzers=[analyzer], interventions=[bundle], n_agents=5000, duration=duration)
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


if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    timer = sc.timer()

    test_default_eligibility_and_coverage_bundled()

    sc.heading("Total:")
    timer.toc()

    if do_plot:
        plt.show()