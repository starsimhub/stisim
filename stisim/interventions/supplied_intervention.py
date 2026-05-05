from typing import Callable

import starsim as ss

from stisim.interventions.supplies import Supplies


class SuppliedIntervention(ss.Intervention):
    def __init__(self,
                 name: str,
                 eligibilities: dict[str, list[Callable]],
                 coverages: dict[str, list[float]] = None,
                 supplies: Supplies=None, pars=None, *args, **kwargs):
        """
        A SuppliedIntervention can represent 1+ sets of decisions that represent an 'event' containing 1+ 'decisions',
        each with its own set of eligibility criteria and corresponding coverages.

        Representation of a generic distribution intervention with an arbitrary number Supply/Product available for
        1+ number of intervention 'decisions'. Each decision is associated with 1+ eligibility groups at corresponding
        target coverages 0.0-1.0 . Distribution during any decision will occur to agents that are eligible and
        determined to be care seeking. Default care-seeking is all eligible agents. If coverage constraints are hit for
        one or more eligibility groups for a decision, distributions will not push any coverage beyond the target
        coverage. However, there is currently no mechanism for removing agents from an intervention type (e.g. PrEP), so
        coverage can float over target coverage if the eligibility group shrinks over time.

        If finite quantities of a product are modeled in the Supplies and quantity availability restriction
        is reached, the remaining supply counts will be distributed to eligibility groups proportional to the
        "coverage gap" of the groups. For example:
            groupA has target coverage 0.8 and actual coverage 0.6 (gap: 0.2)
            groupB has target coverage 0.5 and actual coverage 0.2 (gap 0.3)

            groupA will receive 40% (0.2 / (0.2 + 0.3)) of the remaining, limited PreP supply
            groupB will receive 60% (0.3 / (0.2 + 0.3)) of the remaining, limited PreP supply.

        For example, a SuppliedIntervention representing a doctor office visit that combines the potential for a HIV
        test and PrEP distribution may utilize the following:

        eligibilities = {'hiv_test': [elig_func1], 'prep': [elig_func2, elig_func3]}
        coverages     = {'hiv_test': [cov1],       'prep': [cov2,       cov3]}

        The step() function would likely use these as follows:
        First, check for agents eligible for an HIV test, delivering an appropriate number to satisfy the target
        coverage as best as possible.

        Second, after any test results have come back, the doctor offers PrEP to eligible agents. There are two different
        criteria here defined by elig_func2 and elig_func3, each with their on coverage targets. elig_func2 might
        specify pregnant women and elig_func3 might specify female sex workers. Note that agents could potentially
        qualify under more than one criterion. Agents matching each criterion will be offered PrEP up to corresponding
        coverage targets.

        Args:
            name: (str) Intervention name. Useful for identification of and disambiguation of interventions.
            eligibilities: (dict[str, list[Callable]]) Defines groups of agents who are eligible at various decision
                point(s) in the intervention. The first layer of key: [funcs] represents organization by decisions where
                the key is the name of the decision for use in step(). The [funcs] represent the functions that
                identify eligible agents.
            coverages: (dict[str, list[Callable]]) The target fraction of each corresponding eligibility group that will
                enter (be covered) at each decision in the intervention. Default: 1.0 per group for all decisions.
            supplies: (Supplies): a Supplies object representing a discrete quantity of 1+ Supply of Products. Can
                be used to in conjunction with or independently of coverage-limiting distribution. Default: No Supplies
                (or limits of them)
            pars: (dict) intervention/module parameters to set


        """
        super().__init__(*args, **kwargs)
        self.name = name
        self.supplies = Supplies() if supplies is None else supplies
        self.eligibilities = eligibilities
        self.coverages = self.default_coverages if coverages is None else coverages

        # if self.eligibilities.keys() != self.coverages.keys():
        #     raise ValueError(f"The number of eligibility types {len(eligibilities.keys())} must equal the number of "
        #                      f"provided coverage types: {len(self.coverages.keys())}")

        if len(self.eligibilities) != len(self.coverages):
            raise ValueError(f"The number of eligibility groups {len(self.eligibilities)} must equal the number of "
                             f"provided coverage targets: {len(self.coverages)}")

        self.accrued_cost = 0  # for supply use, distribution, and anything else a subclass wants to model
        # self.allow_reuptake = False  # disabled for now pending further clarification/research request

        self.update_pars(pars, **kwargs)

        states = [
            # TODO: ensure this defaults to False
            ss.BoolState('on_intervention'),  # T/F: are agents actively on/in this intervention in some way?
        ]
        self.define_states(*states)


    @property
    def _name_to_index(self):
        """the intervention index of this intervention"""
        return self.sim.interventions.keys().index(self.name)

    @property
    def default_coverages(self):
        return [1.0 for _ in range(len(self.eligibilities))]

    def determine_care_seeking(self, eligible):
        """defaults to all eligible agents. Override if desired to restrict or tailor this."""
        return eligible.uids

    def use(self, prod_name, quantity):
        """record cost to intervention for supplies used in addition to using up the quantity specified"""
        remaining_quantity, cost = self.supplies.use(prod_name=prod_name, quantity=quantity)
        self.accrued_cost += cost
        return remaining_quantity, cost

    def calc_supply_distribution(self, offer_pools: list, cur_coverages: list, target_coverages: list,
                                 n_eligibles: list, n_supply):
        # agents who can be offered to are those who are:
        # seeking care AND not on PrEP of any kind
        # ... OR are looking to reuptake this intervention
        offer_pool_size = [len(pool) for pool in offer_pools]

        coverage_gaps = [(coverage - cur_coverages[i]) for i, coverage in enumerate(target_coverages)]
        coverage_gaps = [max(gap, 0) for gap in coverage_gaps]  # we ignore negative gaps. They can happen due to pop changes, but should not affect OTHER coverage targets.
        total_coverage_gaps = sum(coverage_gaps)
        if total_coverage_gaps > 0:
            cov_gap_proportion = [cov_gap / total_coverage_gaps for cov_gap in coverage_gaps]
            n_gaps = [cov_gap * n_eligibles[i] for i, cov_gap in enumerate(coverage_gaps)]

            # This is the maximum *number* of distributions per group that can occur, taking into account:
            # - count gap to target coverage
            # - offer pool size
            # - coverage gap -based proportionality of the remaining supply
            max_supply_dist = [int(min(n_gap, offer_pool_size[i], cov_gap_proportion[i] * n_supply))
                               for i, n_gap in enumerate(n_gaps)]
        else:
            max_supply_dist = [0] * len(coverage_gaps)  # nothing to distribute, all target coverages met
        return max_supply_dist

    def distribute(self, offer_pools, cur_coverages, target_coverages, n_eligibles, n_supply, dist_func):
        import random  # TODO: change this to use a ss-based selection

        # determine how many agents will be distributed to in each offer pool
        supply_dists = self.calc_supply_distribution(offer_pools=offer_pools,
                                                     cur_coverages=cur_coverages, target_coverages=target_coverages,
                                                     n_eligibles=n_eligibles, n_supply=n_supply)

        all_selected = ss.arrays.uids()
        for i, offer_pool in enumerate(offer_pools):
            selected = random.sample(population=list(offer_pool), k=supply_dists[i])
            selected = ss.arrays.uids(selected)

            if len(selected) > 0:
                dist_func(uptakers=selected)
                all_selected = all_selected.union(selected)
        return all_selected