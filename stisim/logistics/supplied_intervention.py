import abc
from typing import Callable

import starsim as ss

from stisim.logistics.supplies import Supplies


class SuppliedIntervention(ss.Intervention, metaclass=abc.ABCMeta):
    """
    Abstract base for distribution interventions that hand out supply-limited (or unlimited) Products to eligible,
    care-seeking agents.

    This class provides the building blocks for a distribution intervention but is NOT runnable on its own: it does
    not implement step(), so it cannot be instantiated directly. Concrete subclasses MUST implement step() to
    orchestrate, each timestep, the flow:

        eligibilities -> determine_care_seeking -> calc_supply_distribution -> distribute -> use

    Subclassing contract:
        - step(self): REQUIRED. Drive the per-timestep distribution using the building-block methods below.
        - determine_care_seeking(self, eligible): OPTIONAL override; defaults to all eligible agents.
        - the dist_func passed to distribute(): a callback that applies the intervention's effect to selected
          agents and draws down supply via use() (see distribute() for its signature).
        - the on_intervention BoolState is defined here for subclasses to track active participation; the base
          neither reads nor writes it.

    Provided building blocks: determine_care_seeking(), calc_supply_distribution(), distribute(), use().
    """

    def __init__(self,
                 name: str,
                 eligibilities: list[Callable],
                 coverages: list[float] = None,
                 supplies: Supplies=None, pars=None, *args, **kwargs):
        """
        A SuppliedIntervention an intervention 'event' 1+ eligibility criteria and corresponding target coverages.
        Supplies can either be finite or infinite, depending on use case and usage costs are accrued

        Representation of a generic distribution intervention with an arbitrary number Supply/Product available for
        distribution. Each eligibility group is associated with its own target coverage at corresponding 0.0-1.0 .
        Distribution will occur to agents that are eligible and determined to be care seeking. Default care-seeking is
        all eligible agents. If coverage constraints are hit for one or more eligibility groups, distributions will not
        push any coverage beyond the target coverage. However, there is currently no mechanism for removing agents from
        an intervention type (e.g. PrEP), so coverage can float over target coverage if the eligibility group shrinks
        over time.

        If finite quantities of a product are modeled in the Supplies and quantity availability restriction
        is reached, the remaining supply counts will be distributed to eligibility groups proportional to the
        "coverage gap" of the groups. For example:
            groupA has target coverage 0.8 and actual coverage 0.6 (gap: 0.2)
            groupB has target coverage 0.5 and actual coverage 0.2 (gap 0.3)

            groupA will receive 40% (0.2 / (0.2 + 0.3)) of the remaining, limited supply
            groupB will receive 60% (0.3 / (0.2 + 0.3)) of the remaining, limited supply.

        Args:
            name: (str) Intervention name. Useful for identification of and disambiguation of interventions.
            eligibilities: (list[Callable]) Defines groups of agents who are eligible for the intervention. The
                [funcs] represent the functions that identify eligible agents and take one argument: sim .
            coverages: (list[float]) The target fraction of each corresponding eligibility group to cover by by the
                intervention. Default: 1.0 per eligibility function.
            supplies: (Supplies): a Supplies object representing a discrete quantity of 1+ Supply of Products. Can
                be used to in conjunction with or independently of coverage-limiting distribution. The same Supplies
                object may be passed to more than one SuppliedIntervention, in which case it is a shared resource:
                quantity drawdown and pooled cost (Supplies.accrued_cost) are aggregated across all sharers. Default:
                No Supplies (or limits of them)
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

        # Cost attributable to THIS intervention's own use() calls (supply use, distribution, and anything else a
        # subclass wants to model). When self.supplies is shared across interventions, this is only this
        # intervention's share; query self.supplies.accrued_cost for the shared pool total across all sharers.
        self.accrued_cost = 0
        # self.allow_reuptake = False  # disabled for now pending further clarification/research request

        self.update_pars(pars, **kwargs)

        states = [
            # TODO: ensure this defaults to False
            ss.BoolState('on_intervention'),  # T/F: are agents actively on/in this intervention in some way?
        ]
        self.define_states(*states)

    @abc.abstractmethod
    def step(self):
        """
        Advance the intervention by one timestep. REQUIRED in subclasses; the sim calls this each step.

        A typical implementation: identify eligible agents via self.eligibilities, narrow to care-seekers with
        determine_care_seeking(), size per-group allocations with calc_supply_distribution(), then pass the offer
        pools to distribute() with a dist_func that applies the effect and draws down supply via self.use().
        """
        raise NotImplementedError

    @property
    def default_coverages(self):
        return [1.0 for _ in range(len(self.eligibilities))]

    def determine_care_seeking(self, eligible):
        """defaults to all eligible agents. Override if desired to restrict or tailor this."""
        return eligible.uids

    def use(self, prod_name, quantity):
        """
        Use the specified quantity of a product, accruing its cost to this intervention and returning
        (remaining_quantity, cost).

        The cost is added to self.accrued_cost (this intervention's own spend). If self.supplies is shared with other
        interventions, those interventions' costs are not reflected here; query self.supplies.accrued_cost for the
        pooled total across all sharers.
        """
        remaining_quantity, cost = self.supplies.use(prod_name=prod_name, quantity=quantity)
        self.accrued_cost += cost
        return remaining_quantity, cost

    def calc_supply_distribution(self, offer_pools: list, cur_coverages: list, target_coverages: list,
                                 n_eligibles: list, n_supply):
        # offer_pools are the agents eligible to be offered the product this step (seeking care and not already on
        # this intervention, or seeking to re-uptake it); assembled by the caller. Here we only size each pool.
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
        """
        Select recipients across the offer pools and apply the intervention's effect to them.

        Args:
            offer_pools: (list[ss.uids]) per eligibility group, the agents who may be offered the product this
                step (typically care-seekers not already covered, plus any seeking to re-uptake). The caller
                (a subclass's step()) is responsible for assembling these.
            cur_coverages: (list[float]) current coverage per eligibility group.
            target_coverages: (list[float]) target coverage per eligibility group.
            n_eligibles: (list[int]) number of eligible agents per group.
            n_supply: (int | float) units available to distribute this step (np.inf for unlimited).
            dist_func: (Callable) callback applied to each group's selected recipients, invoked as
                dist_func(uptakers=<ss.uids>) -> None. It applies the product's effect to the given agents (e.g.
                setting state, flipping on_intervention) and is expected to draw down supply via self.use().

        Returns:
            ss.uids of all agents selected across all offer pools.
        """
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
