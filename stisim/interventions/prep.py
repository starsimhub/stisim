import random
from typing import Iterable, Callable

import starsim as ss

from stisim.interventions.prep_manager import PrepManager
from stisim.interventions.supplied_intervention import SuppliedIntervention
from stisim.interventions.supplies import Supplies


class SuppliedPrep(SuppliedIntervention):

    default_eligibilities = [lambda sim: sim.networks.structuredsexual.fsw & (~sim.diseases.hiv.infected)]

    def __init__(self,
                 name: str,
                 eligibilities: Iterable[Callable] = None,
                 coverages: Iterable[float] = None,
                 supplies: Supplies = None,
                 pars=None, *args, **kwargs):
        """
        Representation of a standard PrEP distribution intervention with only one PrEP Supply/Product. The distribution
        of PrEP can be to 1+ eligibility groups at any target coverage 0.0-1.0 . Distribution will occur to agents
        that are eligible and determined to be care seeking. If coverage constraints are hit for one or more
        eligibility groups, distributions will not push any coverage beyond the target coverage. However, there is
        currently no mechanism for removing agents from PrEP, so coverage can float over target coverage if the
        eligibility group shrinks over time.

        If finite quantities of PrEP are modeled in the Supply of the PrEP product and quantity availablity restriction
        is reached, the remaining doses of PrEP will be distributed to eligibility groups proportional to the
        "coverage gap" of the groups. For example:
            groupA has target coverage 0.8 and actual coverage 0.6 (gap: 0.2)
            groupB has target coverage 0.5 and actual coverage 0.2 (gap 0.3)

            groupA will receive 40% (0.2 / (0.2 + 0.3)) of the remaining, limited PreP supply
            groupB will receive 60% (0.3 / (0.2 + 0.3)) of the remaining, limited PreP supply.


        Supplies object is expected to contain exactly one PrEP Supply of a Product.

        Args:
            name: (str) Intervention name. Useful for identification of and disambiguation of interventions.
            eligibilities: (list[callable]) Defines groups of agents who are eligible to enter the intervention.
                Default: one group: all HIV- FSW
            coverages: (list[float]) The target fraction of each corresponding eligibility group that will enter
                (be covered) by the intervention. Default: 1.0 per group.
            pars: (dict) intervention/module parameters to set
        """
        super().__init__(supplies=supplies, *args, **kwargs)
        self.name = name

        # This particular Intervention currently expects only one prep supply/product
        self.prep_supply = self.supplies.get_supply(prod_type='prep')[0]  # TODO: use of string 'prep' not very clean

        # setting eligibility default to uninfected FSW @100% coverage
        self.eligibilities = self.default_eligibilities if eligibilities is None else eligibilities
        self.coverages =     self.default_coverages     if coverages     is None else coverages

        if len(self.eligibilities) != len(self.coverages):
            raise ValueError(f"The number of eligibility groups {len(self.eligibilities)} must equal the number of "
                             f"provided coverage targets: {len(self.coverages)}")

        self.update_pars(pars, **kwargs)

        # Actual values for p will be overridden as needed
        self.care_seeking_dist = ss.bernoulli(p=0)
        self.care_seeking_dist.init(force=True)
        self.base_care_seeking_rate = 0.8  # TODO, arbitrary, probably a param?

    @property
    def default_coverages(self):
        return [1.0 for _ in range(len(self.eligibilities))]

    # def start_step(self):
    def update_efficacy(self):
        """identify agents on THIS prep intervention+product combo and update efficacy based on time"""
        sim = self.sim
        product = self.prep_supply.product
        PrepManager.update_eff(sim=sim, to_update=self.on_this_prep, product=product)
        return

    def step(self):
        sim = self.sim
        hiv = sim.diseases.hiv
        supply = self.prep_supply
        product = supply.product

        self.update_efficacy()

        # identify new uptake & reuptake
        eligibilities = [self.check_eligibility(eligibility=eligibility, return_uids=False)
                         for eligibility in self.eligibilities]

        # determine current prep coverage levels. Assuming ANY prep counts here (not just this intervention).
        cur_coverages = [len((hiv.on_prep & eligible).uids) / len(eligible.uids) for eligible in eligibilities]

        # determine the maximum PrEP distribution based on coverage targets of eligible populations
        # (absent supply constraints or care seeking)
        n_to_reach_coverage = [int((self.coverages[i] - cur_coverages[i]) * len(eligible.uids))
                       for i, eligible in enumerate(eligibilities)]
        n_to_reach_coverages = sum(n_to_reach_coverage)

        # determine which eligible agents are seeking care
        care_seeking = [self.determine_care_seeking(eligible=eligible) for eligible in eligibilities]

        # agents who can be offered to are those who are:
        # seeking care AND not on PrEP of any kind
        # ... OR are looking to reuptake this intervention
        offer_pools = [(seeking & (~hiv.on_prep).uids) |
                       (eligibilities[i].uids & self.will_reuptake.uids)
                       for i, seeking in enumerate(care_seeking)]
        offer_pool_size = [len(pool) for pool in offer_pools]
        n_seeking = sum(offer_pool_size)

        # The number of agents who go on PrEP will be the minimum constrained by: coverage, seeking behavior, supply
        n_offer = min(n_to_reach_coverages, n_seeking, supply.quantity)
        # determine the supply rationing ratio (to be applied proportionally to target coverages), if supply limited
        if n_offer > 0:
            # figure how who to give PrEP to
            if n_offer <= supply.quantity:
                coverages = self.coverages # no supply limits. Select agents for PrEP based on seeking & coverages, only
            else:
                # apply supply limits proportional to specified group coverages.
                # Computing the curtailed-coverage that is achievable based on supply (short of target coverage)
                supply_fraction = supply.quantity / n_offer
                coverages = [(coverage - cur_coverages[i]) * supply_fraction + cur_coverages[i]
                            for i, coverage in enumerate(self.coverages)]

            # recompute n_to_reach_coverage(s), since coverages MAY now be lower (due to supply)
            n_to_reach_coverage = [int((coverages[i] - cur_coverages[i]) * len(eligible.uids))
                                   for i, eligible in enumerate(eligibilities)]
            # if coverage has overshot (for any valid reason), do not put more agents on
            n_to_reach_coverage = [max(0, n) for n in n_to_reach_coverage]

            # Now determine if we are coverage or seeking limited -- PER POOL
            n_offer_to_pool = []
            for i, pool in enumerate(offer_pools):
                if n_to_reach_coverage[i] >= offer_pool_size[i]:
                    # seeking limited, all seekers go onto PrEP
                    n_offer_to_pool.append(offer_pool_size[i])
                else:
                    # coverage limited, we will have to select agents to go onto PrEP
                    n_offer_to_pool.append(n_to_reach_coverage[i])

            # Now select agents to get PrEP from the offer pools up to computed seeking/coverage/supply limit
            for i, offer_pool in enumerate(offer_pools):
                n = n_offer_to_pool[i]
                selected = random.sample(population=list(offer_pool), k=n)
                selected = ss.arrays.uids(selected)

                # classify the agents going onto PrEP
                reuptakers = selected & self.will_reuptake.uids
                uptakers = selected - reuptakers

                if len(uptakers) > 0:
                    PrepManager.uptake(sim=sim, to_uptake=uptakers, prep_intervention=self, product=product,
                                       allow_reuptake=self.allow_reuptake, use_supplies=True)
                if len(reuptakers) > 0:
                    PrepManager.reuptake(sim=sim, to_uptake=reuptakers, prep_intervention=self, product=product,
                                         use_supplies=True)

    @property
    def _name_to_index(self):
        """the intervention index of this intervention"""
        return self.sim.interventions.keys().index(self.name)


    @property
    def on_this_prep(self):
        """on the PrEP from THIS intervention or not"""
        hiv = self.sim.diseases.hiv
        return hiv.on_prep & (hiv.prep_source == self._name_to_index)

    @property
    def dropped_this_prep(self):
        """whether or not agents have dropped out of THIS intervention's PrEP by choice this ti"""
        ti = self.ti
        hiv = self.sim.diseases.hiv
        return (hiv.prep_source == self._name_to_index) & (hiv.ti_prep_drop == ti)

    @property
    def will_reuptake(self):
        """whether or not agents are ready to try for reuptake of this intervention's product this ti"""
        ti = self.ti
        hiv = self.sim.diseases.hiv
        return self.on_this_prep & (hiv.ti_prep_end == ti) & (hiv.prep_n_reuptake > 0)

    def determine_care_seeking(self, eligible):
        """determines whether eligible agents seek care for this intervention. Modifies param: eligible for return"""
        hiv = self.sim.diseases.hiv

        # TODO: we need to ensure this algorithm is appropriate; it is a stand-in for now.
        n_eligible = len(eligible.uids)
        self.care_seeking_dist.pars.p = hiv.care_seeking[eligible] * min(self.base_care_seeking_rate, 1.0)
        seeking_draw = self.care_seeking_dist.rvs(n=n_eligible)
        uids_seeking_care = eligible.uids[seeking_draw]

        return uids_seeking_care
