import random
from typing import Callable

import starsim as ss

from stisim.interventions.prep_manager import PrepManager
from stisim.interventions.supplied_intervention import SuppliedIntervention
from stisim.interventions.supplies import Supplies


class SuppliedPrep(SuppliedIntervention):

    default_eligibilities = [lambda sim: sim.networks.structuredsexual.fsw & (~sim.diseases.hiv.infected)]

    def __init__(self,
                 name: str,
                 eligibilities: list[Callable] = None,
                 coverages: list[float] = None,
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
        # setting eligibility default to uninfected FSW @100% coverage
        eligibilities = self.default_eligibilities if eligibilities is None else eligibilities
        super().__init__(name=name, supplies=supplies, eligibilities=eligibilities, coverages=coverages, *args, **kwargs)

        # This particular Intervention currently expects only one prep supply/product
        self.prep_supply = self.supplies.get_supply(prod_type='prep')[0]  # TODO: use of string 'prep' not very clean

        self.update_pars(pars, **kwargs)

        # Actual values for p will be overridden as needed
        self.care_seeking_dist = ss.bernoulli(p=0)
        self.care_seeking_dist.init(force=True)
        self.base_care_seeking_rate = 0.8  # TODO, arbitrary, probably a param?

    # def start_step(self):
    def update_efficacy(self):
        """identify agents on THIS prep intervention+product combo and update efficacy based on time"""
        sim = self.sim
        product = self.prep_supply.product
        PrepManager.update_eff(sim=sim, to_update=self.on_this_prep, product=product)

    def step(self):
        sim = self.sim
        hiv = sim.diseases.hiv
        supply = self.prep_supply
        product = supply.product

        self.update_efficacy()


        # identify new uptake
        eligibilities = self.eligibilities
        eligibilities = [self.check_eligibility(eligibility=eligibility, return_uids=False)
                         for eligibility in eligibilities]
        n_eligibles = [len(eligibilities[i].uids) for i in range(len(eligibilities))]

        # determine current prep coverage levels. Assuming ANY prep counts here (not just this intervention).
        cur_coverages = [len((hiv.on_prep & eligible).uids) / len(eligible.uids) for eligible in eligibilities]

        # determine which eligible agents are seeking care
        care_seeking = [self.determine_care_seeking(eligible=eligible) for eligible in eligibilities]

        # agents who can be offered to are those who are:
        # seeking care AND not on PrEP of any kind
        # ... OR are looking to reuptake this intervention
        offer_pools = [(seeking & (~hiv.on_prep).uids)
                       # | (eligibilities[i].uids & self.will_reuptake.uids)
                       for i, seeking in enumerate(care_seeking)]

        dist_func = lambda uptakers: PrepManager.uptake(sim=sim, to_uptake=uptakers, prep_intervention=self,
                                                        product=product, use_supplies=True)

        # Now select agents to get PrEP from the offer pools up to computed seeking/coverage/supply limit and distribute
        target_coverages = self.coverages
        selected = self.distribute(offer_pools=offer_pools,
                                   cur_coverages=cur_coverages, target_coverages=target_coverages,
                                   n_eligibles=n_eligibles, n_supply=supply.quantity, dist_func=dist_func)

    @property
    def on_this_prep(self):
        """on the PrEP from THIS intervention or not"""
        hiv = self.sim.diseases.hiv
        return hiv.on_prep & (hiv.prep_source == self._name_to_index)

    def determine_care_seeking(self, eligible):
        """determines whether eligible agents seek care for this intervention. Modifies param: eligible for return"""
        hiv = self.sim.diseases.hiv

        # TODO: we need to ensure this algorithm is appropriate; it is a stand-in for now.
        n_eligible = len(eligible.uids)
        self.care_seeking_dist.pars.p = hiv.care_seeking[eligible] * min(self.base_care_seeking_rate, 1.0)
        seeking_draw = self.care_seeking_dist.rvs(n=n_eligible)
        uids_seeking_care = eligible.uids[seeking_draw]

        return uids_seeking_care

    # Disabled pending researcher input
    # @property
    # def dropped_this_prep(self):
    #     """whether or not agents have dropped out of THIS intervention's PrEP by choice this ti"""
    #     ti = self.ti
    #     hiv = self.sim.diseases.hiv
    #     return (hiv.prep_source == self._name_to_index) & (hiv.ti_prep_drop == ti)

    # Disabled pending researcher input
    # @property
    # def will_reuptake(self):
    #     """whether or not agents are ready to try for reuptake of this intervention's product this ti"""
    #     ti = self.ti
    #     hiv = self.sim.diseases.hiv
    #     return self.on_this_prep & (hiv.ti_prep_end == ti) & (hiv.prep_n_reuptake > 0)
