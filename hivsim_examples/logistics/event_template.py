from typing import Callable

import numpy as np
import starsim as ss

from stisim import SuppliedIntervention

# template variables and their settings
# product_name: {product_name}
# product_quantity: {product_quantity}
# event_name_downcased: {event_name_downcased}


class SuppliedInterventionSubclass(SuppliedIntervention):
    def distribute_event(self, uptakers: ss.uids, product_name: str, product_quantity_per_agent: int):
        # apply the effect of "sample event" to the uptakers, modifying their state

        # User note: define the effect of this intervention/event in here.

        #
        # This is an empty sample distribution function
        #

        # this prod_name matches the product_name for this event in the step() function.
        # this quantity (LHS of product) matches the product_quantity_per_agent for this event in the step() function.
        # Both are derived from the 'Supplies' dict key of the event, e.g.: 'Supplies': {'lenacapavir': 1}
        if product_name is not None:
            self.supplies.use(prod_name=product_name, quantity=product_quantity_per_agent*len(uptakers))

        results = np.zeros(len(uptakers), dtype=bool)  # sample return, but could be other dtype, always same length as uptakers arg
        return results

    def determine_care_seeking(self, eligible):
        # Filter method that can be used to apply health seeking behavior to eventB before distributing eventB.

        # User note: define any logic for determining a care seeking subset (if any) in here.

        #
        # This is an empty sample care seeking method that applies no filtering. All events should utilize a copy of
        # this method by default.
        #

        return eligible

    def step(self):
        offer_pools = [self.check_eligibility(eligibility=eligibility, return_uids=False) for eligibility in self.eligibilities]
        care_seeking = [self.determine_care_seeking(eligible=offer_pool) for offer_pool in offer_pools]
        # User note: if you wish to add target_coverage limits to this step (fraction of incoming agents), do so here.
        #  Alternatively, additional constraints can be applied (e.g., limiting global PrEP coverage of arbitrary
        #  demographics, if this were a PrEP distribution event. Doing may also require creating subsets of the offer
        #  pools, above, for passing to the event execution, below.
        target_coverages = None  # default
        cur_coverages = None  # default

        # User note: if the eligible fraction of 1+ offer_pool is less than everyone in the pool, adjust here.
        n_eligibles = [len(offer_pool) for offer_pool in offer_pools]

        # From the 'Supplies' dict key of the event, e.g.: 'Supplies': {'lenacapavir': 1}
        product_name = {product_name}
        product_quantity_per_agent = {product_quantity}

        event_selected, event_results = self.event(offer_pools=care_seeking, n_eligibles=n_eligibles,
                                                   product_name=product_name, product_quantity_per_agent=product_quantity_per_agent,
                                                   dist_func=self.distribute_event,
                                                   cur_coverages=cur_coverages, target_coverages=target_coverages)

        self.parent.{event_name_downcased}_selected = event_selected   # TODO: update to use setattr to allow consistent use of repr() in .format() call
        self.parent.{event_name_downcased}_results = event_results  # TODO: update to use setattr to allow consistent use of repr() in .format() call
        return
        # End template.
