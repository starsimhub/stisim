import starsim as ss

from stisim.interventions.supplies import Supplies


class SuppliedIntervention(ss.Intervention):
    def __init__(self, supplies: Supplies=None, *args, **kwargs):
        """

        Args:
            supplies: (Supplies): a Supplies object representing a discrete quantity of 1+ Supply of Products. Can
                be used to in conjunction with or independently of coverage-limiting distribution. Default: No Supplies
                (or limits of them)
        """
        super().__init__(*args, **kwargs)
        self.supplies = Supplies() if supplies is None else supplies
        self.accrued_cost = 0  # for supply use, distribution, and anything else a subclass wants to model
        self.allow_reuptake = False  # disabled for now pending further clarification/research request

    def use(self, prod_name, quantity):
        """record cost to intervention for supplies used in addition to using up the quantity specified"""
        remaining_quantity, cost = self.supplies.use(prod_name=prod_name, quantity=quantity)
        self.accrued_cost += cost
        return remaining_quantity, cost