import starsim as ss

from stisim.interventions.prep_manager import PrepManager


class Product:
    def __init__(self, name, type, delivery_mode, cost, eff_by_ti, rel_eff_by_adherence):
        self.name = name # lenacapavir, ...
        self.type = type  # prep, art, ...
        self.delivery_mode = delivery_mode  # e.g., shot, pill, topical, ... TODO: not currently used
        self.cost = cost # per dose/item
        self.eff_by_ti = eff_by_ti
        self.rel_eff_by_adherence = rel_eff_by_adherence  # TODO: not used

    def efficacy_at_ti(self, ti):
        """efficacy at a given time index after initiation. 0 beyond product durability"""
        # TODO, usage of this should accept different formats
        eff = self.eff_by_ti[ti] if ti < len(self.eff_by_ti) else 0
        return eff

    @property
    def max_durability(self):
        # TODO, usage of this should accept different formats
        return len(self.eff_by_ti)


class Supply:
    class InsufficientSupplyException(ValueError):
        pass

    def __init__(self, quantity, product):
        self.product = product
        self.quantity = quantity
        self.accrued_cost = 0

    def use(self, quantity):
        remaining = self.quantity - quantity
        if remaining < 0:
            raise self.InsufficientSupplyException(f"Attempted to use {quantity} {self.product.name}. "
                                                   f"Only {self.quantity} exists.")
        self.quantity = remaining

        cost = quantity * self.product.cost
        self.accrued_cost += cost

        return cost

    def add(self, quantity):
        self.quantity += quantity
        return self.quantity


class Supplies:
    class MissingSupplyException(ValueError):
        pass

    class DuplicateSupplyException(ValueError):
        pass

    def __init__(self, supplies=None):
        self.supplies = [] if supplies is None else supplies
        self._supplies_by_name = {supply.product.name: supply for supply in self.supplies}
        if len(self._supplies_by_name) < len(self.supplies):
            raise self.DuplicateSupplyException(f"Two or more products with the same name were added to a Supplies. "
                                                f"Product names must be unique.")
        self._supplies_by_type = {}
        for supply in self.supplies:
            typ = supply.product.type
            if typ not in self._supplies_by_type:
                self._supplies_by_type[typ] = []
            self._supplies_by_type[typ].append(supply)

    @property
    def accrued_cost(self):
        cost = sum([supply.accrued_cost for supply in self.supplies])
        return cost

    @property
    def product_names(self):
        return self._supplies_by_name.keys()

    def has_product(self, name):
        return name in self._supplies_by_name

    def get_quantity(self, name):
        return self.get_supply(prod_name=name).quantity

    def get_product(self, name):
        return self.get_supply(prod_name=name).product

    def get_supply(self, prod_name=None, prod_type=None):
        if not ((prod_name is None) ^ (prod_type is None)):
            raise ValueError(f"Exactly one of prod_name and prod_type must be specified for method "
                             f"Supplies.get_supply()")
        if prod_type is None:
            # prod_name was specified: a single product will be returned
            requested = self._supplies_by_name.get(prod_name, None)
            if requested is None:
                raise self.MissingSupplyException(f"Supply of product of name: {prod_name} requested but is "
                                                  f"not included.")
        else:
            # prod_type was specified: a list of 1+ products will be returned
            requested = self._supplies_by_type.get(prod_name, None)
            if requested is None:
                raise self.MissingSupplyException(f"Supplies of product of type: {prod_type} requested but are "
                                                  f"not included.")
        return requested

    def __repr__(self):
        strs = [f"{supply.product.name}: {supply.quantity}" for supply in self.supplies]
        result = '\n'.join(strs)
        return result


class SuppliedIntervention(ss.Intervention):
    def __init__(self, supplies: Supplies=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.supplies = Supplies() if supplies is None else supplies
        self.accrued_cost = 0  # for supply use, distribution, and anything else a subclass wants to model

    def use(self, product_name, quantity):
        """record cost to intervention for supplies used in addition to using up the quantity specified"""
        supply = self.supplies.get_supply(prod_name=product_name)
        cost = supply.use(quantity=quantity)
        self.accrued_cost += cost


"""
Questions:
- how to model reuptake: currently, this is a "per-durability" check
-- consider a person to have dropped PrEP at end of durability (even if they come back for more PrEP)?
-- Supply-based vs coverage-based & uptake speed: (n_doses + "event"-style ramp up) vs. grab-everyone-and-put-them-on-now
--- do we need to worry about uptake ramping, all doses unrealistically delivered in 1 ti? Is this a historical vs. future detail in a country model?
-- how to model willingness? Using hiv.care_seeking as likelihood of being willing to continue through this intervention.
- Define: prep drop vs prep durability end? (necessary?)
"""
class SuppliedPrep(SuppliedIntervention):

    def __init__(self, name, allow_reuptake=False, eligibility=None, pars=None, *args, **kwargs):
        # currently expects exactly ONE supply object in supplies

        # setting eligibility default to uninfected FSW
        if eligibility is None:
            eligibility = lambda sim: sim.networks.structuredsexual.fsw & (~sim.diseases.hiv.infected)
        super().__init__(eligibility=eligibility, *args, **kwargs)

        self.define_pars(
            eff_intervention=0.8
        )
        self.update_pars(pars, **kwargs)

        self.name = name  # prep intervention object disambiguation ("clinic_checkup", "community_test_and_treat", ...)
        self.allow_reuptake = allow_reuptake
        self.supply = self.supplies.get_supply(prod_type='prep')[0]  # This particular Intervention expects only one

    def start_step(self):
        # identify agents on THIS prep intervention and update efficacy based on time
        sim = self.sim
        product = self.supply.product
        PrepManager.update_eff(sim=sim, to_update=self.on_this_prep, product=product)

    def step(self):
        # NOTE: this is a supply/effectiveness-driven, not a coverage-driven, intervention. They have overlaps.
        sim = self.sim
        supply = self.supply
        product = supply.product

        # identify new uptake & reuptake

        # eligible agents are those targeted, including those who reached the end of this intervention product durability this ti
        eligible = self.check_eligibility()
        care_seeking = self.apply_care_seeking(eligible=eligible)  # determine which eligible agents are seeking care

        # offer to those who:
        # 1) are seeking care AND did not drop this intervention this ti
        # 2) are covered by this intervention, at durability end, and returning for reuptake
        offer_pool = (care_seeking & ~self.dropped_this_prep) | self.will_reuptake

        # program chance to find an eligible, care seeking agent
        offered = self.apply_effectiveness(seeking=offer_pool)

        # apply supply constraints
        uptake = self.apply_supply_limit(offered=offered, limit=supply.quantity)

        # apply new uptake and reuptake, record usage/costs
        PrepManager.uptake(sim=sim, to_uptake=uptake, prep_intervention=self, product=product,
                           allow_reuptake=self.allow_reuptake)

        # now check for reuptakers that failed to reup due to supply constraints. Force them to drop this prep now.
        # Note that self.will_reuptake now excludes those who successfully executed uptake above.
        PrepManager.dropout(sim=sim, to_drop=self.will_reuptake)

    @property
    def on_this_prep(self):
        """on the PrEP from THIS intervention or not"""
        hiv = self.sim.diseases.hiv
        return hiv.on_prep & (hiv.prep_source == self.name)

    @property
    def dropped_this_prep(self):
        """whether or not agents have dropped out of THIS intervention's PrEP by choice this ti"""
        ti = self.ti
        hiv = self.sim.diseases.hiv
        return (hiv.prep_source == self.name) & (hiv.ti_prep_drop == ti)

    @property
    def will_reuptake(self):
        """whether or not agents are ready to try for reuptake of this intervention's product this ti"""
        ti = self.ti
        hiv = self.sim.diseases.hiv
        return self.on_this_prep & (hiv.ti_prep_end == ti)  # assumes non-reuptakers have dropped (hiv.prep_reuptake is redundant check)

    def apply_care_seeking(self, eligible):
        """determines whether agents seek care ((re-)enter this intervention)"""
        hiv = self.sim.diseases.hiv
        n_eligible = len(eligible.uids)
        dist = ss.bernoulli(p=hiv.care_seeking[eligible])
        care_seeking = dist.rvs(n=n_eligible)
        return care_seeking

    def apply_effectiveness(self, seeking):
        """determiness whether agents are reached by this intervention"""
        n_seeking = len(seeking.uids)
        dist = ss.bernoulli(p=self.pars.eff_intervention)
        reached = dist.rvs(n=n_seeking)
        return reached

    def apply_supply_limit(self, offered, limit):
        """
        basic function to determine who is offered the intervention if a supply constraint is reached.
        This basic version is just random, but could be overridden with other metrics.
        """
        import random

        uids = offered.uids
        n_offered = len(uids)
        if n_offered > limit:
            chosen = random.sample(population=uids, k=limit)
            offered = chosen
        return offered
