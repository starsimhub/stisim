import starsim as ss

"""
utilizes these attributes on hiv module:


on_prep: bool, actively on effective prep?
prep_source: str, the intervention name that delivered prep most recently
prep_product: str, the product name most recently used
ti_prep_start: int, ti of most recent prep start.
ti_prep_end: int, ti when most recent prep is no longer active. Also ti for potential reuptake/start new product.
eff_prep: float, current effectiveness of prep
prep_reuptake: bool, whether an agent will try to reuptake prep at ti_prep_end
ti_prep_drop: int, ti of most recent prep dropout

rel_trans: (float) prep modifies this by its effectiveness rating
healthcare_engagement: (float) likelihood to follow-through with a recommended healthcare step
"""

class PrepManager(ss.Intervention):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.interventions = []


    def step(self):
        pass

    def start_step(self):
        """
        This is the core requirement for utilizing a PrepManager intervention; dropping agents so they can reup on
        another intervention and/or product without a ti with no coverage (unless they choose or are forced to do so).
        """
        ti = self.ti
        sim = self.sim
        hiv = self.sim.diseases.hiv

        # For all agents on prep, determine who is at the end of their current regimen and handle dropouts
        # may_reup = hiv.prep_reuptake & (hiv.ti_prep_end == ti)
        to_drop = (~hiv.prep_reuptake) & (hiv.ti_prep_end == ti)

        # these agents definitely drop out now
        self.dropout(sim=sim, to_drop=to_drop)

    @staticmethod
    def uptake(sim, to_uptake, prep_intervention, product, allow_reuptake, use_supplies=False):
        """
        agents defined by to_uptake will start the intervention, use a supply . use_supplies=False ignores supply usage
        and constraints. Agents will also decide if they will reup at durability end, if reuptake is allowed.
        """
        ti = sim.ti
        hiv = sim.diseases.hiv

        hiv.on_prep[to_uptake] = True
        hiv.prep_source[to_uptake] = prep_intervention.name
        hiv.ti_prep_start[to_uptake] = ti
        hiv.ti_prep_end[to_uptake] = product.max_durability
        eff_prep = product.efficacy_at_ti(ti=0)
        hiv.rel_trans[to_uptake] = hiv.rel_trans[to_uptake] * (1 - eff_prep)
        hiv.eff_prep[to_uptake] = eff_prep

        # determine per-agent if they will try to reuptake/reenroll/refill
        n_uptake = len(to_uptake.uids)
        if allow_reuptake:
            # TODO: are reuptakers more likely to return than new agents are to go the first time?
            dist = ss.bernoulli(p=hiv.care_seeking[to_uptake])
            hiv.prep_reuptake[to_uptake] = dist.rvs(n=n_uptake)  # one draw per to_uptake

        # update product supply
        # supply.use(quantity=n_uptake)  # reduces supply, records reduction at this ti
        if use_supplies:
            prep_intervention.use(product_name=product.name, quantity=n_uptake)

    @staticmethod
    def dropout(sim, to_drop):
        """agents defined by to_drop will immediately drop off of PrEP"""
        ti = sim.ti
        hiv = sim.diseases.hiv

        hiv.on_prep[to_drop] = False
        hiv.rel_trans[to_drop] = hiv.rel_trans[to_drop] / (1 - hiv.eff_prep[to_drop])
        hiv.eff_prep[to_drop] = None
        hiv.ti_prep_drop[to_drop] = ti

    @staticmethod
    def update_eff(sim, to_update, product):
        """ updates the efficacy of PrEP for all specified agents to the time-dependent efficacy of the given product"""
        ti = sim.ti
        hiv = sim.diseases.hiv

        tis_on_prep = ti - hiv.ti_prep_start[to_update]
        for ti_on_prep in range(1, product.max_durability-1):   # TODO: measured in ti
            cohort = tis_on_prep == ti_on_prep

            new_eff = product.efficacy_at_ti(ti=ti_on_prep)  # TODO: update efficacy_at_ti to be 0 at-and-beyond max_durability?
            old_eff = product.efficacy_at_ti(ti=ti_on_prep-1)

            # rescale rel_trans based on altered efficacy
            hiv.rel_trans[cohort] = hiv.rel_trans[cohort] / (1 - old_eff) * (1 - new_eff)
            hiv.eff_prep[cohort] = new_eff  # informational
