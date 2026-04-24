import starsim as ss


class PrepManager(ss.Intervention):
    """
    Helper class that primarily acts to kick agents off of PrEP before any other PrEP-related intervention runs step().
    This ensures that agents can switch (if eligible, care seeking, ...) between PrEP products and/or product sources
    with no coverage time gap. The same assurance applies for agents simply re-seeking the same product/source at the
    end of their prior product durability time. Related, common PrEP functionality related to uptake/efficacy updates
    are contained in here for centralization.
    """

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

        # For all agents on ANY prep, determine who is at the end of their current regimen and handle dropouts
        to_drop = (hiv.prep_n_reuptake < 1) & (hiv.ti_prep_end == ti)

        # these agents definitely drop out now
        n_to_drop = len(to_drop.uids)
        if n_to_drop > 0:
            self.dropout(sim=sim, to_drop=to_drop)

    @staticmethod
    def uptake(sim, to_uptake, prep_intervention, product, allow_reuptake=False, use_supplies=False):
        """
        agents defined by to_uptake will start the intervention, use a supply . use_supplies=False ignores supply usage
        and constraints. Agents will also decide if they will reup at durability end, if reuptake is allowed.
        """
        ti = sim.ti
        hiv = sim.diseases.hiv

        hiv.on_prep[to_uptake] = True
        hiv.prep_source[to_uptake] = prep_intervention._name_to_index
        hiv.prep_product[to_uptake] = product.id
        hiv.ti_prep_start[to_uptake] = ti
        hiv.ti_prep_end[to_uptake] = ti + product.max_durability
        eff_prep = product.efficacy_at_ti(ti=0)
        hiv.rel_sus[to_uptake] = hiv.rel_sus[to_uptake] * (1 - eff_prep)
        hiv.prep_eff[to_uptake] = eff_prep

        n_uptake = len(to_uptake)

        # Disabled pending further clarification/research request
        # if allow_reuptake:
        #     # Determine the max number of times an agent will return to refill/reenroll/reuptake from the same
        #     # intervention/supply
        #     # dist = ss.Dist(p=hiv.care_seeking[to_uptake])
        #     # dist.init(force=True)
        #     # hiv.prep_n_reuptake[to_uptake] = int(dist.rvs(n=n_uptake) * 10)  # one draw per to_uptake
        #
        #     hiv.prep_n_reuptake[to_uptake] = 3  # TODO: stand-in, only
        # else:
        hiv.prep_n_reuptake[to_uptake] = 0

        # update product supply
        if use_supplies:
            prep_intervention.use(prod_name=product.name, quantity=n_uptake)

    # Disabled pending further clarification/research request
    # @staticmethod
    # def reuptake(sim, to_uptake, prep_intervention, product, use_supplies=False):
    #     """
    #     agents defined by to_uptake will start the intervention, use a supply . use_supplies=False ignores supply usage
    #     and constraints. Agents will also decide if they will reup at durability end, if reuptake is allowed.
    #     """
    #     ti = sim.ti
    #     hiv = sim.diseases.hiv
    #
    #     hiv.on_prep[to_uptake] = True
    #     hiv.prep_source[to_uptake] = prep_intervention._name_to_index
    #     hiv.prep_product[to_uptake] = product.id
    #     hiv.ti_prep_start[to_uptake] = ti
    #     hiv.ti_prep_end[to_uptake] = ti + product.max_durability
    #     eff_prep = product.efficacy_at_ti(ti=0)
    #     # current_eff = hiv.prep_eff[to_uptake]
    #     hiv.rel_sus[to_uptake] = hiv.rel_sus[to_uptake] * (1 - eff_prep)
    #     hiv.prep_eff[to_uptake] = eff_prep
    #
    #     # use up a reuptake
    #     n_uptake = len(to_uptake)
    #     hiv.prep_n_reuptake[to_uptake] = hiv.prep_n_reuptake[to_uptake] - 1
    #
    #     # update product supply
    #     # supply.use(quantity=n_uptake)  # reduces supply, records reduction at this ti
    #     if use_supplies:
    #         prep_intervention.use(prod_name=product.name, quantity=n_uptake)

    @staticmethod
    def dropout(sim, to_drop):
        """agents defined by to_drop will immediately drop off of PrEP"""
        ti = sim.ti
        hiv = sim.diseases.hiv

        hiv.on_prep[to_drop] = False
        # hiv.rel_sus[to_drop] = hiv.rel_sus[to_drop] / (1 - hiv.prep_eff[to_drop])  # TODO: verify this is not needed during dropout due to supply constraints.
        hiv.prep_eff[to_drop] = None
        hiv.ti_prep_drop[to_drop] = ti
        hiv.prep_n_reuptake[to_drop] = 0

    @staticmethod
    def update_eff(sim, to_update, product):
        """
        updates the efficacy of PrEP for all specified agents to the time-index--dependent efficacy of the given
        product
        """
        ti = sim.ti
        hiv = sim.diseases.hiv
        to_update_uids = to_update.uids

        tis_on_prep = ti - hiv.ti_prep_start
        for ti_on_prep in range(product.max_durability):
            cohort = (tis_on_prep == ti_on_prep).uids & to_update_uids

            new_eff = product.efficacy_at_ti(ti=ti_on_prep)
            # old_eff = product.efficacy_at_ti(ti=ti_on_prep-1)

            # rescale rel_sus based on new efficacy
            hiv.rel_sus[cohort] = hiv.rel_sus[cohort] * (1-new_eff)
            hiv.prep_eff[cohort] = new_eff
        return
