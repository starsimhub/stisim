class Product:  # TODO: determine sublassing potential with ss.Product + features that can be pushed there, too
    def __init__(self, name, type, delivery_mode, cost, eff_by_ti, rel_eff_by_adherence=None):
        """
        A product that can be distributed to agents in interventions. They can be anything that modifies agent state
        or may be quantity limited (examples: doses of lenacapavir PrEP, "seats" in an educational campaign,
        doctor appointment slots)

        Args:
            name: (str) identifier of the product
            type: (str) category of the product, used to determine if 1+ Products have a similar function (for example,
                type "prep" applies to oral and injectable PrEP).
            delivery_mode: (str) How the product is delivered, used to distinguish Products that may have the same type
                but are delivered differently in intervention-relevant ways (for example, pill vs. shot)
            cost: (float) a per dose/distribution usage cost (units: currency)
            eff_by_ti: (list[float]) A effectivity durability map of the product in time, by time indices since
                delivery.
        """
        self.name = name # lenacapavir, ...
        self.type = type  # prep, art, ...
        self.delivery_mode = delivery_mode  # e.g., shot, pill, topical, educational, ...
        self.cost = cost # per dose/item
        self.eff_by_ti = eff_by_ti # TODO, usage of this should accept different formats than just by ti. Probably a Time/Date class-based representation?
        self.id = id(self)

    def efficacy_at_ti(self, ti):
        """efficacy at a given time index after initiation. 0 beyond product durability"""
        eff = self.eff_by_ti[ti] if ti < len(self.eff_by_ti) else 0
        return eff

    @property
    def max_durability(self):
        """The maximum number of timesteps for which the product has any efficacy"""
        return len(self.eff_by_ti)
