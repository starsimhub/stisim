from stisim.logistics.delivery_mode import DeliveryMode
from stisim.logistics.product_category import ProductCategory


class Product:  # TODO: determine sublassing potential with ss.Product + features that can be pushed there, too
    def __init__(self, name: str, category: ProductCategory, delivery_mode: DeliveryMode, cost: float, eff_by_ti: list[float]):
        """
        A product that can be distributed to agents in interventions. They can be anything that modifies agent state
        or may be quantity limited (examples: doses of lenacapavir PrEP, "seats" in an educational campaign,
        doctor appointment slots)

        Args:
            name: (str) identifier of the product
            category: (ProductCategory) category of the product, used to determine if 1+ Products have a similar
                function (for example, ProductCategory.PREP applies to oral and injectable PrEP). Must be a
                ProductCategory member; anything else (including its raw string value) raises ValueError at
                construction.
            delivery_mode: (DeliveryMode) How the product is delivered, used to distinguish Products that may
                have the same category but are delivered differently in intervention-relevant ways (for example,
                pill vs. shot). Must be a DeliveryMode member; anything else (including its raw string value) raises
                ValueError at construction.
            cost: (float) a per dose/distribution usage cost (units: currency)
            eff_by_ti: (list[float]) A effectivity durability map of the product in time, by time indices since
                delivery.
        """
        self.name = name # lenacapavir, ...
        self.category = self._validate(category, ProductCategory, 'category')  # ProductCategory.PREP, ...
        self.delivery_mode = self._validate(delivery_mode, DeliveryMode, 'delivery_mode')  # DeliveryMode.SHOT, ...
        self.cost = cost # per dose/item
        self.eff_by_ti = eff_by_ti # TODO, usage of this should accept different formats than just by ti. Probably a Time/Date class-based representation?

    @staticmethod
    def _validate(value, enum_cls, field):
        """Require an enum member; raw strings (and anything else) are rejected with a helpful ValueError."""
        if not isinstance(value, enum_cls):
            valid = [member.value for member in enum_cls]
            raise ValueError(f"{field} must be a {enum_cls.__name__} member, got {value!r}. Valid options: {valid}.")
        return value

    def efficacy_at_ti(self, ti):
        """efficacy at a given time index after initiation. 0 beyond product durability"""
        if ti < 0:
            raise IndexError(f"ti must be a non-negative time index since delivery, got {ti}.")
        eff = self.eff_by_ti[ti] if ti < len(self.eff_by_ti) else 0
        return eff

    @property
    def max_durability(self):
        """The maximum number of timesteps for which the product has any efficacy"""
        return len(self.eff_by_ti)
