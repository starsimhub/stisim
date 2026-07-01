"""Define the Product class: a distributable item with time-varying efficacy and a per-unit cost."""
import warnings

from stisim.logistics.delivery_mode import DeliveryMode
from stisim.logistics.product_category import ProductCategory

__all__ = ["Product"]


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
            cost: (float) a per dose/distribution usage cost (units: currency). Negative costs are permitted but
                emit a warning, as they are unusual.
            eff_by_ti: (list[float]) A effectivity durability map of the product in time, by time indices since
                delivery. Must be non-empty and every entry must be a fraction in [0.0, 1.0].
        """
        self.name = name # lenacapavir, ...
        self.id = id(self)
        self.category = self._validate(category, ProductCategory, 'category')  # ProductCategory.PREP, ...
        self.delivery_mode = self._validate(delivery_mode, DeliveryMode, 'delivery_mode')  # DeliveryMode.SHOT, ...

        if cost < 0:
            warnings.warn(f"Product cost is negative ({cost}); this is unusual but permitted.")
        self.cost = cost # per dose/item

        if len(eff_by_ti) == 0:
            raise ValueError("eff_by_ti must be non-empty; a product needs at least one efficacy value.")
        for i, eff in enumerate(eff_by_ti):
            if not (0.0 <= eff <= 1.0):
                raise ValueError(f"eff_by_ti entries must be in [0.0, 1.0], got {eff!r} at index {i}.")
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

    def __repr__(self):
        return (f"Product(name={self.name!r}, category={self.category.name}, "
                f"delivery_mode={self.delivery_mode.name}, cost={self.cost})")
