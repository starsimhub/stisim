"""Define Supplies: an immutable, name-keyed collection of Supply objects."""
from stisim.logistics import ProductCategory
from stisim.logistics.supply import Supply

__all__ = ["Supplies"]


class Supplies:
    class MissingSupplyException(ValueError):
        pass

    class DuplicateSupplyException(ValueError):
        pass

    def __init__(self, supplies: list[Supply] = None):
        """
        Represents a cache of 1+ Supply of Products. The mix of Supply/Products is arbitrary, for example, the provided
        Supply objects could contain oral PrEP, injectable PrEP, condoms, needles, and appointment slots.

        Args:
            supplies: The individual Supply objects containing Products that will be distributed within an
                intervention.
        """
        # A Supplies is immutable after construction: the contained Supply objects are fixed here and the
        # _supplies_by_name/_supplies_by_category indexes are built once, so they can never desync. We copy the
        # passed-in list so later mutation of the caller's list does not affect us, and expose the contents only
        # through the read-only `supplies` property. (The Supply objects themselves still mutate as they are used.)
        self._supplies = [] if supplies is None else list(supplies)
        self._supplies_by_name = {supply.product.name: supply for supply in self._supplies}
        if len(self._supplies_by_name) < len(self._supplies):
            raise self.DuplicateSupplyException(f"Two or more products with the same name were added to a Supplies. "
                                                f"Product names must be unique.")
        self._supplies_by_category = {}
        for supply in self._supplies:
            category = supply.product.category
            if category not in self._supplies_by_category:
                self._supplies_by_category[category] = []
            self._supplies_by_category[category].append(supply)

    @property
    def supplies(self) -> list[Supply]:
        """The contained Supply objects (a copy; Supplies is immutable after construction)."""
        return list(self._supplies)

    def __len__(self) -> int:
        return len(self._supplies)

    def __iter__(self):
        """Iterate over the contained product names (Supplies is keyed by product name)."""
        return iter(self._supplies_by_name)

    def __contains__(self, prod_name) -> bool:
        return prod_name in self._supplies_by_name

    def __getitem__(self, prod_name) -> Supply:
        """supplies[prod_name] -> the Supply for that product name (raises MissingSupplyException if absent)."""
        return self.get_supply(prod_name)

    def use(self, prod_name, quantity):
        """use a specified quantity of a product, returning the remaining quantity and usage cost"""
        supply = self.get_supply(prod_name=prod_name)
        remaining_quantity, cost = supply.use(quantity=quantity)
        return remaining_quantity, cost

    @property
    def accrued_cost(self):
        """
        The summed cost of all contained Supply objects.

        A Supplies object may be shared across more than one SuppliedIntervention (it can be passed to several at
        construction), so this is the total cost of the pool across every intervention that has drawn on it. For the
        cost attributable to a single intervention, see SuppliedIntervention.accrued_cost instead.
        """
        cost = sum([supply.accrued_cost for supply in self._supplies])
        return cost

    @property
    def product_names(self) -> list[str]:
        return list(self._supplies_by_name)

    def has_product(self, prod_name) -> bool:
        return prod_name in self

    def get_quantity(self, prod_name) -> int:
        return self.get_supply(prod_name=prod_name).quantity

    def get_product(self, prod_name):
        return self.get_supply(prod_name=prod_name).product

    def get_supply(self, prod_name: str) -> Supply:
        """obtain the single contained Supply object for the given product name"""
        requested = self._supplies_by_name.get(prod_name, None)
        if requested is None:
            raise self.MissingSupplyException(f"Supply of product of name: {prod_name} requested but is not included.")
        return requested

    def get_supplies(self, prod_category: ProductCategory) -> list[Supply]:
        """obtain the list of contained Supply objects sharing the given product category"""
        if not isinstance(prod_category, ProductCategory):
            raise ValueError(f"prod_category must be a ProductCategory member, got {prod_category!r}.")
        requested = self._supplies_by_category.get(prod_category, None)
        if requested is None:
            raise self.MissingSupplyException(f"Supplies of category: {prod_category} requested but are not included.")
        return requested

    def __repr__(self):
        strs = [f"{supply.product.name}: {supply.quantity}" for supply in self._supplies]
        result = '\n'.join(strs)
        return result
