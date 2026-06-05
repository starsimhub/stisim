from stisim.logistics import ProductCategory
from stisim.logistics.supply import Supply


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
        # NOTE: `supplies` is stored by reference, and the _supplies_by_name/_supplies_by_category indexes below are
        # built once here. Mutating the passed-in list (or self.supplies) after construction would desync the
        # indexes. This is fine given intended usage (supplies are fixed at construction); if mutation is ever
        # needed, add add_supply/remove_supply methods that rebuild the indexes, or copy with list(supplies).
        self.supplies = [] if supplies is None else supplies
        self._supplies_by_name = {supply.product.name: supply for supply in self.supplies}
        if len(self._supplies_by_name) < len(self.supplies):
            raise self.DuplicateSupplyException(f"Two or more products with the same name were added to a Supplies. "
                                                f"Product names must be unique.")
        self._supplies_by_category = {}
        for supply in self.supplies:
            category = supply.product.category
            if category not in self._supplies_by_category:
                self._supplies_by_category[category] = []
            self._supplies_by_category[category].append(supply)

    def use(self, prod_name, quantity):
        """use a specified quantity of a product, returning the remaining quantity and usage cost"""
        supply = self.get_supply(prod_name=prod_name)
        remaining_quantity, cost = supply.use(quantity=quantity)
        return remaining_quantity, cost

    @property
    def accrued_cost(self):
        """the summed cost of all contained Supply objects"""
        cost = sum([supply.accrued_cost for supply in self.supplies])
        return cost

    @property
    def product_names(self):
        return self._supplies_by_name.keys()

    def has_product(self, prod_name):
        return prod_name in self._supplies_by_name

    def get_quantity(self, prod_name):
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
        strs = [f"{supply.product.name}: {supply.quantity}" for supply in self.supplies]
        result = '\n'.join(strs)
        return result
