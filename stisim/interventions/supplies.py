from typing import Iterable

from stisim.interventions.supply import Supply


class Supplies:
    class MissingSupplyException(ValueError):
        pass

    class DuplicateSupplyException(ValueError):
        pass

    def __init__(self, supplies: Iterable[Supply] = None):
        """
        Represents a cache of 1+ Supply of Products. The mix of Supply/Products is arbitrary, for example, the provided
        Supply objects could contain oral PrEP, injectable PrEP, condoms, needles, and appointment slots.

        Args:
            supplies: The individual Supply objects containing Products that will be distributed within an
                intervention.
        """
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

    def get_supply(self, prod_name=None, prod_type=None):
        """obtain a contained Supply object either by name or type"""
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
            requested = self._supplies_by_type.get(prod_type, None)
            if requested is None:
                raise self.MissingSupplyException(f"Supplies of product of type: {prod_type} requested but are "
                                                  f"not included.")
        return requested

    def __repr__(self):
        strs = [f"{supply.product.name}: {supply.quantity}" for supply in self.supplies]
        result = '\n'.join(strs)
        return result
