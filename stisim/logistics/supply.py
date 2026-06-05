import numpy as np

from stisim.logistics.product import Product


class Supply:
    class InsufficientSupplyException(ValueError):
        pass

    def __init__(self, quantity: int, product: Product):
        """
        Represents a distinct quantity of a Product. For example, two Supply objects of the same Product can coexist
        if they are acquired from different sources.

        Args:
            quantity: (int) How much of the product is available in this supply. Must be a finite integer count, or
                np.inf to model an effectively unlimited/infinite supply. Passing np.inf means use() will never raise
                InsufficientSupplyException and the quantity remains np.inf, while usage cost still accrues normally.
                Arbitrary (non-integer) floats are not permitted.
            product: (Product) The product for distribution in intervention(s).
        """
        # quantity is a discrete count of items, so only an int is valid; np.inf is the one allowed exception, used
        # as the sentinel for an unlimited supply. Reject everything else (e.g. 3.5, 2.0) up front.
        if not (isinstance(quantity, int) or quantity == np.inf):
            raise ValueError(f"Supply quantity must be an int (finite count) or np.inf (unlimited), "
                             f"got {quantity!r}.")
        self.product = product
        self.quantity = quantity
        self.accrued_cost = 0

    def use(self, quantity):
        """uses a specified quantity of the product, returning the remaining quantity and cost of usage"""
        remaining = self.quantity - quantity
        if remaining < 0:
            raise self.InsufficientSupplyException(f"Attempted to use {quantity} {self.product.name}. "
                                                   f"Only {self.quantity} exists.")
        self.quantity = remaining

        cost = quantity * self.product.cost
        self.accrued_cost += cost

        return self.quantity, cost

    def add(self, quantity):
        """increase the available quantity of the product"""
        self.quantity += quantity
        return self.quantity
