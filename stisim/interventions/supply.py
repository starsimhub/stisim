from stisim.interventions.product import Product


class Supply:
    class InsufficientSupplyException(ValueError):
        pass

    def __init__(self, quantity: int, product: Product):
        """
        Represents a distinct quantity of a Product. For example, two Supply objects of the same Product can coexist
        if they are acquired from different sources.

        Args:
            quantity: (int) How much of the product is available in this supply.
            product: (Product) The product for distribution in intervention(s).
        """
        self.product = product
        self.quantity = quantity
        self.accrued_cost = 0

    def use(self, quantity):
        """uses a specified quantity of the product, returning the quantity remaining and cost of usage"""
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
