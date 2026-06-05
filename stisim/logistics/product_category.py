from enum import Enum


class ProductCategory(Enum):
    """
    Category of a Product, used to determine if 1+ Products have a similar function (for example,
    ProductCategory.PREP applies to both oral and injectable PrEP).

    Pass members directly (ProductCategory.PREP); raw string values are not accepted.
    """
    PREP = 'prep'
    ART = 'art'
    BARRIER = 'barrier'
