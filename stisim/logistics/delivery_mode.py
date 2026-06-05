from enum import Enum

__all__ = ["DeliveryMode"]


class DeliveryMode(Enum):
    """
    How a Product is delivered, used to distinguish Products that may have the same category but are delivered
    differently in intervention-relevant ways (for example, pill vs. shot).

    Pass members directly (DeliveryMode.PILL); raw string values are not accepted.
    """
    SHOT = 'shot'
    PILL = 'pill'
    TOPICAL = 'topical'
    EDUCATIONAL = 'educational'
