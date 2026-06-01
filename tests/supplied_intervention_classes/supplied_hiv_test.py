from typing import Callable

from tests.supplied_intervention_classes.supplied_intervention import SuppliedIntervention
from tests.supplied_intervention_classes.supplies import Supplies


class SuppliedHIVTest(SuppliedIntervention):
    def __init__(self,
                 name: str,
                 eligibilities: list[Callable] = None,
                 coverages: list[float] = None,
                 supplies: Supplies = None, pars=None, *args, **kwargs):
        super().__init__(name=name, supplies=supplies, eligibilities=eligibilities, coverages=coverages, *args, **kwargs)
