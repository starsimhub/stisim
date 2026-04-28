"""
Unit tests for Supplies class

Tests to ensure appropriate behavior of Supplies, including quantity updates on
use, cost accrual, overuse prevention, and retrieval of Supply objects by product
name or type (but not both simultaneously).
"""

import matplotlib.pyplot as plt
import pytest
import sciris as sc
import sys

from pathlib import Path

from stisim.interventions.product import Product
from stisim.interventions.supply import Supply
from stisim.interventions.supplies import Supplies

tests_directory = Path(__file__).resolve().parent
sys.path.append(str(tests_directory))

verbose = False
do_plot = False
sc.options(interactive=False)


def _make_product(name='pill', type='prep', cost=10.0):
    return Product(name=name, type=type, delivery_mode='pill', cost=cost, eff_by_ti=[1.0])


def _make_supply(name='pill', type='prep', quantity=50, cost=10.0):
    return Supply(quantity=quantity, product=_make_product(name=name, type=type, cost=cost))


@sc.timer()
def test_supplies_initializes_empty():
    sc.heading("Ensuring Supplies initializes without error when no supplies are provided.")

    s = Supplies()

    assert s.supplies == [], f"Expected empty list, got {s.supplies}"
    assert s.accrued_cost == 0, f"Expected accrued_cost 0, got {s.accrued_cost}"


@sc.timer()
def test_supplies_initializes_with_supplies():
    sc.heading("Ensuring Supplies stores the provided Supply objects.")

    supply = _make_supply()
    s = Supplies(supplies=[supply])

    assert supply in s.supplies, "Expected provided Supply to be stored in Supplies"


@sc.timer()
def test_supplies_raises_on_duplicate_product_name():
    sc.heading("Ensuring Supplies raises DuplicateSupplyException when two supplies share a product name.")

    supply_a = _make_supply(name='oral_prep')
    supply_b = _make_supply(name='oral_prep')

    with pytest.raises(Supplies.DuplicateSupplyException):
        Supplies(supplies=[supply_a, supply_b])


@sc.timer()
def test_use_decrements_quantity():
    sc.heading("Ensuring use() decrements the quantity of the named product.")

    supply = _make_supply(name='oral_prep', quantity=40)
    s = Supplies(supplies=[supply])
    s.use(prod_name='oral_prep', quantity=15)

    assert s.get_quantity('oral_prep') == 25, \
        f"Expected quantity 25 after using 15 from 40, got {s.get_quantity('oral_prep')}"


@sc.timer()
def test_use_returns_remaining_and_cost():
    sc.heading("Ensuring use() returns (remaining_quantity, cost).")

    supply = _make_supply(name='oral_prep', quantity=20, cost=5.0)
    s = Supplies(supplies=[supply])
    remaining, cost = s.use(prod_name='oral_prep', quantity=4)

    assert remaining == 16, f"Expected remaining 16, got {remaining}"
    assert cost == 20.0, f"Expected cost 20.0 (4 * 5.0), got {cost}"


@sc.timer()
def test_use_raises_on_overuse():
    sc.heading("Ensuring use() raises Supply.InsufficientSupplyException when quantity exceeds supply.")

    supply = _make_supply(name='oral_prep', quantity=10)
    s = Supplies(supplies=[supply])

    with pytest.raises(Supply.InsufficientSupplyException):
        s.use(prod_name='oral_prep', quantity=11)


@sc.timer()
def test_use_does_not_modify_state_on_overuse():
    sc.heading("Ensuring quantity and accrued_cost are unchanged after a failed overuse attempt.")

    supply = _make_supply(name='oral_prep', quantity=10, cost=3.0)
    s = Supplies(supplies=[supply])
    s.use(prod_name='oral_prep', quantity=4)  # quantity=6, accrued_cost=12

    try:
        s.use(prod_name='oral_prep', quantity=7)
    except Supply.InsufficientSupplyException:
        pass

    assert s.get_quantity('oral_prep') == 6, \
        f"Expected quantity 6 (unchanged after failed use), got {s.get_quantity('oral_prep')}"
    assert s.accrued_cost == 12.0, \
        f"Expected accrued_cost 12.0 (unchanged after failed use), got {s.accrued_cost}"


@sc.timer()
def test_accrued_cost_sums_all_supplies():
    sc.heading("Ensuring accrued_cost reflects the sum of costs across all contained Supply objects.")

    supply_a = _make_supply(name='oral_prep', type='prep', quantity=50, cost=2.0)
    supply_b = _make_supply(name='injectable_prep', type='prep', quantity=50, cost=5.0)
    s = Supplies(supplies=[supply_a, supply_b])

    s.use(prod_name='oral_prep', quantity=3)   # cost += 6.0
    s.use(prod_name='injectable_prep', quantity=2)   # cost += 10.0

    assert s.accrued_cost == 16.0, f"Expected total accrued_cost 16.0, got {s.accrued_cost}"


@sc.timer()
def test_get_supply_by_name_returns_correct_supply():
    sc.heading("Ensuring get_supply(prod_name=...) returns the Supply for the specified product name.")

    supply = _make_supply(name='oral_prep')
    s = Supplies(supplies=[supply])
    result = s.get_supply(prod_name='oral_prep')

    assert result is supply, "Expected get_supply to return the exact Supply object for 'oral_prep'"


@sc.timer()
def test_get_supply_by_name_raises_if_missing():
    sc.heading("Ensuring get_supply(prod_name=...) raises MissingSupplyException for an unknown name.")

    s = Supplies(supplies=[_make_supply(name='oral_prep')])

    with pytest.raises(Supplies.MissingSupplyException):
        s.get_supply(prod_name='injectable_prep')


@sc.timer()
def test_get_supply_by_type_returns_correct_supplies():
    sc.heading("Ensuring get_supply(prod_type=...) returns all Supply objects of the specified type.")

    supply = _make_supply(name='oral_prep', type='prep')
    s = Supplies(supplies=[supply])
    result = s.get_supply(prod_type='prep')

    assert isinstance(result, list), f"Expected a list, got {type(result)}"
    assert len(result) == 1, f"Expected exactly one prep Supply to be returned, but was: {len(result)}"
    assert supply in result, "Expected get_supply by type to include the matching Supply"


@sc.timer()
def test_get_supply_by_type_returns_all_matching():
    sc.heading("Ensuring get_supply(prod_type=...) returns every Supply sharing that type.")

    supply_a = _make_supply(name='oral_prep', type='prep')
    supply_b = _make_supply(name='injectable_prep', type='prep')
    supply_c = _make_supply(name='condom', type='barrier')
    s = Supplies(supplies=[supply_a, supply_b, supply_c])

    result = s.get_supply(prod_type='prep')
    assert isinstance(result, list), f"Expected a list, got {type(result)}"
    assert len(result) == 2, f"Expected 2 prep supplies, got {len(result)}"
    assert supply_a in result, "Expected oral_prep in prep supplies"
    assert supply_b in result, "Expected injectable_prep in prep supplies"
    assert supply_c not in result, "Expected condom not in prep supplies"

    result = s.get_supply(prod_type='barrier')
    assert isinstance(result, list), f"Expected a list, got {type(result)}"
    assert len(result) == 1, f"Expected 1 condom supplies, got {len(result)}"
    assert supply_a not in result, "Expected oral_prep not in barrier supplies"
    assert supply_b not in result, "Expected injectable_prep not in barrier supplies"
    assert supply_c in result, "Expected condom in barrier supplies"


@sc.timer()
def test_get_supply_by_type_raises_if_missing():
    sc.heading("Ensuring get_supply(prod_type=...) raises MissingSupplyException for an unknown type.")

    s = Supplies(supplies=[_make_supply(name='oral_prep', type='prep')])

    with pytest.raises(Supplies.MissingSupplyException):
        s.get_supply(prod_type='art')


@sc.timer()
def test_get_supply_raises_if_both_specified():
    sc.heading("Ensuring get_supply() raises ValueError when both prod_name and prod_type are provided.")

    s = Supplies(supplies=[_make_supply(name='oral_prep', type='prep')])

    with pytest.raises(ValueError):
        s.get_supply(prod_name='oral_prep', prod_type='prep')


@sc.timer()
def test_get_supply_raises_if_neither_specified():
    sc.heading("Ensuring get_supply() raises ValueError when neither prod_name nor prod_type is provided.")

    s = Supplies(supplies=[_make_supply()])

    with pytest.raises(ValueError):
        s.get_supply()


if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    timer = sc.timer()

    test_supplies_initializes_empty()
    test_supplies_initializes_with_supplies()
    test_supplies_raises_on_duplicate_product_name()
    test_use_decrements_quantity()
    test_use_returns_remaining_and_cost()
    test_use_raises_on_overuse()
    test_use_does_not_modify_state_on_overuse()
    test_accrued_cost_sums_all_supplies()
    test_get_supply_by_name_returns_correct_supply()
    test_get_supply_by_name_raises_if_missing()
    test_get_supply_by_type_returns_correct_supplies()
    test_get_supply_by_type_returns_all_matching()
    test_get_supply_by_type_raises_if_missing()
    test_get_supply_raises_if_both_specified()
    test_get_supply_raises_if_neither_specified()

    sc.heading("Total:")
    timer.toc()

    if do_plot:
        plt.show()
