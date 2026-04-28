"""
Unit tests for Supply class

Tests to ensure appropriate behavior of Supply, including initialization,
quantity updates on use and add, cost accrual, and overuse prevention.
"""

import matplotlib.pyplot as plt
import pytest
import sciris as sc
import sys

from pathlib import Path

from stisim.interventions.product import Product
from stisim.interventions.supply import Supply

tests_directory = Path(__file__).resolve().parent
sys.path.append(str(tests_directory))

verbose = False
do_plot = False
sc.options(interactive=False)

def _make_product(cost=10.0):
    return Product(name='test_product', type='prep', delivery_mode='pill', cost=cost, eff_by_ti=[1.0])


@sc.timer()
def test_supply_initializes_correctly():
    sc.heading("Ensuring Supply stores quantity and initializes accrued_cost to zero.")

    product = _make_product(cost=5.0)
    supply = Supply(quantity=100, product=product)

    assert supply.quantity == 100, f"Expected initial quantity 100, got {supply.quantity}"
    assert supply.accrued_cost == 0, f"Expected initial accrued_cost 0, got {supply.accrued_cost}"
    assert supply.product is product, "Expected supply.product to reference the provided Product"


@sc.timer()
def test_use_decrements_quantity():
    sc.heading("Ensuring use() decrements quantity by the amount used.")

    supply = Supply(quantity=50, product=_make_product())
    supply.use(20)

    assert supply.quantity == 30, f"Expected quantity 30 after using 20 from 50, got {supply.quantity}"


@sc.timer()
def test_use_returns_remaining_and_cost():
    sc.heading("Ensuring use() returns (remaining_quantity, cost_of_usage).")

    product = _make_product(cost=4.0)
    supply = Supply(quantity=10, product=product)
    remaining, cost = supply.use(3)

    assert remaining == 7, f"Expected remaining 7, got {remaining}"
    assert cost == 12.0, f"Expected cost 12.0 (3 * 4.0), got {cost}"


@sc.timer()
def test_use_accrues_cost():
    sc.heading("Ensuring use() adds the usage cost to accrued_cost.")

    product = _make_product(cost=7.5)
    supply = Supply(quantity=100, product=product)
    supply.use(4)

    assert supply.accrued_cost == 30.0, f"Expected accrued_cost 30.0 (4 * 7.5), got {supply.accrued_cost}"


@sc.timer()
def test_multiple_uses_accumulate_cost():
    sc.heading("Ensuring repeated use() calls accumulate cost in accrued_cost.")

    product = _make_product(cost=2.0)
    supply = Supply(quantity=100, product=product)
    supply.use(5)
    supply.use(3)
    supply.use(10)

    expected_cost = (5 + 3 + 10) * 2.0
    assert supply.accrued_cost == expected_cost, \
        f"Expected accumulated cost {expected_cost}, got {supply.accrued_cost}"
    assert supply.quantity == 82, f"Expected quantity 82 after using 18 from 100, got {supply.quantity}"


@sc.timer()
def test_use_entire_supply_does_not_raise():
    sc.heading("Ensuring use() of the entire available quantity does not raise an exception.")

    supply = Supply(quantity=25, product=_make_product())
    remaining, _ = supply.use(25)

    assert remaining == 0, f"Expected remaining 0 after using full supply, got {remaining}"
    assert supply.quantity == 0, f"Expected quantity 0, got {supply.quantity}"


@sc.timer()
def test_use_raises_on_overuse():
    sc.heading("Ensuring use() raises InsufficientSupplyException when requested quantity exceeds supply.")

    supply = Supply(quantity=10, product=_make_product())

    with pytest.raises(Supply.InsufficientSupplyException):
        supply.use(11)


@sc.timer()
def test_use_does_not_modify_state_on_overuse():
    sc.heading("Ensuring use() leaves quantity and accrued_cost unchanged when overuse is attempted.")

    product = _make_product(cost=3.0)
    supply = Supply(quantity=10, product=product)
    supply.use(5)  # legitimate use first: quantity=5, accrued_cost=15

    try:
        supply.use(6)  # overuse attempt
    except Supply.InsufficientSupplyException:
        pass

    assert supply.quantity == 5, \
        f"Expected quantity 5 (unchanged after failed use), got {supply.quantity}"
    assert supply.accrued_cost == 15.0, \
        f"Expected accrued_cost 15.0 (unchanged after failed use), got {supply.accrued_cost}"


@sc.timer()
def test_add_increases_quantity_and_returns_new_quantity():
    sc.heading("Ensuring add() increases the available quantity.")

    supply = Supply(quantity=10, product=_make_product())
    result = supply.add(40)

    assert supply.quantity == 50, f"Expected quantity 50 after adding 40 to 10, got {supply.quantity}"
    assert result == 50, f"Expected return value 50, got {result}"


@sc.timer()
def test_add_does_not_affect_accrued_cost():
    sc.heading("Ensuring add() does not change accrued_cost.")

    product = _make_product(cost=5.0)
    supply = Supply(quantity=10, product=product)
    supply.use(4)  # accrued_cost = 20.0

    supply.add(100)

    assert supply.accrued_cost == 20.0, \
        f"Expected accrued_cost 20.0 (unchanged by add), got {supply.accrued_cost}"


if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    timer = sc.timer()

    test_supply_initializes_correctly()
    test_use_decrements_quantity()
    test_use_returns_remaining_and_cost()
    test_use_accrues_cost()
    test_multiple_uses_accumulate_cost()
    test_use_entire_supply_does_not_raise()
    test_use_raises_on_overuse()
    test_use_does_not_modify_state_on_overuse()
    test_add_increases_quantity_and_returns_new_quantity()
    test_add_does_not_affect_accrued_cost()

    sc.heading("Total:")
    timer.toc()

    if do_plot:
        plt.show()
