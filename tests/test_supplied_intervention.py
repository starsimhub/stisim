"""
Unit tests for SuppliedIntervention class

Tests to ensure appropriate behavior of SuppliedIntervention, with a focus on
accrued cost tracking as Supply objects are consumed via use().
"""
import matplotlib.pyplot as plt
import sys
from pathlib import Path

import pytest
import sciris as sc

from stisim.interventions.product import Product
from stisim.interventions.supply import Supply
from stisim.interventions.supplies import Supplies
from stisim.interventions.supplied_intervention import SuppliedIntervention

tests_directory = Path(__file__).resolve().parent
sys.path.append(str(tests_directory))

verbose = False
do_plot = False
sc.options(interactive=False)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_product(name='oral_prep', cost=10.0):
    return Product(name=name, type='prep', delivery_mode='pill', cost=cost, eff_by_ti=[1.0])


def make_supply(product, quantity=100):
    return Supply(quantity=quantity, product=product)


def make_intervention(*supplies):
    return SuppliedIntervention(supplies=Supplies(list(supplies)))


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@sc.timer()
def test_use_records_cost_for_single_supply():
    sc.heading("Ensuring a single use() call correctly accrues cost on the intervention.")

    product = make_product(cost=10.0)
    supply = make_supply(product, quantity=50)
    intervention = make_intervention(supply)

    quantity_used = 3
    intervention.use(prod_name=product.name, quantity=quantity_used)

    expected_cost = quantity_used * product.cost
    assert intervention.accrued_cost == expected_cost, \
        f"Expected accrued_cost {expected_cost} after using {quantity_used} units at cost {product.cost}, " \
        f"got {intervention.accrued_cost}"


@sc.timer()
def test_use_accumulates_cost_across_multiple_calls():
    sc.heading("Ensuring repeated use() calls accumulate correctly in accrued_cost.")

    product = make_product(cost=5.0)
    supply = make_supply(product, quantity=100)
    intervention = make_intervention(supply)

    intervention.use(prod_name=product.name, quantity=4)
    intervention.use(prod_name=product.name, quantity=6)

    expected_cost = (4 + 6) * product.cost
    assert intervention.accrued_cost == expected_cost, f"Expected accumulated accrued_cost {expected_cost}, got {intervention.accrued_cost}"


@sc.timer()
def test_use_reflects_costs_across_multiple_contained_supplies():
    sc.heading("Ensuring accrued_cost sums costs across multiple Supply objects.")

    product_a = make_product(name='product_a', cost=10.0)
    product_b = make_product(name='product_b', cost=3.0)
    supply_a = make_supply(product_a, quantity=50)
    supply_b = make_supply(product_b, quantity=50)
    intervention = make_intervention(supply_a, supply_b)

    qty_a, qty_b = 2, 5
    intervention.use(prod_name=product_a.name, quantity=qty_a)
    intervention.use(prod_name=product_b.name, quantity=qty_b)

    expected_cost = qty_a * product_a.cost + qty_b * product_b.cost
    assert intervention.accrued_cost == expected_cost, f"Expected combined accrued_cost {expected_cost}, got {intervention.accrued_cost}"


@sc.timer()
def test_use_returns_remaining_quantity_and_cost():
    sc.heading("Ensuring use() returns correct remaining quantity and per-call cost.")

    product = make_product(cost=7.0)
    supply = make_supply(product, quantity=20)
    intervention = make_intervention(supply)

    remaining, cost = intervention.use(prod_name=product.name, quantity=5)

    assert remaining == 15, f"Expected 15 remaining, got {remaining}"
    assert cost == 35.0, f"Expected per-call cost 35.0, got {cost}"


@sc.timer()
def test_use_raises_on_insufficient_supply():
    sc.heading("Ensuring use() raises InsufficientSupplyException when quantity exceeds available supply.")

    product = make_product()
    supply = make_supply(product, quantity=2)
    intervention = make_intervention(supply)

    with pytest.raises(Supply.InsufficientSupplyException):
        intervention.use(prod_name=product.name, quantity=10)


if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    timer = sc.timer()

    test_use_records_cost_for_single_supply()
    test_use_accumulates_cost_across_multiple_calls()
    test_use_reflects_costs_across_multiple_contained_supplies()
    test_use_returns_remaining_quantity_and_cost()
    test_use_raises_on_insufficient_supply()

    sc.heading("Total:")
    timer.toc()

    if do_plot:
        plt.show()
