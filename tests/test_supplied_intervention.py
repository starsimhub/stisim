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
import starsim as ss

from stisim.logistics.product import Product
from stisim.logistics import ProductCategory, DeliveryMode
from stisim.logistics.supply import Supply
from stisim.logistics.supplies import Supplies
from stisim.logistics.supplied_intervention import SuppliedIntervention

tests_directory = Path(__file__).resolve().parent
sys.path.append(str(tests_directory))

verbose = False
do_plot = False
sc.options(interactive=False)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

test_eligibilities = [lambda sim: sim.people.female == True]

class ConcreteSuppliedIntervention(SuppliedIntervention):
    """Minimal concrete subclass so the abstract base's behavior (use(), cost accrual) can be exercised."""
    def step(self):
        pass


class UsingSuppliedIntervention(SuppliedIntervention):
    """Concrete subclass that uses fixed quantities of named products each step (for in-sim tests).

    `usage` maps product name -> quantity used per step. `n_steps_run` counts how many times step() ran, so a
    test can compute expected remaining quantities independently of the Supplies accounting under test.
    """
    def __init__(self, usage, **kwargs):
        super().__init__(**kwargs)
        self.usage = dict(usage)
        self.n_steps_run = 0

    def step(self):
        for prod_name, quantity in self.usage.items():
            self.use(prod_name=prod_name, quantity=quantity)
        self.n_steps_run += 1


def make_product(name='oral_prep', cost=10.0):
    return Product(name=name, category=ProductCategory.PREP, delivery_mode=DeliveryMode.PILL, cost=cost, eff_by_ti=[1.0])


def make_supply(product, quantity=100):
    return Supply(quantity=quantity, product=product)


def make_intervention(*supplies):
    return ConcreteSuppliedIntervention(supplies=Supplies(list(supplies)), name='intervention_name',
                                        eligibilities=test_eligibilities)


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


@sc.timer()
def test_abstract_base_cannot_be_instantiated():
    sc.heading("Ensuring the abstract SuppliedIntervention base cannot be instantiated without a step() override.")

    with pytest.raises(TypeError):
        SuppliedIntervention(name='abstract', eligibilities=test_eligibilities)


@sc.timer()
def test_out_of_range_coverage_raises():
    sc.heading("Ensuring a coverage outside [0.0, 1.0] raises ValueError at construction.")

    with pytest.raises(ValueError):
        ConcreteSuppliedIntervention(name='x', eligibilities=test_eligibilities, coverages=[1.5])
    with pytest.raises(ValueError):
        ConcreteSuppliedIntervention(name='x', eligibilities=test_eligibilities, coverages=[-0.3])


@sc.timer()
def test_calc_supply_distribution_raises_on_length_mismatch():
    sc.heading("Ensuring calc_supply_distribution raises ValueError when the parallel lists differ in length.")

    intervention = make_intervention()

    with pytest.raises(ValueError):
        intervention.calc_supply_distribution(
            offer_pools=[[], []],
            cur_coverages=[0.0],            # length 1: mismatched with the others (length 2)
            target_coverages=[1.0, 1.0],
            n_eligibles=[10, 10],
            n_supply=5,
        )


@sc.timer()
def test_calc_supply_distribution_raises_on_empty_lists():
    sc.heading("Ensuring calc_supply_distribution raises ValueError when the (aligned) lists are empty.")

    intervention = make_intervention()

    with pytest.raises(ValueError):
        intervention.calc_supply_distribution(
            offer_pools=[],
            cur_coverages=[],
            target_coverages=[],
            n_eligibles=[],
            n_supply=5,
        )


@sc.timer()
def test_calc_supply_distribution_accepts_aligned_lists():
    sc.heading("Ensuring calc_supply_distribution returns a per-group result for correctly aligned lists.")

    intervention = make_intervention()

    result = intervention.calc_supply_distribution(
        offer_pools=[[], []],
        cur_coverages=[0.0, 0.0],
        target_coverages=[1.0, 1.0],
        n_eligibles=[10, 10],
        n_supply=5,
    )

    assert len(result) == 2, f"Expected one distribution count per group (2), got {len(result)}"


@sc.timer()
def test_shared_supplies_accounting_after_sim_run():
    sc.heading("Ensuring a Supplies shared by two interventions accounts for remaining quantities after a sim run.")

    # Three finite supplies. Both interventions use A; intervention 1 also uses B; intervention 2 also uses C.
    A0, B0, C0 = 100_000, 100_000, 100_000
    usage1 = {'A': 2, 'B': 3}   # intervention 1: shared A + exclusive B
    usage2 = {'A': 5, 'C': 11}   # intervention 2: shared A + exclusive C

    shared = Supplies([
        make_supply(make_product('A'), quantity=A0),
        make_supply(make_product('B'), quantity=B0),
        make_supply(make_product('C'), quantity=C0),
    ])
    i1 = UsingSuppliedIntervention(usage=usage1, name='iv1', eligibilities=test_eligibilities, supplies=shared)
    i2 = UsingSuppliedIntervention(usage=usage2, name='iv2', eligibilities=test_eligibilities, supplies=shared)

    sim = ss.Sim(n_agents=50, dur=5, dt=1, interventions=[i1, i2])
    sim.run()

    # ss.Sim deep-copies its inputs, so the runtime objects are NOT i1/i2; read them back from the sim. The copy
    # preserves the shared-by-reference Supplies, which this test guards.
    ri1, ri2 = sim.interventions[0], sim.interventions[1]
    assert ri1.supplies is ri2.supplies, "Expected the shared Supplies reference to survive sim initialization"

    s = ri1.supplies
    n1, n2 = ri1.n_steps_run, ri2.n_steps_run
    assert n1 > 0 and n2 > 0, f"Expected both interventions to have stepped, got n1={n1}, n2={n2}"

    expected_A = A0 - (usage1['A'] * n1 + usage2['A'] * n2)  # A drawn down by BOTH interventions
    expected_B = B0 - (usage1['B'] * n1)                     # B drawn down by intervention 1 only
    expected_C = C0 - (usage2['C'] * n2)                     # C drawn down by intervention 2 only

    assert s.get_quantity('A') == expected_A, \
        f"Shared product A: expected {expected_A} remaining, got {s.get_quantity('A')}"
    assert s.get_quantity('B') == expected_B, \
        f"Product B (intervention 1 only): expected {expected_B} remaining, got {s.get_quantity('B')}"
    assert s.get_quantity('C') == expected_C, \
        f"Product C (intervention 2 only): expected {expected_C} remaining, got {s.get_quantity('C')}"
    return sim

if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    timer = sc.timer()

    test_use_records_cost_for_single_supply()
    test_use_accumulates_cost_across_multiple_calls()
    test_use_reflects_costs_across_multiple_contained_supplies()
    test_use_returns_remaining_quantity_and_cost()
    test_use_raises_on_insufficient_supply()
    test_abstract_base_cannot_be_instantiated()
    test_out_of_range_coverage_raises()
    test_calc_supply_distribution_raises_on_length_mismatch()
    test_calc_supply_distribution_raises_on_empty_lists()
    test_calc_supply_distribution_accepts_aligned_lists()
    test_shared_supplies_accounting_after_sim_run()

    sc.heading("Total:")
    timer.toc()

    if do_plot:
        plt.show()
