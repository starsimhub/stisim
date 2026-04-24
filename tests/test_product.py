"""
Unit tests for Product class

Tests to ensure appropriate behavior of Product, including initialization and
efficacy lookup at valid and out-of-range time indices.
"""

import pytest
import sciris as sc

from stisim.interventions.product import Product


@sc.timer()
def test_product_attributes_stored():
    sc.heading("Ensuring Product stores all constructor arguments correctly.")

    eff_by_ti = [1.0, 0.8, 0.6, 0.4]
    p = Product(name='lenacapavir', type='prep', delivery_mode='shot', cost=25.0, eff_by_ti=eff_by_ti)

    assert p.name == 'lenacapavir', f"Expected name 'lenacapavir', got '{p.name}'"
    assert p.type == 'prep', f"Expected type 'prep', got '{p.type}'"
    assert p.delivery_mode == 'shot', f"Expected delivery_mode 'shot', got '{p.delivery_mode}'"
    assert p.cost == 25.0, f"Expected cost 25.0, got {p.cost}"
    assert p.eff_by_ti == eff_by_ti, f"Expected eff_by_ti {eff_by_ti}, got {p.eff_by_ti}"


@sc.timer()
def test_product_ids_are_unique():
    sc.heading("Ensuring each Product instance receives a unique id.")

    p1 = Product(name='a', type='prep', delivery_mode='pill', cost=1.0, eff_by_ti=[1.0])
    p2 = Product(name='a', type='prep', delivery_mode='pill', cost=1.0, eff_by_ti=[1.0])

    assert p1.id != p2.id, f"Expected unique ids for separate Product instances, but both are {p1.id}"


@sc.timer()
def test_max_durability_matches_eff_by_ti_length():
    sc.heading("Ensuring max_durability equals len(eff_by_ti).")

    eff_by_ti = [1.0, 0.5, 0.0]
    p = Product(name='x', type='prep', delivery_mode='shot', cost=5.0, eff_by_ti=eff_by_ti)

    assert p.max_durability == len(eff_by_ti), \
        f"Expected max_durability {len(eff_by_ti)}, got {p.max_durability}"


@sc.timer()
def test_efficacy_at_ti_contained_index():
    sc.heading("Ensuring efficacy_at_ti returns the correct value for an index within range.")

    eff_by_ti = [1.0, 0.8, 0.6, 0.4]
    p = Product(name='oral_prep', type='prep', delivery_mode='pill', cost=10.0, eff_by_ti=eff_by_ti)

    assert p.efficacy_at_ti(0) == 1.0, f"Expected efficacy 1.0 at ti=0, got {p.efficacy_at_ti(0)}"
    assert p.efficacy_at_ti(2) == 0.6, f"Expected efficacy 0.6 at ti=2, got {p.efficacy_at_ti(2)}"
    assert p.efficacy_at_ti(3) == 0.4, f"Expected efficacy 0.4 at ti=3 (last index), got {p.efficacy_at_ti(3)}"


@sc.timer()
def test_efficacy_at_ti_out_of_range_index():
    sc.heading("Ensuring efficacy_at_ti returns 0 for an index beyond product durability.")

    eff_by_ti = [1.0, 0.8, 0.6, 0.4]
    p = Product(name='oral_prep', type='prep', delivery_mode='pill', cost=10.0, eff_by_ti=eff_by_ti)

    assert p.efficacy_at_ti(4) == 0, \
        f"Expected efficacy 0 at ti=4 (one past end), got {p.efficacy_at_ti(4)}"
    assert p.efficacy_at_ti(100) == 0, \
        f"Expected efficacy 0 at ti=100 (far past end), got {p.efficacy_at_ti(100)}"


if __name__ == '__main__':
    sc.heading("Total:")
    timer = sc.timer()

    test_product_attributes_stored()
    test_product_ids_are_unique()
    test_max_durability_matches_eff_by_ti_length()
    test_efficacy_at_ti_contained_index()
    test_efficacy_at_ti_out_of_range_index()

    timer.toc()
