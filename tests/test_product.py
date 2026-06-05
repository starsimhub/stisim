"""
Unit tests for Product class

Tests to ensure appropriate behavior of Product, including initialization and
efficacy lookup at valid and out-of-range time indices.
"""
import matplotlib.pyplot as plt
import pytest
import sciris as sc
import sys

from pathlib import Path

from stisim.logistics.product import Product
from stisim.logistics import ProductCategory, DeliveryMode

tests_directory = Path(__file__).resolve().parent
sys.path.append(str(tests_directory))


@sc.timer()
def test_product_attributes_stored():
    sc.heading("Ensuring Product stores all constructor arguments correctly.")

    eff_by_ti = [1.0, 0.8, 0.6, 0.4]
    p = Product(name='lenacapavir', category=ProductCategory.PREP, delivery_mode=DeliveryMode.SHOT, cost=25.0, eff_by_ti=eff_by_ti)

    assert p.name == 'lenacapavir', f"Expected name 'lenacapavir', got '{p.name}'"
    assert p.category is ProductCategory.PREP, f"Expected category ProductCategory.PREP, got {p.category!r}"
    assert p.delivery_mode is DeliveryMode.SHOT, f"Expected delivery_mode DeliveryMode.SHOT, got {p.delivery_mode!r}"
    assert p.cost == 25.0, f"Expected cost 25.0, got {p.cost}"
    assert p.eff_by_ti == eff_by_ti, f"Expected eff_by_ti {eff_by_ti}, got {p.eff_by_ti}"


@sc.timer()
def test_max_durability_matches_eff_by_ti_length():
    sc.heading("Ensuring max_durability equals len(eff_by_ti).")

    eff_by_ti = [1.0, 0.5, 0.0]
    p = Product(name='x', category=ProductCategory.PREP, delivery_mode=DeliveryMode.SHOT, cost=5.0, eff_by_ti=eff_by_ti)

    assert p.max_durability == len(eff_by_ti), \
        f"Expected max_durability {len(eff_by_ti)}, got {p.max_durability}"


@sc.timer()
def test_efficacy_at_ti_contained_index():
    sc.heading("Ensuring efficacy_at_ti returns the correct value for an index within range.")

    eff_by_ti = [1.0, 0.8, 0.6, 0.4]
    p = Product(name='oral_prep', category=ProductCategory.PREP, delivery_mode=DeliveryMode.PILL, cost=10.0, eff_by_ti=eff_by_ti)

    assert p.efficacy_at_ti(0) == 1.0, f"Expected efficacy 1.0 at ti=0, got {p.efficacy_at_ti(0)}"
    assert p.efficacy_at_ti(2) == 0.6, f"Expected efficacy 0.6 at ti=2, got {p.efficacy_at_ti(2)}"
    assert p.efficacy_at_ti(3) == 0.4, f"Expected efficacy 0.4 at ti=3 (last index), got {p.efficacy_at_ti(3)}"


@sc.timer()
def test_efficacy_at_ti_out_of_range_index():
    sc.heading("Ensuring efficacy_at_ti returns 0 for an index beyond product durability.")

    eff_by_ti = [1.0, 0.8, 0.6, 0.4]
    p = Product(name='oral_prep', category=ProductCategory.PREP, delivery_mode=DeliveryMode.PILL, cost=10.0, eff_by_ti=eff_by_ti)

    assert p.efficacy_at_ti(4) == 0, \
        f"Expected efficacy 0 at ti=4 (one past end), got {p.efficacy_at_ti(4)}"
    assert p.efficacy_at_ti(100) == 0, \
        f"Expected efficacy 0 at ti=100 (far past end), got {p.efficacy_at_ti(100)}"


@sc.timer()
def test_enum_members_accepted_directly():
    sc.heading("Ensuring Product accepts and stores ProductCategory/DeliveryMode members.")

    p = Product(name='lenacapavir', category=ProductCategory.PREP, delivery_mode=DeliveryMode.SHOT,
                cost=25.0, eff_by_ti=[1.0])

    assert p.category is ProductCategory.PREP, f"Expected ProductCategory.PREP, got {p.category!r}"
    assert p.delivery_mode is DeliveryMode.SHOT, f"Expected DeliveryMode.SHOT, got {p.delivery_mode!r}"


@sc.timer()
def test_invalid_category_raises():
    sc.heading("Ensuring a non-member category (including a raw string) raises ValueError at construction.")

    with pytest.raises(ValueError):
        Product(name='x', category='prep', delivery_mode=DeliveryMode.PILL, cost=1.0, eff_by_ti=[1.0])


@sc.timer()
def test_invalid_delivery_mode_raises():
    sc.heading("Ensuring a non-member delivery_mode (including a raw string) raises ValueError at construction.")

    with pytest.raises(ValueError):
        Product(name='x', category=ProductCategory.PREP, delivery_mode='carrier_pigeon', cost=1.0, eff_by_ti=[1.0])


@sc.timer()
def test_empty_eff_by_ti_raises():
    sc.heading("Ensuring an empty eff_by_ti raises ValueError at construction.")

    with pytest.raises(ValueError):
        Product(name='x', category=ProductCategory.PREP, delivery_mode=DeliveryMode.PILL, cost=1.0, eff_by_ti=[])


@sc.timer()
def test_out_of_range_efficacy_raises():
    sc.heading("Ensuring eff_by_ti entries outside [0.0, 1.0] raise ValueError at construction.")

    with pytest.raises(ValueError):
        Product(name='x', category=ProductCategory.PREP, delivery_mode=DeliveryMode.PILL, cost=1.0, eff_by_ti=[1.5])
    with pytest.raises(ValueError):
        Product(name='x', category=ProductCategory.PREP, delivery_mode=DeliveryMode.PILL, cost=1.0, eff_by_ti=[0.8, -0.2])


@sc.timer()
def test_negative_cost_warns_but_is_permitted():
    sc.heading("Ensuring a negative cost emits a warning but is still accepted.")

    with pytest.warns(Warning):
        p = Product(name='x', category=ProductCategory.PREP, delivery_mode=DeliveryMode.PILL, cost=-3.0, eff_by_ti=[1.0])
    assert p.cost == -3.0, f"Expected negative cost to be stored, got {p.cost}"


if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    timer = sc.timer()

    test_product_attributes_stored()
    test_max_durability_matches_eff_by_ti_length()
    test_efficacy_at_ti_contained_index()
    test_efficacy_at_ti_out_of_range_index()
    test_enum_members_accepted_directly()
    test_invalid_category_raises()
    test_invalid_delivery_mode_raises()
    test_empty_eff_by_ti_raises()
    test_out_of_range_efficacy_raises()
    test_negative_cost_warns_but_is_permitted()

    sc.heading("Total:")
    timer.toc()

    if do_plot:
        plt.show()