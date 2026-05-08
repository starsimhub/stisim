"""
Minimal unit tests for sti.route_pars — flat-kwarg routing into per-category
buckets (sim, sti, nw, dem, connector).
"""
import pytest
import stisim as sti


def test_basic_routing():
    """Flat kwargs route to the correct categories."""
    routed = sti.route_pars(start=2010, debut_f=18, init_prev=0.05, verbose=False)
    assert routed.sim['start'] == 2010
    assert routed.nw['debut_f'] == 18
    assert routed.sti['init_prev'] == 0.05
    return routed


def test_broadcast():
    """A key valid in 2+ categories is broadcast to all of them (e.g. dt)."""
    routed = sti.route_pars(dt=0.5, verbose=False)
    assert routed.sim['dt'] == 0.5
    assert routed.sti['dt'] == 0.5
    assert routed.connector['dt'] == 0.5
    return routed


def test_strict_raises_on_unknown():
    """strict=True raises on unknown pars; non-strict stashes them."""
    with pytest.raises(ValueError):
        sti.route_pars(asdf=1, strict=True, verbose=False)
    routed = sti.route_pars(asdf=1, strict=False, verbose=False)
    assert routed.unmatched == {'asdf': 1}
    return routed


def test_pre_categorized():
    """Pre-categorized dicts merge into their bucket (incl. connector_pars)."""
    routed = sti.route_pars(
        sim_pars=dict(start=2010),
        nw_pars=dict(debut_f=18),
        connector_pars=dict(rel_sus_hiv_syph=3.5),
        verbose=False,
    )
    assert routed.sim['start'] == 2010
    assert routed.nw['debut_f'] == 18
    assert routed.connector['rel_sus_hiv_syph'] == 3.5
    return routed


def test_connector_bucket():
    """Flat connector pars route into the connector bucket."""
    routed = sti.route_pars(rel_sus_hiv_syph=3.5, verbose=False)
    assert routed.connector['rel_sus_hiv_syph'] == 3.5
    return routed


if __name__ == '__main__':
    test_basic_routing()
    test_broadcast()
    test_strict_raises_on_unknown()
    test_pre_categorized()
    test_connector_bucket()
