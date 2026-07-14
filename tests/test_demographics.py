"""
Tests for stisim demographics modules.
"""

import pandas as pd
import starsim as ss
import stisim as sti


def _make_sim(rel_migration=1.0, n_migrants=2000, n_agents=2000, start=2000, stop=2010):
    data = pd.DataFrame({'Time': [start, stop], 'Value': [n_migrants, n_migrants]})
    migration = sti.Migration(migration_data=data, pars=dict(rel_migration=rel_migration))
    sim = ss.Sim(n_agents=n_agents, start=start, stop=stop, demographics=[migration],
                 networks=[], diseases=[], verbose=0)
    sim.run()
    return sim


def test_rel_migration_scales_immigration():
    """ rel_migration should scale up net immigration proportionally """
    sim1 = _make_sim(rel_migration=1.0, n_migrants=2000)
    sim2 = _make_sim(rel_migration=2.0, n_migrants=2000)
    total1 = sim1.demographics.migration.results.new_migrants.sum()
    total2 = sim2.demographics.migration.results.new_migrants.sum()
    assert total1 > 0
    assert total2 == 2 * total1


def test_rel_migration_scales_emigration():
    """ rel_migration should scale up net emigration (negative values) proportionally """
    sim1 = _make_sim(rel_migration=1.0, n_migrants=-500)
    sim2 = _make_sim(rel_migration=2.0, n_migrants=-500)
    total1 = sim1.demographics.migration.results.new_migrants.sum()
    total2 = sim2.demographics.migration.results.new_migrants.sum()
    assert total1 < 0
    assert total2 == 2 * total1


def test_rel_migration_zero_disables_migration():
    """ rel_migration=0 should produce no migrants """
    sim = _make_sim(rel_migration=0.0, n_migrants=2000)
    assert sim.demographics.migration.results.new_migrants.sum() == 0


if __name__ == '__main__':
    test_rel_migration_scales_immigration()
    test_rel_migration_scales_emigration()
    test_rel_migration_zero_disables_migration()
