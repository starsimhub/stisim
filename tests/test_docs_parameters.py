"""
Check that the user-guide parameter tables match the code, and that the parameter table can be built from them.
"""
import sys
import sciris as sc
import stisim as sti

docs = sti.root / 'docs'
sys.path.insert(0, str(sti.root / 'stisim/ai/plugin/skills/parameter-table/scripts'))
import parameter_table

# Internal pars that are placeholders or set at runtime, so are not documented
internal = {'beta', 'dt', 'include_care', 'dur_post_art', 'p_hiv_death', 'birth_outcome_keys', 'p_spontaneous',
            '_f_partnership_taper_offset', 'age_diffs', 'dur_dist', 'match_dist', 'p_actually_looking', 'p_hi_risk', 'p_lo_risk',
            'init_latent_prev'} # A constructor argument of sti.Syphilis (spread over latency), not a parameter


def test_parameter_docs_coverage():
    sc.heading('Every module par is in the user-guide tables, and every documented par exists')
    rows = parameter_table.parse(docs=docs)
    for cls, (undocumented, unknown) in parameter_table.check_coverage(rows, exclude=internal).items():
        undocumented = [k for k in undocumented if not k.endswith('_dist')]
        unknown = [k for k in unknown if k not in internal]
        assert not undocumented, f'{cls} pars missing from the docs: {undocumented}'
        assert not unknown, f'{cls} pars documented but not in the code: {unknown}'
    return rows


def test_parameter_table_sim():
    sc.heading('Parameter table from a sim marks values that differ from the defaults')
    sim = sti.Sim(diseases=['hiv', 'syph'], n_agents=200, start=2000, stop=2001, hiv=dict(beta_m2f=0.02))
    rows = parameter_table.add_sim(parameter_table.parse(docs=docs), sim, label='Test')
    vals = {(r['modname'], r['name']): r['values']['Test'] for r in rows}
    assert vals[('hiv', 'beta_m2f')] == '0.02*'
    assert vals[('hiv', 'eff_condom')] == '0.9'
    assert vals[('ct', 'beta_m2f')] == 'not in sim'
    md = parameter_table.to_markdown(rows)
    assert '| Test |' in md.split('\n')[2]
    return md


if __name__ == '__main__':
    test_parameter_docs_coverage()
    test_parameter_table_sim()
