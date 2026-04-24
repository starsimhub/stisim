# Coinfection connectors

Connectors model interactions between co-circulating diseases by modifying susceptibility and transmissibility when an agent is coinfected.

| File | Connectors |
|------|-----------|
| `hiv_sti.py` | `hiv_syph`, `hiv_tv`, `hiv_ng`, `hiv_ct`, `hiv_bv` -- HIV coinfection effects on each STI and vice versa |
| `gud_syph.py` | `gud_syph` -- Genital ulcer disease interaction with syphilis |

Each connector extends `ss.Connector` and overrides `step()` to adjust `rel_sus` and `rel_trans` arrays on the relevant disease modules each timestep.
