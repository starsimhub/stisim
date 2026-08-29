# Diseases

STIsim includes the following disease modules:

| Disease | Class | Alias | Model type |
|---------|-------|-------|------------|
| [HIV](hiv.md) | `sti.HIV` | `'hiv'` | CD4-based progression |
| [Syphilis](syphilis.md) | `sti.Syphilis` | `'syphilis'` / `'syph'` | Staged (primary/secondary/latent/tertiary) |
| [Chlamydia](chlamydia.md) | `sti.Chlamydia` | `'ct'` | SEIS |
| [Gonorrhea](gonorrhea.md) | `sti.Gonorrhea` | `'ng'` | SEIS |
| [Trichomoniasis](trichomoniasis.md) | `sti.Trichomoniasis` | `'tv'` | SEIS with persistence |
| [Bacterial vaginosis](bv.md) | `sti.BV` | `'bv'` | CST-based microbiome |
| [Genital ulcer disease](gud.md) | `sti.GUD` | `'gud'` | Simple SIS |

## Infected vs. infectious

The SEIS diseases (chlamydia, gonorrhea, trichomoniasis) follow starsim's `ss.SEIR` naming: `exposed` and `infectious` are the literal E and I compartments, and `infected` is derived as E plus I. So `n_infected` and `prevalence` count everyone carrying the infection, including those still in the latent period, while transmission depends on `infectious` alone. Use `n_infectious` if you want the transmitting compartment on its own.

Correspondingly, `ti_exposed` is the time of acquisition and `ti_infectious` the time of becoming infectious; there is no `ti_infected` on an SEIS disease. Incidence (`new_infections`, `incidence`) is counted at acquisition, off `ti_exposed`.

HIV, syphilis and BV set `infected` at acquisition already, so the same reading applies to them: `infected` means "has the infection", not "is transmitting".

All diseases can be passed to `sti.Sim` by name (string alias) or as module instances. When passed as a string, default parameters are used. To customize, either pass parameters via `sti_pars` or create the module directly:

```python
# By name with defaults
sim = sti.Sim(diseases='hiv')

# By name with custom parameters
sim = sti.Sim(diseases=['hiv', 'ng'], sti_pars=dict(hiv=dict(init_prev=0.1)))

# As module instances
sim = sti.Sim(diseases=[sti.HIV(init_prev=0.1), sti.Gonorrhea()])
```
