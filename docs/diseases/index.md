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

All diseases can be passed to `sti.Sim` by name (string alias) or as module instances. When passed as a string, default parameters are used. To customize, either pass parameters via `sti_pars` or create the module directly:

```python
# By name with defaults
sim = sti.Sim(diseases='hiv')

# By name with custom parameters
sim = sti.Sim(diseases=['hiv', 'ng'], sti_pars=dict(hiv=dict(init_prev=0.1)))

# As module instances
sim = sti.Sim(diseases=[sti.HIV(init_prev=0.1), sti.Gonorrhea()])
```
