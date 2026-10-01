# Coinfection connectors

Connectors model interactions between co-circulating diseases. They run each timestep and adjust agents' relative susceptibility (`rel_sus`) and relative transmissibility (`rel_trans`) on a target disease module, based on their state in another disease.

## Why connectors?

A multi-STI simulation is not just a sum of independent diseases — biological evidence supports several coinfection effects:

- **HIV ↔ ulcerative STIs (syphilis, GUD)**: Ulcers increase HIV acquisition and transmission risk.
- **HIV ↔ non-ulcerative STIs (gonorrhea, chlamydia, trichomoniasis)**: Inflammation increases HIV transmissibility; HIV-induced immunosuppression can prolong STI infection.
- **HIV ↔ BV**: Dysbiosis is associated with elevated HIV acquisition risk.
- **GUD ↔ syphilis**: Syphilitic chancres are a primary cause of GUD; the connector links syphilis stage to GUD prevalence.

Each effect is implemented as a connector class so that users can mix and match.

## Available connectors

| Class | Diseases | Effect |
|-------|----------|--------|
| `sti.hiv_syph` | HIV ↔ syphilis | Bidirectional. Syphilis raises HIV `rel_sus` (default 2.67×) and `rel_trans` (1.2×); HIV/AIDS state can modify syphilis `rel_sus` / `rel_trans` |
| `sti.hiv_tv` | trichomoniasis → HIV | Trichomoniasis infection raises HIV `rel_sus` (default 1.5×) |
| `sti.hiv_ng` | gonorrhea → HIV | Gonorrhea infection raises HIV `rel_sus` (default 1.2×) and `rel_trans` (1.2×) |
| `sti.hiv_ct` | chlamydia → HIV | Chlamydia infection raises HIV `rel_sus` (default 1× — placeholder) |
| `sti.hiv_bv` | BV → HIV | CST-IV state raises HIV `rel_sus` (default 2×) and `rel_trans` (2×) |
| `sti.gud_syph` | syphilis → GUD | Syphilis stage drives GUD prevalence |

Most connectors are unidirectional (the STI affects HIV but not vice versa); `hiv_syph` is the exception.

## Parameters

Multipliers are applied to `rel_sus` or `rel_trans` of the affected disease while the coinfection is present. The Source column follows the conventions of the [parameter reference](../parameters.qmd).

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `rel_sus_hiv_syph` | 2.67 | HIV susceptibility multiplier with symptomatic syphilis (`sti.hiv_syph`) | Not recorded |
| `rel_trans_hiv_syph` | 1.2 | HIV transmissibility multiplier with symptomatic syphilis (`sti.hiv_syph`) | Not recorded |
| `rel_sus_syph_hiv` | 1 | Syphilis susceptibility multiplier with HIV (CD4 < 500) (`sti.hiv_syph`) | Not applicable |
| `rel_sus_syph_aids` | 1 | Syphilis susceptibility multiplier with AIDS (CD4 < 200) (`sti.hiv_syph`) | Not applicable |
| `rel_trans_syph_hiv` | 1 | Syphilis transmissibility multiplier with HIV (CD4 < 500) (`sti.hiv_syph`) | Not applicable |
| `rel_trans_syph_aids` | 1 | Syphilis transmissibility multiplier with AIDS (CD4 < 200) (`sti.hiv_syph`) | Not applicable |
| `rel_sus_hiv_tv` | 1.5 | HIV susceptibility multiplier with trichomoniasis (`sti.hiv_tv`) | [WHO trichomoniasis fact sheet](https://www.who.int/news-room/fact-sheets/detail/trichomoniasis); [Masha 2019](https://doi.org/10.1136/sextrans-2018-053713) |
| `rel_sus_hiv_ng` | 1.2 | HIV susceptibility multiplier with gonorrhea (`sti.hiv_ng`) | Not recorded |
| `rel_trans_hiv_ng` | 1.2 | HIV transmissibility multiplier with gonorrhea (`sti.hiv_ng`) | Not recorded |
| `rel_sus_hiv_ct` | 1 | HIV susceptibility multiplier with chlamydia (`sti.hiv_ct`) | Assumption (placeholder) |
| `rel_sus_hiv_bv` | 2 | HIV susceptibility multiplier with CST IV (`sti.hiv_bv`) | Assumption (placeholder) |
| `rel_trans_hiv_bv` | 2 | HIV transmissibility multiplier with CST IV (`sti.hiv_bv`) | Assumption (placeholder) |
| `rel_sus_syph_gud` | 2 | Syphilis susceptibility multiplier with GUD (`sti.gud_syph`) | Not recorded |
| `rel_trans_syph_gud` | 2 | Syphilis transmissibility multiplier with GUD (`sti.gud_syph`) | Not recorded |
| `rel_sus_gud_syph` | 2 | GUD susceptibility multiplier with syphilis (`sti.gud_syph`) | Not recorded |
| `rel_trans_gud_syph` | 2 | GUD transmissibility multiplier with syphilis (`sti.gud_syph`) | Not recorded |

## Adding connectors

Connectors are passed explicitly when constructing an `sti.Sim`. By default, STIsim automatically wires connectors if `connectors=True` (the default). Alternatively, provide an explicit list of connectors instead:

```python
hiv = sti.HIV()
syph = sti.Syphilis()
sim = sti.Sim(
    diseases=[hiv, syph],
    connectors=[sti.hiv_syph(hiv_module=hiv, syphilis_module=syph, rel_sus_hiv_syph=3.0)],
)
```

## Writing a custom connector

A connector is a subclass of `ss.Connector` that overrides `step()` to update `rel_sus` / `rel_trans` arrays:

```python
class my_connector(ss.Connector):
    def __init__(self, disease_a, disease_b, **kwargs):
        super().__init__()
        self.a = disease_a
        self.b = disease_b
        self.define_pars(rel_sus_a_b=2.0)
        self.update_pars(**kwargs)

    def step(self):
        infected_b = self.b.infected
        self.a.rel_sus[infected_b] = self.pars.rel_sus_a_b
        self.a.rel_sus[~infected_b] = 1.0
```

Connectors run after disease `step()` and before transmission, so changes to `rel_sus` / `rel_trans` take effect on the same timestep.

