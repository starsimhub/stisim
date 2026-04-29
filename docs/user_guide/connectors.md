# Coinfection connectors

Connectors model interactions between co-circulating diseases. They run each timestep and adjust agents' relative susceptibility (`rel_sus`) and relative transmissibility (`rel_trans`) on a target disease module, based on their state in another disease.

## Why connectors?

A multi-STI simulation is not just a sum of independent diseases — biological evidence supports several coinfection effects:

- **HIV ↔ ulcerative STIs (syphilis, GUD)**: ulcers increase HIV acquisition and transmission risk.
- **HIV ↔ non-ulcerative STIs (gonorrhea, chlamydia, trichomoniasis)**: inflammation increases HIV transmissibility; HIV-induced immunosuppression can prolong STI infection.
- **HIV ↔ BV**: dysbiosis is associated with elevated HIV acquisition risk.
- **GUD ↔ syphilis**: syphilitic chancres are a primary cause of GUD; the connector links syphilis stage to GUD prevalence.

Each effect is implemented as a connector class so that users can mix and match.

## Available connectors

| Class | Diseases | Effect |
|-------|----------|--------|
| `sti.hiv_syph` | HIV ↔ syphilis | Syphilis raises HIV `rel_sus` (default 2.67×) and `rel_trans` (1.2×); HIV/AIDS can modify syphilis acquisition and transmission |
| `sti.hiv_tv` | HIV ↔ trichomoniasis | Symmetric coinfection multipliers |
| `sti.hiv_ng` | HIV ↔ gonorrhea | Symmetric coinfection multipliers |
| `sti.hiv_ct` | HIV ↔ chlamydia | Symmetric coinfection multipliers |
| `sti.hiv_bv` | HIV ↔ bacterial vaginosis | BV state modifies HIV susceptibility |
| `sti.gud_syph` | GUD ↔ syphilis | Syphilis stage drives GUD prevalence |

## Default connectors

If you pass `connectors=True` (or omit it) when creating an `sti.Sim`, the appropriate connectors are added automatically based on which diseases are in the simulation. To suppress this behavior, pass `connectors=False` or supply your own list:

```python
sim = sti.Sim(
    diseases=['hiv', 'syphilis'],
    connectors=[sti.hiv_syph(hiv_module=..., syphilis_module=..., rel_sus_hiv_syph=3.0)],
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

See [the source](https://github.com/starsimhub/stisim/tree/main/stisim/connectors) for full implementations of the built-in connectors.
