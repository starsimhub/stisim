# Analyzers

Analyzers observe the simulation each timestep and store derived results without
changing dynamics. STIsim ships analyzers for coinfection, sex-work transmission,
and network/partnership structure, on top of Starsim's analyzer framework. Attach
them via the `analyzers=[...]` argument to `sti.Sim`.

```python
sim = sti.Sim(
    diseases=['hiv', 'syphilis'],
    analyzers=[sti.coinfection_stats(), sti.PartnershipFormationAnalyzer()],
)
```

## Available analyzers

| Analyzer | What it tracks |
|----------|----------------|
| `sti.coinfection_stats` | Coinfection statistics for two diseases (conditional prevalences). |
| `sti.sw_stats` | New infections and transmissions among sex workers and their clients. |
| `sti.art_coverage` | ART coverage (number and proportion) by sex and age bin. |
| `sti.RelationshipDurations` | Relationship duration distributions in a `StructuredSexual` network. |
| `sti.NetworkDegree` | Lifetime partner-count distributions. |
| `sti.TimeBetweenRelationships` | Gaps between successive relationships. |
| `sti.partner_age_diff` | Age differences between sexual partners. |
| `sti.DebutAge` | Proportion of agents sexually debuted by age. |
| `sti.PartnershipFormationAnalyzer` | Partnership formation per network, gender, and age bin. Relies on networks recording the timestep each relationship ends. |

The grouped-result analyzers build on `result_grouper`, which provides conditional
probability utilities for stratified results.

> **Stub** — expand with usage examples for each analyzer, the result keys they
> write, and how to plot them. See the API reference for
> [`analyzers`](../api/analyzers.qmd) and the
> [Results tutorial](../tutorials/tut_results.qmd) for custom analysis patterns.
