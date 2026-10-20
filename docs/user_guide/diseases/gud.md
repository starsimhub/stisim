# Genital ulcer disease (GUD)

GUD is a simple SIS infection model used as a syndromic proxy for ulcerative STIs. It can serve as a standalone module or as a placeholder for use with connectors (e.g., to model the interaction between GUD and HIV transmission).

**Class:** `sti.GUD` | **Alias:** `'gud'` | **Base class:** `ss.Infection`

## States and transitions

```
     ┌─────────────┐           ┌─────────────┐
     │ Susceptible │──────────▶│  Infected   │
     │             │◀──────────│  (~3 mo)    │
     └─────────────┘ recovery  └─────────────┘
```

## Parameters

Pars of `sti.GUD`. The Source column follows the conventions of the [parameter reference](../../parameters.qmd).

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `dur_inf` | lognorm(3 mo, 1 mo) | Duration of infection | Not recorded |
| `beta_m2f` | 1.0 | Per-act male-to-female transmission probability | Assumption (placeholder) |
| `rel_beta_f2m` | 1.0 | Female-to-male per-act probability relative to `beta_m2f` | Assumption (placeholder) |
| `beta_m2c` | 0 | Not used: no mother-to-child transmission | Not applicable |
| `beta_m2m` | 0 | Per-act male-to-male transmission probability (used on an `msm` network) | Not applicable |
| `eff_condom` | 1.0 | Efficacy of condoms | Not recorded |
| `init_prev` | 0.0 | Initial prevalence | Not applicable |
