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

| Parameter | Default | Description |
|-----------|---------|-------------|
| `dur_inf` | lognorm(3 mo, 1 mo) | Duration of infection |
| `beta` | 1.0 | Transmission probability (placeholder) |
| `init_prev` | 0.0 | Initial prevalence |
