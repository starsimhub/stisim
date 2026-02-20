# Structured sexual network

**Class:** `sti.StructuredSexual` | **Base class:** `ss.SexualNetwork`

The `StructuredSexual` network is STIsim's primary sexual contact network. It is created automatically when you make a sim.

## Overview

Agents are assigned a **risk group** (low, medium, high) at initialization. Risk group determines partnership patterns: concurrency, relationship type, and duration. Partnerships form each timestep based on age preferences, dissolve after their sampled duration, and are classified as stable, casual, one-time, or sex work.

## Risk groups

| Risk group | Description | Female share | Male share |
|-----------|-------------|-------------|------------|
| 0 (low) | Marry and remain monogamous | 85% | 80% |
| 1 (medium) | Marry then divorce, or have concurrent partners | 14% | 18% |
| 2 (high) | Never marry | 1% | 2% |

## Partnership types

How a partnership is classified depends on the risk groups of both partners:

- **Stable**: Both partners in the same risk group and pass a stability check (RG0: 90%, RG1: 50%, RG2: 0%)
- **Casual**: Mismatched risk groups, 50% chance (otherwise one-time)
- **One-time**: Single-timestep contact
- **Sex work**: FSW-client contacts, one-time duration

## Concurrency

The number of concurrent partners an agent seeks depends on their sex and risk group:

| | Low (RG0) | Medium (RG1) | High (RG2) |
|---|-----------|-------------|------------|
| Female | 0.0001 | 0.01 | 0.1 |
| Male | 0.0001 | 0.2 | 0.5 |

Values represent the parameter of a Poisson distribution (+1), so a male in RG2 typically seeks 1-3 concurrent partners.

## Age mixing

Women sample a preferred partner age from a normal distribution. Parameters vary by female age group and risk group:

| Female age group | RG0 (mean, std) | RG1 | RG2 |
|-----------------|-----------------|-----|-----|
| Teens (<20) | (7, 3) | (6, 3) | (5, 1) |
| Young (20-25) | (8, 3) | (7, 3) | (5, 2) |
| Adult (25+) | (8, 3) | (7, 3) | (5, 2) |

Values are the age difference (male minus female) in years. Partners are matched by sorting both groups by age.

## Sex work

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fsw_shares` | 5% | Proportion of women who are sex workers |
| `client_shares` | 12% | Proportion of men who are clients |
| `sw_seeking_rate` | 1/month | Rate at which clients seek FSWs |

## Coital acts and condoms

| Parameter | Default | Description |
|-----------|---------|-------------|
| `acts` | lognorm(80/yr, 30/yr) | Coital acts per partnership per year |
| `condom_data` | None | Condom use probability (can be stratified by risk group) |

Transmission probability per partnership accounts for both condom use and number of acts:

P(transmission) = 1 - [1 - beta * (1 - eff_condom)]^(acts * p_condom) * [1 - beta]^(acts * (1 - p_condom))

## Other networks

STIsim also includes:

- **`ss.MaternalNet`**: Mother-to-child transmission network (created automatically)
- **`sti.AgeMatchedMSM`**: Men-who-have-sex-with-men network with age-based matching
- **`sti.PriorPartners`**: Stores former partners for contact tracing (opt-in via `recall_prior=True`)
