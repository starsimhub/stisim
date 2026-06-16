# Structured sexual network

**Class:** `sti.StructuredSexual` | **Base class:** `sti.MFNetwork`

The `StructuredSexual` network is STIsim's primary sexual contact network. It is created automatically when you make a sim, and bundles heterosexual partnerships with sex work on a single edge list.

## Architecture

STIsim's sexual networks are layered. Each network below extends `sti.BaseNetwork`, which provides the shared machinery: sexual debut, coital acts, condom-use data, and prior-partner recall.

| Class | Models | Notes |
|-------|--------|-------|
| `sti.BaseNetwork` | shared infrastructure | debut, acts, condom data, `recall_prior` |
| `sti.MFNetwork` | heterosexual partnerships | risk groups, concurrency, age mixing |
| `sti.SWNetwork` | sex work | FSW–client edges |
| `sti.StructuredSexual` | MF + sex work | the default; extends `MFNetwork` and adds sex-work edges |

`StructuredSexual` keeps the historical API, so `sim.networks.structuredsexual` exposes `fsw`, `client`, risk groups, and concurrency directly. For modular use — heterosexual partnerships without sex work, or sex work alone — use `sti.MFNetwork` or `sti.SWNetwork` directly.

## Overview

Agents are assigned a **risk group** (low, medium, high) at initialization. Risk group determines partnership patterns: concurrency, relationship type, and duration. Partnerships form each timestep based on age preferences, dissolve after their sampled duration, and are classified as stable, casual, one-time, or sex work.

## Risk groups

| Risk group | Description | Female share | Male share |
|-----------|-------------|-------------|------------|
| 0 (low) | Marry and remain monogamous | 85% | 80% |
| 1 (medium) | Marry then divorce, or have concurrent partners | 14% | 18% |
| 2 (high) | Never marry | 1% | 2% |

Shares are set by `prop_f0`, `prop_m0` (low) and `prop_f2`, `prop_m2` (high); the medium group is the remainder.

## Partnership types

How a partnership is classified depends on the risk groups of both partners:

- **Stable**: Both partners in the same risk group and pass a stability check (`p_matched_stable`: RG0 90%, RG1 50%, RG2 0%)
- **Casual**: Mismatched risk groups, 50% chance (`p_mismatched_casual`; otherwise one-time)
- **One-time**: Single-timestep contact
- **Sex work**: FSW-client contacts, one-time duration

## Concurrency

The number of concurrent partners an agent seeks depends on their sex and risk group:

| | Low (RG0) | Medium (RG1) | High (RG2) |
|---|-----------|-------------|------------|
| Female | 0.0001 | 0.01 | 0.1 |
| Male | 0.0001 | 0.2 | 0.5 |

Each value is the target *mean* concurrency. By default `concurrency_dist` is a Poisson (`lam` = the value above) and the realized count is `rvs + 1`, so every active agent seeks at least one partner. Set `concurrency_dist` to `ss.nbinom` for an overdispersed alternative — it is reparameterized from the same mean automatically.

## Age mixing

Women sample a preferred partner age from a normal distribution. Parameters vary by female age group and risk group:

| Female age group | RG0 (mean, std) | RG1 | RG2 |
|-----------------|-----------------|-----|-----|
| Teens (<20) | (7, 3) | (6, 3) | (5, 1) |
| Young (20-25) | (8, 3) | (7, 3) | (5, 2) |
| Adult (25+) | (8, 3) | (7, 3) | (5, 2) |

Values are the age difference (male minus female) in years.

### Pair matching

How sampled preferences are turned into pairs is set by `match_method` (default `'closest_age_tapered_seeking'`). See [Pair matching](matching.md) for the algorithms and how to supply your own.

## Sex work

FSW and client status are lifetime fates (`fsw_shares`, `client_shares`) gated by a per-agent active window — an entry age (`age_sw_start`) and duration (`dur_sw`) — so an agent counts as currently FSW only while inside that window.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fsw_shares` | 5% | Lifetime proportion of women who are sex workers |
| `client_shares` | 12% | Lifetime proportion of men who are clients |
| `sw_seeking_rate` | 1/month | Rate at which clients seek FSWs |
| `fsw_mf_conc_mult` | 1.0 | Multiplier on FSW non-sex-work concurrency (values <1 mean active sex workers have fewer non-sex-work partners) |

## Coital acts and condoms

| Parameter | Default | Description |
|-----------|---------|-------------|
| `acts` | lognorm(80/yr, 30/yr) | Coital acts per partnership per year |
| `condom_data` | None | Condom use probability (scalar, or stratified by risk-group pair / sex-work edges) |

Transmission probability per partnership accounts for both condom use and number of acts:

P(transmission) = 1 - [1 - beta * (1 - eff_condom)]^(acts * p_condom) * [1 - beta]^(acts * (1 - p_condom))

## Other networks

STIsim also includes:

- **`sti.MFNetwork`** / **`sti.SWNetwork`**: the heterosexual and sex-work layers, usable standalone (see [Architecture](#architecture)).
- **`ss.MaternalNet`**: Mother-to-child transmission network (created automatically).
- **MSM networks**: men-who-have-sex-with-men contact, with several matching strategies — see [MSM networks](msm.md).
- **`sti.PriorPartners`**: Stores former partners for contact tracing. Add it alongside an MF network and set `recall_prior=True`.
