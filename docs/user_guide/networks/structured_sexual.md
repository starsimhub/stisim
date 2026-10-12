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

## Coital acts and condoms

Coital acts per partnership are set by `acts`, and condom use by `condom_data` (see [Parameters](#parameters)). Transmission probability per partnership accounts for both condom use and number of acts:

P(transmission) = 1 - [1 - beta * (1 - eff_condom)]^(acts * p_condom) * [1 - beta]^(acts * (1 - p_condom))

## Parameters

Pars of `sti.StructuredSexual`, including those inherited from `sti.MFNetwork`, `sti.SWNetwork`, and `sti.BaseNetwork`. The Source column follows the conventions of the [parameter reference](../../parameters.qmd).

### Sexual debut

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `debut_f` | lognorm(20 yr, 3 yr) | Age at sexual debut, females | Not recorded |
| `debut_m` | lognorm(21 yr, 3 yr) | Age at sexual debut, males | Not recorded |

### Risk groups

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `n_risk_groups` | 3 | Number of risk groups | Not applicable |
| `prop_f0` | 0.85 | Proportion of females in the low risk group | Not recorded |
| `prop_m0` | 0.80 | Proportion of males in the low risk group | Not recorded |
| `prop_f2` | 0.01 | Proportion of females in the high risk group | Not recorded |
| `prop_m2` | 0.02 | Proportion of males in the high risk group | Not recorded |
| `f_age_group_bins` | teens: 0-20, young: 20-25, adult: 25+ | Female age groups used to stratify age mixing and partnership durations | Not applicable |

### Concurrency

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `concurrency_dist` | poisson | Distribution of the number of additional concurrent partners sought; its mean is set by the values below | Not applicable |
| `f0_conc` | 0.0001 | Mean concurrency, low-risk females | Not recorded |
| `f1_conc` | 0.01 | Mean concurrency, medium-risk females | Not recorded |
| `f2_conc` | 0.1 | Mean concurrency, high-risk females | Not recorded |
| `m0_conc` | 0.0001 | Mean concurrency, low-risk males | Not recorded |
| `m1_conc` | 0.2 | Mean concurrency, medium-risk males | Not recorded |
| `m2_conc` | 0.5 | Mean concurrency, high-risk males | Not recorded |

### Partnership formation and duration

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `age_diff_pars` | teens: (7, 3), (6, 3), (5, 1); young and adult: (8, 3), (7, 3), (5, 2) | Preferred partner age difference (male minus female, years), normal (mean, sd) by female age group and risk group | Not recorded |
| `match_method` | 'closest_age_tapered_seeking' | Pair-matching algorithm (see Pair matching) | Not applicable |
| `f_partnership_taper_cut` | 55 | Age above which females stop seeking new partnerships | Not recorded |
| `p_pair_form` | 0.5 | Probability that a matched pair forms a partnership | Not recorded |
| `p_matched_stable` | [0.9, 0.5, 0] | Probability that a same-risk-group partnership is stable rather than casual, by risk group | Not recorded |
| `p_mismatched_casual` | [0.5, 0.5, 0.5] | Probability that a mixed-risk-group partnership is casual rather than one-time, by risk group | Not recorded |
| `stable_dur_pars` | RG0: [100 yr, 1 yr]; RG1: teens [8 yr, 2 yr], young [10 yr, 3 yr], adult [12 yr, 3 yr]; RG2: effectively zero (1e-4 mo) | Duration of stable partnerships, lognormal [mean, sd], by female age group and risk group | Not recorded |
| `casual_dur_pars` | [1 yr, 3 yr] for all groups | Duration of casual partnerships, lognormal [mean, sd], by female age group and risk group | Not recorded |

### Sex work

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `fsw_shares` | 0.05 | Lifetime proportion of women who are sex workers | Not recorded |
| `client_shares` | 0.12 | Lifetime proportion of men who are clients | Not recorded |
| `age_sw_start` | normal(20 yr, 3 yr) | Age at entry into sex work | Not recorded |
| `dur_sw` | lognorm(5 yr, 3 yr) | Duration of sex work | Not recorded |
| `age_client_start` | normal(25 yr, 5 yr) | Age at becoming a client | Not recorded |
| `dur_client` | lognorm(10 yr, 5 yr) | Duration as a client | Not recorded |
| `sw_seeking_rate` | 1/month | Probability per month that an active client seeks an FSW | Not recorded |
| `sw_intensity` | uniform(0, 1) | Per-step working intensity of each FSW, weighting client allocation | Not recorded |
| `sw_beta` | 1 | Multiplier on transmission probability for sex-work edges | Not applicable |
| `fsw_mf_conc_mult` | 1.0 | Multiplier on FSW non-sex-work concurrency; values below 1 give active sex workers fewer non-sex-work partners | Not applicable |

### Coital acts and condoms

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| `acts` | lognorm(80/yr, 30/yr) | Coital acts per partnership per year | Not recorded |
| `condom_data` | None | Condom use probability: scalar, or time series stratified by risk-group pair and sex-work edges | Not applicable |
| `condom_smoothness` | None | Smoothing applied when interpolating time-varying `condom_data` | Not applicable |
| `recall_prior` | False | Whether former partners are stored for partner notification (requires `sti.PriorPartners`) | Not applicable |
## Other networks

STIsim also includes:

- **`sti.MFNetwork`** / **`sti.SWNetwork`**: the heterosexual and sex-work layers, usable standalone (see [Architecture](#architecture)).
- **`ss.MaternalNet`**: Mother-to-child transmission network (created automatically).
- **MSM networks**: men-who-have-sex-with-men contact, with several matching strategies — see [MSM networks](msm.md).
- **`sti.PriorPartners`**: Stores former partners for contact tracing. Add it alongside an MF network and set `recall_prior=True`.
