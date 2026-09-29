# Data sources for HIV intervention design

Where to find the data each intervention needs. Prefer machine- readable downloads and APIs over manual transcription when they exist. Where they don't, the manual path is documented.

## By intervention

### ART

| Source | Provides | Access |
|---|---|---|
| **UNAIDS AIDSInfo** (`aidsinfo.unaids.org`) | Country-year ART coverage (number and proportion of PLHIV on ART), by sex and broad age groups when country reports separately. Also viral load suppression. | Web download (CSV / Excel via the Advanced Data Search). No public REST API as of 2026. |
| **PHIA-family surveys** (SHIMS, ZIMPHIA, MPHIA, ZAMPHIA, etc.) | Cross-sectional ART coverage among PLHIV, by age × sex, at a single survey year. Also VLS. | Per-country PDF reports on `phia.icap.columbia.edu`; micro-data on request. |
| **National HIV program reports / Spectrum files** | Country-year ART on-treatment counts, VLS. | Ministry of Health websites, or via UNAIDS Spectrum outputs (usually national estimation team). |
| **PEPFAR MER (Monitoring, Evaluation, Reporting)** (`data.pepfar.gov`) | Country-quarter programmatic counts on ART, by age × sex, by facility / district / country. Only for PEPFAR-supported programs. | Panorama Spotlight web tool; some datasets available as bulk downloads. |

**Format expected by `sti.ART`**: DataFrame with columns `Year`, `Gender`, `AgeBin` (optional), and either `n_art` (counts) or `p_art` (proportions of PLHIV on ART). Both reference analyses build this from combined UNAIDS + PHIA + national program data.

### HIV testing

| Source | Provides | Access |
|---|---|---|
| **DHS Statcompiler** (`statcompiler.com`, `api.dhsprogram.com`) | Ever-tested-for-HIV and recent-testing rates in DHS-surveyed countries, by age × sex × wealth. Historical rounds go back to the 1990s. | REST API documented at `https://api.dhsprogram.com/`; see `network-data/references/dhs-statcompiler.md` in this plugin for the pattern. |
| **PHIA-family surveys** | Recent testing, awareness of status, and by-population testing coverage. | PDF reports; micro-data on request. |
| **UNAIDS AIDSInfo** | National annual HIV testing coverage (people tested, people newly diagnosed). | Web download. |
| **PEPFAR MER (HTS_TST)** | Programmatic testing counts by modality, population, facility, quarter. | `data.pepfar.gov` Spotlight; bulk downloads for some indicators. |
| **National HIV program HTS reports** | Testing yields, positivity, modality-specific breakdowns. | Ministry of Health websites; national HIV council annual reports. |

**Format expected by `sti.HIVTest`**: `test_prob_data` as an annual probability aligned to `years`. Reference analyses construct this from ramps informed by DHS + UNAIDS + program reports — see `references/testing.md` for the pattern.

**Modality-specific data considerations:**

- **General-population testing** — DHS is the strongest source; PHIA triangulates.
- **Key-population testing** — IBBS (Integrated Bio-Behavioral Surveys) and PEPFAR MER (disaggregated by key population).
- **Symptomatic / low-CD4 testing** — program reports on provider-initiated testing at clinical encounters; rarely disaggregated to CD4 explicitly, so the low-CD4 rate is usually an assumption grounded in a plausible fraction of testing.
- **ANC testing** — PEPFAR MER (PMTCT_STAT), national PMTCT program reports.
- **Infant testing** — PEPFAR MER (PMTCT_EID for early infant diagnosis).

### VMMC

| Source | Provides | Access |
|---|---|---|
| **DHS Statcompiler** | Male circumcision prevalence (self-reported), by age. Includes traditional + programmatic in DHS-surveyed countries. Indicator group: "Sexual behavior / Circumcision". | REST API. |
| **PHIA-family surveys** | Male circumcision prevalence measured or self-reported, by age × region. | PDF reports; micro-data. |
| **PEPFAR MER (VMMC_CIRC)** | Program procedures performed, by age × country × quarter. Excludes traditional. | `data.pepfar.gov` Spotlight; bulk downloads. |
| **National VMMC program reports** | Cumulative procedures, coverage targets vs achieved. | Ministry of Health websites. |

**Format expected by `sti.VMMC`**: DataFrame with `Year`, `Gender`, `AgeBin` (optional), and either `p_vmmc` (prevalence) or `n_vmmc` (counts). Because `sti.VMMC` treats coverage as a stock target (prevalence), *procedure count* data (a flow) needs conversion — sum cumulatively and divide by the male population size in each stratum to get a prevalence proxy. DHS and PHIA give prevalence directly and are the simpler starting point.

### PrEP

| Source | Provides | Access |
|---|---|---|
| **PEPFAR MER (PrEP_NEW, PrEP_CT, PrEP_CURR)** | Number initiating, on treatment, and continuing, by age × country × quarter. Only PEPFAR-supported programs. | `data.pepfar.gov` Spotlight; bulk downloads. |
| **UNAIDS AIDSInfo** | Number of people on PrEP by country-year (limited data, recent years). | Web download. |
| **National HIV program reports** | PrEP program scale-up, targets, populations served. | Ministry of Health websites; PEPFAR Country Operating Plans. |
| **Published PrEP scale-up studies** | Setting-specific program evaluations with detailed coverage, adherence, persistence. | PubMed / journals. |

**Format expected by `sti.Prep`**: same shape as ART/VMMC coverage. For efficacy and adherence, trial estimates (iPrEx, PROUD) inform `prep_eff`; setting-specific adherence studies inform `prep_adh`. Efficacy and adherence are almost always from literature, not from program data.

## Cross-cutting sources

### Existing HIVsim country analyses

The two reference analyses ship their own compiled data:

- **`hivsim_eswatini/data/`** — SHIMS-based ART/VLS by age × sex, DHS-based network parameters, PEPFAR-based VMMC (`vmmc_coverage.csv`).
- **`hivsim_zim/data/`** — ZimPHIA-based initial prevalence, UNAIDS- based ART (`p_art.csv`), national program VMMC (`n_vmmc.csv`).

For a new country analysis, browsing these gives concrete templates for what the CSVs should look like.

### API-accessible sources — summary

Only DHS Statcompiler has a documented REST API with confirmed programmatic access (see `network-data/references/dhs-statcompiler.md` for the pattern). UNAIDS AIDSInfo, PEPFAR Panorama, and PHIA reports are all web-only downloads for now. That means:

- DHS-covered indicators (testing, circumcision, network behavior) can be scripted via `requests`.
- ART, VLS, VMMC counts, PrEP counts from UNAIDS / PEPFAR / national programs require manual downloads and one-off processing scripts, ideally committed alongside the derived CSV so the transformation is reproducible (see `stisim_vddx_zim/process_ihme_data.py` for the pattern — a `process_*.py` script alongside the derived output, raw extract kept out of the repo).

### When to prefer which source

- **Multi-country consistency** matters → DHS (standardized across DHS-surveyed countries) or UNAIDS (all countries with HIV reporting).
- **Recent measured HIV outcomes** matter → PHIA (single-year, laboratory-confirmed, but only in countries where PHIA has run).
- **Program granularity** matters → PEPFAR MER (age × facility × quarter, but PEPFAR-supported programs only) or national program reports.
- **Historical scale-up trajectory** matters → UNAIDS AIDSInfo for the annualized country series.
- **Cross-check between measured and modeled** → run two sources through the workflow independently and compare — the two reference analyses do this implicitly by combining sources per indicator.

## Recording the source

For every data-fed intervention parameter, record in `analysis-brief.md`:
- Which source(s) the data came from.
- The extraction date and version if the source is API-backed (e.g. DHS survey year and API pull date).
- Any pre-processing applied (aggregation, gap-filling, interpolation).
- The `data/` filename and the `process_*.py` script (if any) that produced it.

That is the audit trail that lets a cold reader — or the calibrator six months later — reconstruct where the numbers came from.
