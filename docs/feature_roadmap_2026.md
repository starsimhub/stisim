# STIsim Feature Map & Roadmap: Q1 2026 -- Q2 2027

## Context

22 upcoming HIVsim/STIsim research projects (from `docs/HIVsim projects.xlsx`), 39 open GitHub issues, and cross-cutting model hardening goals. This plan delivers: per-project feature maps showing what exists vs what's needed, a project-driven quarterly roadmap (2-week sprints, 2-3 person team), and GitHub Projects board setup.

## Implementation Steps

1. **Write the full roadmap into `docs/feature_planning.md`** -- replaces the current stub
2. **Set up GitHub Projects board** -- "STIsim Roadmap 2026-2027" (via `gh` CLI, separate step)

---

## 1. Priority 1 Projects: Feature Maps

### P1.1: Syphilis diagnostics Zimbabwe -- COMPLETED (Q1 2026)
**Lead**: Robyn Stuart | **Repo**: `syph_dx_zim` | **Status**: Manuscript in revision

| Category | Module | Status |
|----------|--------|--------|
| **Diseases** | HIV | CORE |
| | Syphilis (multi-stage: primary -> secondary -> early/late latent -> tertiary) | CORE |
| | GUD placeholder | CORE |
| **Networks** | StructuredSexual (risk groups, FSW/clients, condom dynamics) | CORE |
| | MaternalNet | STARSIM |
| **Connectors** | hiv_syph (bidirectional susceptibility/transmissibility) | CORE |
| **Interventions** | ART (coverage data, CD4-based initiation) | CORE |
| | VMMC | CORE |
| | PrEP (basic) | CORE |
| | HIVTest (FSW-targeted, general pop, low-CD4) | CORE |
| | SyphTest, SyphDx, ProductMix, ANCSyphTest | CORE |
| **Analyzers** | coinfection_stats, sw_stats, partner_age_diff | CORE |
| **Calibration** | Optuna with weighted likelihoods, survival filtering, seed tracking | CORE |
| **Custom (project-local)** | DALY analyzers (incidence + hybrid) | PROJECT |
| | Treatment pathway attribution (GUD syndromic, ANC, KP dual, PLHIV, newborn) | PROJECT |
| | Transmission-by-stage accounting (sexual vs MTC, birth outcomes by stage) | PROJECT |
| | DualTest (cascading HIV/syphilis testing) | PROJECT |
| | Pregnancy risk reduction intervention | PROJECT |
| | Network snapshots | PROJECT |

**Takeaway**: This project exercised nearly the entire STIsim core. The DALY analyzers and treatment pathway attribution are candidates for upstreaming.

---

### P1.2: Novel POC diagnostics for demand generation -- HIGHEST ACTIVE PRIORITY
**Leads**: Robyn Stuart + LSHTM | **Partners**: WHI + LSHTM | **Deliver by**: End Q2 2026
**Pathogens**: NG + CT + TV + syphilis
**Question**: How might novel POC diagnostics generate demand for STI testing, via partner notification or otherwise? Also: how avoidance of unnecessary PN might reduce GBV.

Builds on the `stisim_vddx_zim` project (POC diagnostics for vaginal/urethral discharge) which already uses:

| Category | Module | Status |
|----------|--------|--------|
| **Diseases** | HIV | CORE |
| | Gonorrhea (NG) | CORE |
| | Chlamydia (CT) | CORE |
| | Trichomoniasis (TV) | CORE |
| | SimpleBV | CORE |
| **Networks** | StructuredSexual (risk groups, FSW, condoms) | CORE |
| | MaternalNet | STARSIM |
| **Connectors** | hiv_ng, hiv_ct, hiv_tv | CORE |
| **Interventions** | HIVTest, ART, VMMC, PrEP | CORE |
| | SymptomaticTesting (POC scenario) | CORE |
| | GonorrheaTreatment, STITreatment (CT, TV) | CORE |
| | SyndromicMgmt (differential treatment by cervical infection) | PROJECT (needs upstream) |
| **Analyzers** | sw_stats | CORE |
| | NetworkDegree, RelationshipDurations, DebutAge, partner_age_diff | CORE |
| **Calibration** | Optuna (disease symptomaticity, care-seeking, beta) | CORE |
| **Additional needed** | PartnerNotification (for demand generation scenarios) | CORE |
| | GBV linkage modeling (unnecessary PN -> GBV risk) | **NEW** |
| | Syphilis integration (adding syph to existing NG/CT/TV model) | CORE (disease exists) |

**Gap**: Only real missing feature is **GBV linkage modeling** -- quantifying how better diagnostics (fewer false positives) reduce unnecessary partner notification and downstream GBV risk.

---

### P1.3: STI testing to target HIV prevention (Uganda)
**Lead**: Adam Akullian | **Partners**: IDM / WHI | **Deliver by**: Q3 2026 (intermediary), Q1 2027 (full)
**Initial tasks**: Create repo, project plan, gather data, set up baseline model

| Category | Module | Status |
|----------|--------|--------|
| **Diseases** | HIV | CORE |
| | Syphilis, NG, CT, TV (subset TBD) | CORE |
| **Networks** | StructuredSexual (risk groups, FSW, condoms) | CORE |
| | MaternalNet | STARSIM |
| | PriorPartners (for contact tracing) | CORE |
| **Connectors** | hiv_syph, hiv_ng, hiv_ct, hiv_tv (as needed) | CORE |
| **Interventions** | HIVTest, ART, VMMC | CORE |
| | STITest, SymptomaticTesting, STITreatment | CORE |
| | PartnerNotification | CORE |
| | PrEP (flexible targeting, cascade, adherence) | **NEEDS WORK** (#329-#331) |
| | Cross-intervention referral (STI+ -> HIV test -> PrEP) | EXISTS via eligibility patterns |
| **Analyzers** | sw_stats, coinfection_stats | CORE |
| **MTCT pathway** | ANC HIV testing, breastfeeding Tx, PMTCT | **NEEDS WORK** (milestone 18) |
| **Calibration** | Optuna | CORE |
| **Demographics** | Pregnancy, Deaths, Migration | CORE/STARSIM |

**Gaps**: PrEP system (#329-#334) and MTCT pathway (milestone 18). Cross-intervention referral already exists via eligibility: (1) IntvA sets `self.positive[uids]=True`, IntvB checks eligibility on this, or (2) IntvA directly sets eligibility for B.

---

### P1.4: ANC screening for asymptomatic STIs (SA, Kenya, Nigeria, Ethiopia)
**Lead**: Robyn Stuart | **Partners**: IDM + WHO | **Deliver by**: Q3 2026
**Main gap**: VOI (Value of Information) pipelines

| Category | Module | Status |
|----------|--------|--------|
| **Diseases** | HIV, NG, CT, TV, BV | CORE |
| **Networks** | StructuredSexual, MaternalNet | CORE/STARSIM |
| **Connectors** | hiv_ng, hiv_ct, hiv_tv, hiv_simplebv | CORE |
| **Interventions** | HIVTest, ART | CORE |
| | SyndromicMgmt | PROJECT (needs upstream, #345) |
| | ANCScreen (GA-windowed, multi-disease) | PROJECT (needs upstream) |
| | STIPartnerNotification (ANC-linked) | PROJECT (needs upstream) |
| | GonorrheaTreatment, STITreatment (CT, TV) | CORE |
| **Connectors** | sti_fetal (disease -> fetal health damage) | PROJECT (needs upstream) |
| **Analyzers** | birth_outcome_dalys, intervention_costs | PROJECT (needs upstream) |
| | pregnancy_sti_stats | PROJECT |
| **VOI pipeline** | Paired CRN sims, NMB, EVPI/EVPPI | PROJECT (needs upstream) |
| **Calibration** | Optuna (15 epi params -> 200 posterior sets) | CORE |
| **Demographics** | Pregnancy (with FetalHealth), Deaths | STARSIM |

**Gaps**: The model works end-to-end in `anc_sti_screening` repo. Main upstream needs: SyndromicMgmt, ANCScreen, cost/DALY analyzers, VOI helpers (paired CRN runs, NMB computation, EVPI/EVPPI).

---

### P1.5: HIVsim Eswatini validation (replicate EMOD-HIV)
**Lead**: Daniel Citron | **Deliver by**: Q1 2027 (intermediary goals TBD)
**Core aim**: Show HIVsim can replicate any analysis EMOD-HIV would do

| Category | Module | Status |
|----------|--------|--------|
| **Diseases** | HIV (CD4 tracking, ART states, acute/latent/falling) | CORE |
| **Networks** | StructuredSexual (risk groups, FSW/clients, concurrency) | CORE |
| | MaternalNet | STARSIM |
| **Interventions** | HIVTest (FSW, general, low-CD4) | CORE |
| | ART (coverage data, CD4-based) | CORE |
| | PrEP (needs full system: flexible targeting, LAI, cascade) | **NEEDS WORK** (#329-#334) |
| | VMMC | CORE |
| **Analyzers** | hiv_epi (UNAIDS/PHIA age-sex binning) | PROJECT |
| | NetworkSnapshot | PROJECT |
| | sw_stats | CORE |
| **Calibration** | Optuna (6 params: beta, eff_condom, risk props, concurrency, ART duration) | CORE |
| **Demographics** | Eswatini age dist, ASFR, deaths, migration | CORE |
| **MTCT** | Breastfeeding Tx, ANC HIV testing, perinatal progression | **NEEDS WORK** (milestone 18) |
| **Location system** | `Sim(location='eswatini')` | **PARTIAL** (#150) |

**Gaps**: Full PrEP system, MTCT pathway, location system maturation. The `hivsim_eswatini` repo already has a working baseline model with custom analyzers.

---

### P1.6: DoxyPEP (Emilia)
**Deliver by**: Q4 2026

| Category | Module | Status |
|----------|--------|--------|
| **Diseases** | NG, CT, Syphilis, HIV | CORE |
| **Networks** | StructuredSexual | CORE |
| **Connectors** | hiv_ng, hiv_ct, hiv_syph | CORE |
| **Interventions** | DoxyPEP (post-exposure doxycycline) | **NEW** |
| | GonorrheaTreatment (AMR tracking via rel_treat) | CORE |
| | STITreatment, STITest | CORE |
| **AMR** | Generalized AMR (fitness costs, multi-drug, cross-pathogen) | **NEEDS WORK** |
| **Calibration** | Optuna | CORE |

**Gaps**: DoxyPEP intervention class and generalized AMR framework (currently NG-only).

---

## 2. Feature Readiness Summary

This table shows the % of required modules that already exist for each project:

| Project | Diseases | Networks | Connectors | Interventions | Analyzers | Calibration | Overall |
|---------|----------|----------|------------|---------------|-----------|-------------|---------|
| P1.1 Syph Dx Zim | 100% | 100% | 100% | 100% | 100% | 100% | **COMPLETE** |
| P1.2 POC Diagnostics | 100% | 100% | 100% | 90% | 100% | 100% | **~95%** (GBV linkage missing) |
| P1.3 Uganda HIV | 100% | 100% | 100% | 70% | 100% | 100% | **~85%** (PrEP + MTCT) |
| P1.4 ANC VOI | 100% | 100% | 100% | 80% | 70% | 100% | **~80%** (upstream + VOI) |
| P1.5 Eswatini | 100% | 100% | N/A | 70% | 80% | 100% | **~80%** (PrEP + MTCT + location) |
| P1.6 DoxyPEP | 100% | 100% | 100% | 50% | 100% | 100% | **~75%** (DoxyPEP + AMR) |

**Key message for leadership**: 75-95% of what's needed for every P1 project already exists in core STIsim. The remaining work is primarily: PrEP system buildout (shared across P1.3/P1.5), MTCT pathway (shared across P1.3/P1.5), upstreaming project-local code (P1.4), and two new intervention classes (GBV linkage, DoxyPEP).

---

## 3. Model Hardening (Cross-Cutting, Integrated Throughout)

| Dimension | Current State | Actions | When |
|---|---|---|---|
| **Documentation** | Docs site exists but incomplete | Update API docs with each release | Every release |
| **Test coverage** | HIV sensitivity tests partial (#247-#250, #239, #223, #277) | Complete network + MTCT + U=U tests | Q3-Q4 2026 |
| **Network validation** | Analyzers exist (RelationshipDurations, NetworkDegree, etc.) | Write tests asserting distributions match targets | Q2-Q3 2026 |
| **Multi-country** | Location arg partial (#150) | Country param databases, auto calibration targets | Q1 2027 |
| **Performance** | Known slow spots (#108, #168) | Profile, optimize, reduce test runtime (#362) | Ongoing |
| **Infrastructure** | Coverage parsing local to HIV (#328), interpolation inconsistent (#327) | Extract shared utilities | Q2 2026 |
| **Eligibility bug** | #353 open | Audit and fix all HIV interventions | Q2 2026 Sprint 1 |

---

## 4. Quarterly Roadmap

### Q2 2026 (Apr--Jun): P1.2 Delivery + Foundation -- Release v1.6.0

**Theme**: Deliver P1.2 (POC diagnostics/GBV), upstream project-local code, build shared infrastructure.

**Sprint 1-2 (Apr 7 -- May 2): P1.2 Core + Infrastructure**
- **Stream A** (P1.2): GBV linkage modeling, POC diagnostic scenarios with syphilis, upstream SyndromicMgmt (#345)
- **Stream B** (infrastructure): Shared coverage-targeting (#328), consistent interpolation (#327), fix eligibility (#353), CT eff_condom (#364)

**Sprint 3-4 (May 5 -- May 30): P1.2 Completion + Upstream**
- **Stream A** (P1.2): Integration testing, scenario runs, results analysis
- **Stream B** (upstream): ANCScreen into core, consolidate PN implementations, `current_sw` property (#307)

**Sprint 5-6 (Jun 2 -- Jun 27): PrEP Phase 1 + P1.3 Setup**
- **Stream A** (PrEP): Flexible targeting (#330), oral cascade (#331), adherence testing (#252)
- **Stream B** (P1.3 setup): Create Uganda repo, project plan, gather data, stub baseline model

**Release v1.6.0** (end Jun)

---

### Q3 2026 (Jul--Sep): P1.3/P1.4 Delivery + PrEP/MTCT -- Release v1.7.0

**Theme**: Deliver P1.4 (ANC VOI), advance P1.3 (Uganda), complete PrEP + MTCT.

**Sprint 7-8 (Jul 6 -- Jul 31): PrEP Phase 2 + MTCT Start**
- **Stream A** (PrEP): PrEP on demand (#332), PEP (#333), long-acting injectable (#334)
- **Stream B** (MTCT): ANC HIV testing pathway (#322), breastfeeding transmission (#323)

**Sprint 9-10 (Aug 4 -- Aug 29): MTCT Completion + P1.4 VOI**
- **Stream A** (MTCT): PMTCT/flexible ART (#324), MTCT tracking (#325), perinatal progression (#224)
- **Stream B** (P1.4): VOI pipeline upstream (cost/DALY analyzers, CRN helpers), country params for SA/Kenya

**Sprint 11-12 (Sep 1 -- Sep 26): P1.3 Integration + P1.4 Delivery**
- **Stream A** (P1.3): Uganda baseline model, initial calibration, STI -> HIV prevention workflow
- **Stream B** (P1.4): EVPI/EVPPI analysis for WHO, manuscript

**Release v1.7.0** (end Sep)

---

### Q4 2026 (Oct--Dec): P1.6 DoxyPEP + Validation -- Release v1.8.0

**Theme**: DoxyPEP, AMR expansion, systematic validation, P1.5 groundwork.

**Sprint 13-14 (Oct 5 -- Oct 31): DoxyPEP + POC Products**
- **Stream A** (P1.6): DoxyPEP intervention, AMR interaction modeling
- **Stream B** (dx): POC product library, integration with SyndromicMgmt/ANCScreen

**Sprint 15-16 (Nov 3 -- Nov 28): AMR + Validation**
- **Stream A** (AMR): Generalize beyond NG, drug-specific efficacy, multi-drug tracking
- **Stream B** (validation): HIV sensitivity tests (#247-#250, #239, #223, #277), network validation

**Sprint 17-18 (Dec 1 -- Dec 26): P1.5 Groundwork + Performance**
- **Stream A** (P1.5): Eswatini baseline model, EMOD-HIV comparison scenarios
- **Stream B** (performance): Optimization (#108, #168, #362)

**Release v1.8.0** (end Dec)

---

### Q1 2027 (Jan--Mar): P1.5 Delivery + Location System -- Release v2.0.0

**Theme**: HIVsim Eswatini validation, multi-country, P1.3 completion.

**Sprint 19-20 (Jan 4 -- Jan 29): Location System + P1.5**
- **Stream A** (location): Location argument (#150), data loading (#182), country param sets
- **Stream B** (P1.5): Eswatini calibration, EMOD-HIV comparison

**Sprint 21-22 (Feb 1 -- Feb 26): Validations**
- **Stream A** (P1.5): Eswatini validation write-up, scenario replication
- **Stream B**: Zambia + Kenya validation, age bin (#326), ART logic (#335)

**Sprint 23-24 (Mar 1 -- Mar 26): P1.3 Completion + Polish**
- **Stream A** (P1.3): Uganda model finalization, scenario analysis
- **Stream B**: CD4 analyzer (#340), mortality adjustment (#352), docs, release prep

**Release v2.0.0** (end Mar)

---

## 5. Dependency Map

```
P1.2 (POC diagnostics/GBV, Q2) -- HIGHEST PRIORITY
  |-- SyndromicMgmt upstream (#345)
  |-- GBV linkage modeling (NEW)
  +-- PN refinement (benefits P1.3, P3 projects)

Shared Infrastructure [Q2 Sprint 1-2]
  |-- Coverage logic (#328)
  |-- Interpolation (#327)
  +-- Eligibility fix (#353)

PrEP System (#329-#334) [Q2-Q3]
  +-- P1.3 Uganda HIV Prevention
  +-- P1.5 Eswatini validation
  +-- Maternal PrEP (#360) -> MTCT

MTCT Pathway (#321-#325) [Q3]
  +-- P1.3 Uganda
  +-- P1.5 Eswatini
  +-- MTCT rate validation (#277)

DoxyPEP + AMR [Q4]
  +-- P1.6 (Emilia)

Location System (#150) [Q1 2027]
  +-- P1.5 Eswatini
  +-- Multi-country deployment
```

**Critical path**: SyndromicMgmt upstream -> P1.2 delivery (Q2) -> PrEP Phase 1 -> MTCT -> P1.3/P1.5 integration (Q3-Q1 2027)

---

## 6. Project Tracking: GitHub Projects Board

### Setup
1. **Board**: "STIsim Roadmap 2026-2027" on starsimhub/stisim
   - Columns: Backlog | Sprint Ready | In Progress | In Review | Done
2. **Milestones**: v1.6.0 (Q2), v1.7.0 (Q3), v1.8.0 (Q4), v2.0.0 (Q1 2027)
3. **Labels to add**: `upstream`, `infrastructure`, `validation`, `P1`/`P2`/`P3`, `model-hardening`
4. **Sprint iterations**: 2-week cadence via GitHub Projects Iteration field
5. **Views**: By milestone, by project, by stream (A/B)

### New Issues to Create
- "GBV linkage modeling for unnecessary partner notification" (P1.2)
- "Upstream SyndromicMgmt from anc_sti_screening" (#345 refinement)
- "Upstream ANCScreen from anc_sti_screening"
- "Upstream VOI pipeline components (cost/DALY analyzers, CRN helpers)"
- "DoxyPEP intervention" (P1.6)
- "Generalize AMR framework beyond NG"
- "P1.3 Uganda: repo setup and project plan"
- "P1.5 Eswatini: identify EMOD-HIV comparison scenarios"

---

## 7. Key Files

- [base_interventions.py](stisim/interventions/base_interventions.py) -- destination for SyndromicMgmt, ANCScreen, PN consolidation
- [hiv_interventions.py](stisim/interventions/hiv_interventions.py) -- PrEP rewrite, `parse_coverage()` extraction, ART PMTCT
- [hiv.py](stisim/diseases/hiv.py) -- breastfeeding transmission, perinatal progression, MTCT tracking
- [networks.py](stisim/networks.py) -- `current_sw` property
- [anc_sti_screening/interventions.py](/Users/robynstuart/gf/anc_sti_screening/interventions.py) -- SyndromicMgmt, ANCScreen, STIPartnerNotification source
- [anc_sti_screening/run_voi.py](/Users/robynstuart/gf/anc_sti_screening/run_voi.py) -- VOI pipeline reference
- [anc_sti_screening/analyzers.py](/Users/robynstuart/gf/anc_sti_screening/analyzers.py) -- cost/DALY analyzers to upstream
- [stisim_vddx_zim/model.py](/Users/robynstuart/gf/stisim_vddx_zim/model.py) -- P1.2 base model
- [stisim_vddx_zim/interventions.py](/Users/robynstuart/gf/stisim_vddx_zim/interventions.py) -- SyndromicMgmt source (VDS/UDS variant)
- [hivsim_eswatini/](/Users/robynstuart/gf/hivsim_eswatini/) -- P1.5 baseline model
- [calibration.py](stisim/calibration.py) -- calibration framework
