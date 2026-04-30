"""
PrEP API vignettes — design sandboxing for issue #431 / v1.6 milestone.

Guiding principle: "common things should be simple" (Starsim philosophy).
Eligibility is always explicit — there is no default target population, because
any implicit default (e.g. FSW) is not self-explanatory and will surprise users.

Design decisions captured here:
  - on_prep lives on the HIV disease object (not on each intervention), so that
    multiple concurrent PrEP interventions can coordinate without a manager class.
    Only one shared BoolArr needed: hiv.on_prep. Per-product state (cost, timing)
    stays on the intervention.
  - groups= API (V2b+) allows a single Prep instance to manage multiple target
    populations with separate coverage levels and shared supply awareness.
  - product= (V3+) bundles efficacy + cost; cost tracking enables CE analysis.
  - supply= (V4+) enforces a dose budget; allocation is gap-proportional.

Status per vignette:
  V1, V2a  — run today against rc1.5.4 (Prep class as-is, with explicit eligibility)
  V2b–V5   — API sketches; commented-out or pass-bodied; marked NOTIMPLEMENTED

Clark's PR #432 approach is shown in V3b as an alternative for comparison.
"""

import sciris as sc
import starsim as ss
import stisim as sti


# =============================================================================
# V1 — Simple single-group PrEP
# Runs today. Explicit eligibility required.
# =============================================================================

def v1_fsw():
    """Single group, constant coverage, default oral efficacy (80%)."""
    hiv = sti.HIV()
    net = sti.StructuredSexual()
    prep = sti.Prep(
        coverage=0.5,
        eligibility=lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected,
    )
    sim = ss.Sim(n_agents=500, diseases=hiv, networks=net, interventions=prep, dt=1/12)
    sim.run()
    return sim


def v1b_agyw():
    """Single group, time-varying coverage ramp, explicit AGYW eligibility."""
    hiv = sti.HIV()
    net = sti.StructuredSexual()
    prep = sti.Prep(
        coverage={'year': [2015, 2020], 'value': [0, 0.4]},
        eligibility=lambda sim: (
            sim.people.female & ~sim.diseases.hiv.infected & (sim.people.age < 25)
        ),
    )
    sim = ss.Sim(n_agents=500, diseases=hiv, networks=net, interventions=prep, dt=1/12)
    sim.run()
    return sim


# =============================================================================
# V2 — Multi-group targeting
#
# V2a: multiple Prep instances (works today; name= required to avoid collision).
#   Limitation: if supply is shared across groups, allocation is uncoordinated.
#
# V2b: groups= API on a single Prep instance (NOTIMPLEMENTED).
#   Single instance sees all groups; can enforce shared supply constraints.
# =============================================================================

def v2a_multigroup_instances():
    """
    Multiple Prep instances, one per group. Runs today.
    Each manages its own coverage independently — no cross-group supply awareness.
    """
    hiv = sti.HIV()
    net = sti.StructuredSexual()
    prep_fsw = sti.Prep(
        name='prep_fsw',
        coverage=0.6,
        eligibility=lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected,
    )
    prep_agyw = sti.Prep(
        name='prep_agyw',
        coverage=0.3,
        eligibility=lambda sim: (
            sim.people.female & ~sim.diseases.hiv.infected & (sim.people.age < 25)
        ),
    )
    sim = ss.Sim(n_agents=2000, diseases=hiv, networks=net,
                 interventions=[prep_fsw, prep_agyw], dt=1/12)
    sim.run()
    return sim


def v2b_multigroup_groups_api():
    """
    NOTIMPLEMENTED — groups= API: one Prep instance, multiple target groups.

    Design intent: each group dict has label/eligibility/coverage; the Prep class
    tracks per-group coverage and (when supply= is set) allocates doses gap-proportionally.
    Avoids the need for multiple instances + manual name management.
    """
    # NOTIMPLEMENTED: groups= kwarg on Prep
    #
    # prep = sti.Prep(
    #     groups=[
    #         dict(label='fsw',
    #              eligibility=lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected,
    #              coverage=0.6),
    #         dict(label='agyw',
    #              eligibility=lambda sim: sim.people.female & ~sim.diseases.hiv.infected & (sim.people.age < 25),
    #              coverage=0.3),
    #         dict(label='pbfw',
    #              eligibility=lambda sim: sim.demographics.pregnancy.pregnant & ~sim.diseases.hiv.infected,
    #              coverage=0.5),
    #     ],
    # )
    # sim = ss.Sim(..., interventions=prep); sim.run()
    pass


# =============================================================================
# V3 — Product-aware PrEP (NOTIMPLEMENTED)
#
# Product bundles: name, efficacy, cost (per dose or per person-year), duration.
# Duration is relevant for injectables (lenacapavir: 6-month injection).
# Cost is required for CE analysis.
#
# V3a: product= as a simple dict on Prep (proposed minimal API).
# V3b: Clark's Product/Supply/Supplies class objects (PR #432 API).
# =============================================================================

def v3a_product_dict_api():
    """
    NOTIMPLEMENTED — product= dict on Prep.

    Efficacy from product= overrides eff_prep. Cost tracked in sim.results.
    Duration is optional; if omitted, treated as continuous (oral PrEP behaviour).

    Note on on_prep coordination: with two Prep instances both writing to
    hiv.rel_sus, they must not double-apply. If on_prep is on the HIV object
    (hiv.on_prep), each instance can check ~hiv.on_prep before uptake —
    no PrepManager needed.
    """
    # NOTIMPLEMENTED: product= on Prep
    #
    # oral = sti.Prep(
    #     name='oral',
    #     coverage=0.5,
    #     eligibility=lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected,
    #     product=dict(name='oral_tfv_ftc', efficacy=0.80, cost_per_year=50),
    # )
    # lnc = sti.Prep(
    #     name='lnc',
    #     coverage=0.5,
    #     eligibility=lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected,
    #     product=dict(name='lenacapavir', efficacy=0.999, cost_per_year=40_000,
    #                  duration=ss.dur(6, 'months')),  # injection cycle
    # )
    # sim = ss.Sim(..., interventions=[oral, lnc]); sim.run()
    # print(sim.results.oral.total_cost[-1])
    # print(sim.results.lnc.total_cost[-1])
    pass


def v3b_product_clark_api():
    """
    NOTIMPLEMENTED — Clark's Product/Supply/Supplies API from PR #432.

    Requires stisim.interventions.{product,supply,supplies,prep,prep_manager}.
    Shown here for design comparison; these modules are on branch '330' only.

    Open questions with this approach:
      1. PrepManager must be added as a separate intervention — easy to forget,
         silently wrong if omitted (agents never drop off PrEP).
      2. PrEP states (on_prep, ti_prep_end, prep_source, …) live on hiv.py —
         couples the disease model to intervention internals.
      3. random.sample used for agent selection — breaks CRN (use ss.bernoulli).
      4. Reuptake logic fully disabled (commented out, pending design decision).
      5. base_care_seeking_rate = 0.8 is hardcoded / arbitrary.
    """
    # from stisim.interventions.product import Product
    # from stisim.interventions.supply import Supply
    # from stisim.interventions.supplies import Supplies
    # from stisim.interventions.prep import SuppliedPrep
    # from stisim.interventions.prep_manager import PrepManager
    #
    # product  = Product(name='lenacapavir', type='prep', delivery_mode='injectable',
    #                    cost=40_000, eff_by_ti=[0.999])
    # supply   = Supply(quantity=np.inf, product=product)
    # supplies = Supplies(supplies=[supply])
    #
    # prep_manager = PrepManager()        # hidden dependency — must not be forgotten
    # prep = SuppliedPrep(
    #     name='lnc_fsw',
    #     eligibilities=[lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected],
    #     coverages=[0.6],
    #     supplies=supplies,
    # )
    # sim = ss.Sim(..., interventions=[prep_manager, prep]); sim.run()
    pass


# =============================================================================
# V4 — Supply-constrained allocation (NOTIMPLEMENTED)
#
# Fixed annual dose budget distributed across groups proportional to coverage gap.
# gap_i = (target_coverage_i - actual_coverage_i) * group_size_i
# doses_to_group_i = total_supply * gap_i / sum(gap_j)
#
# This is exactly Clark's allocation logic in SuppliedPrep.step() —
# the question is whether it belongs in Prep or in a Supplies helper.
# =============================================================================

def v4_supply_constrained():
    """
    NOTIMPLEMENTED — fixed supply of lenacapavir doses split gap-proportionally.

    With 500 doses/year and two groups at different coverage levels,
    the scarce supply is directed where the shortfall is largest.
    """
    # NOTIMPLEMENTED: supply= and groups= on Prep
    #
    # prep = sti.Prep(
    #     groups=[
    #         dict(label='fsw',  eligibility=..., coverage=0.6),
    #         dict(label='agyw', eligibility=..., coverage=0.3),
    #     ],
    #     product=dict(name='lenacapavir', efficacy=0.999, cost_per_year=40_000),
    #     supply=500,   # doses per year; gap-proportional allocation across groups
    # )
    # sim = ss.Sim(n_agents=2000, ..., interventions=prep, dt=1/12)
    # sim.run()
    # print('doses used:', prep.doses_used)
    # print('total cost:', prep.total_cost)
    pass


# =============================================================================
# V5 — Kenya lenacapavir CE: price threshold sweep (NOTIMPLEMENTED)
#
# Research question (CEMA collaborator): "At what cost per person-year would
# Kenya be willing to pay for lenacapavir PrEP?"
# Also: evaluate FSW-first vs. AGYW-first vs. combined targeting.
#
# Requires:
#   - product= cost tracking (V3)
#   - infections-averted relative to baseline (analyzer or post-hoc)
#   - Kenya-calibrated parameters (hiv_kenya localization or manual pars)
# =============================================================================

def v5_kenya_price_sweep():
    """
    NOTIMPLEMENTED — sweep lenacapavir cost per person-year, compute ICER vs. oral PrEP.
    Kenya WTP threshold: ~$1,000–2,000/DALY averted (rough estimate).
    """
    # costs_per_year = [50, 200, 500, 1_000, 5_000, 10_000, 40_000]
    #
    # # Baseline: oral PrEP at current coverage
    # baseline = make_kenya_sim(interventions=[
    #     sti.Prep(name='oral', coverage=0.3,
    #              eligibility=..., product=dict(name='oral', efficacy=0.80, cost_per_year=50)),
    # ])
    # baseline.run()
    #
    # results = []
    # for cost in costs_per_year:
    #     sim = make_kenya_sim(interventions=[
    #         sti.Prep(name='lnc', coverage=0.3,
    #                  eligibility=...,
    #                  product=dict(name='lenacapavir', efficacy=0.999, cost_per_year=cost)),
    #     ])
    #     sim.run()
    #     averted = baseline.results.hiv.new_infections[:].sum() - sim.results.hiv.new_infections[:].sum()
    #     incremental_cost = sim.results.lnc.total_cost[-1] - baseline.results.oral.total_cost[-1]
    #     results.append(dict(cost=cost, averted=averted, icer=incremental_cost/averted))
    #
    # # Plot: ICER vs. cost/year; horizontal line at Kenya WTP threshold
    pass


def v5b_kenya_targeting_comparison():
    """
    NOTIMPLEMENTED — fixed lenacapavir supply; compare targeting strategies.
    Same total doses: FSW-only vs. AGYW-only vs. FSW+AGYW vs. PBFW-only.
    """
    # fsw_elig  = lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected
    # agyw_elig = lambda sim: sim.people.female & ~sim.diseases.hiv.infected & (sim.people.age < 25)
    # pbfw_elig = lambda sim: sim.demographics.pregnancy.pregnant & ~sim.diseases.hiv.infected
    #
    # strategies = {
    #     'fsw_only':  [dict(label='fsw',  eligibility=fsw_elig,  coverage=0.8)],
    #     'agyw_only': [dict(label='agyw', eligibility=agyw_elig, coverage=0.4)],
    #     'fsw_agyw':  [dict(label='fsw',  eligibility=fsw_elig,  coverage=0.5),
    #                   dict(label='agyw', eligibility=agyw_elig, coverage=0.3)],
    #     'pbfw_only': [dict(label='pbfw', eligibility=pbfw_elig, coverage=0.6)],
    # }
    # for label, groups in strategies.items():
    #     prep = sti.Prep(groups=groups, product='lenacapavir', supply=500)
    #     sim  = make_kenya_sim(interventions=[prep])
    #     sim.run()
    pass


# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    sc.heading('V1: simple FSW PrEP (runs today)')
    v1_fsw()
    v1b_agyw()

    sc.heading('V2a: multi-group via multiple instances (runs today)')
    v2a_multigroup_instances()

    sc.heading('V2b–V5: API sketches (NOTIMPLEMENTED — see comments)')
