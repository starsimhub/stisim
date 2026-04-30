"""
PrEP API vignettes — design sandboxing for issue #431 / v1.6 milestone.

Each vignette represents a distinct level of complexity and a concrete use case.
V1 runs today. V2–V5 are aspirational API sketches; functions are marked
NOTIMPLEMENTED where stisim doesn't yet support the interface.

The Kenya vignette (V5) is the primary research-question anchor:
"Estimate the price threshold at which Kenya would pay for lenacapavir."
Target groups under consideration: FSW, AGYW, PBFW, MSM.

Clark's PR #432 approach is shown alongside proposed alternatives in V3–V4.
"""

import numpy as np
import sciris as sc
import starsim as ss
import stisim as sti


# =============================================================================
# V1 — Simple FSW PrEP (runs today against rc1.5.4)
# =============================================================================

def v1_simple():
    """
    Simplest possible case. Should work with the current Prep class.
    One group (FSW), one coverage level, default oral efficacy.
    """
    hiv = sti.HIV()
    net = sti.StructuredSexual()
    prep = sti.Prep(coverage=0.5)   # default: targets FSW, eff_prep=0.8
    sim = ss.Sim(n_agents=500, diseases=hiv, networks=net, interventions=prep, dt=1/12)
    sim.run()
    print('V1: infections =', sim.results.hiv.new_infections[:].sum())
    return sim


def v1b_simple_with_eligibility():
    """V1 variant: override default FSW targeting with an explicit lambda."""
    hiv = sti.HIV()
    net = sti.StructuredSexual()
    # Any HIV-negative female under 30 — not just FSW
    prep = sti.Prep(
        coverage=0.4,
        eligibility=lambda sim: sim.people.female & ~sim.diseases.hiv.infected & (sim.people.age < 30),
    )
    sim = ss.Sim(n_agents=500, diseases=hiv, networks=net, interventions=prep, dt=1/12)
    sim.run()
    return sim


# =============================================================================
# V2 — Multi-group targeting (NOTIMPLEMENTED)
#
# Key requirement: different coverage levels per group; coverage gap tracked
# per group separately. Supply (if constrained) allocated proportionally.
#
# Design options:
#   A) groups= list of dicts  →  cleanest user-facing API
#   B) multiple Prep() instances  →  works today, but coverage gap logic is
#      independent per instance (no cross-group supply coordination)
# =============================================================================

def v2a_multigroup_multiple_instances():
    """
    V2 workaround: one Prep per group. No shared supply, but works today.
    Limitation: if supply is fixed across groups, allocation isn't coordinated.
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
        eligibility=lambda sim: sim.people.female & ~sim.diseases.hiv.infected & (sim.people.age < 25),
    )
    sim = ss.Sim(n_agents=2000, diseases=hiv, networks=net,
                 interventions=[prep_fsw, prep_agyw], dt=1/12)
    sim.run()
    return sim


def v2b_multigroup_groups_api():
    """
    groups= API: single Prep instance manages multiple groups with per-group coverage.
    """
    hiv = sti.HIV()
    net = sti.StructuredSexual()
    prep = sti.Prep(
        groups=[
            dict(
                label='fsw',
                eligibility=lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected,
                coverage=0.6,
            ),
            dict(
                label='agyw',
                eligibility=lambda sim: sim.people.female & ~sim.diseases.hiv.infected & (sim.people.age < 25),
                coverage=0.3,
            ),
        ],
    )
    sim = ss.Sim(n_agents=2000, diseases=hiv, networks=net, interventions=prep, dt=1/12)
    sim.run()
    return sim


# =============================================================================
# V3 — Product-aware PrEP (NOTIMPLEMENTED)
#
# Key use cases:
#   - oral TDF/FTC: daily, ~80% efficacy, ~$50/person-year, on/off at any time
#   - CAB-LA (cabotegravir): injection every 2 months, ~90% efficacy, ~$200/year
#   - lenacapavir: injection every 6 months, ~99.9% efficacy, up to $40k/year
#
# Product matters for:
#   1. Cost tracking (CE analysis)
#   2. Duration/re-enrollment (injection ≠ daily pill — can't just drop any time)
#   3. Efficacy profile (lenacapavir: near-complete protection)
#
# Two sub-vignettes showing alternative API designs:
#   V3a: simple product= dict/dataclass on existing Prep
#   V3b: Clark's Product/Supply/Supplies objects (PR #432)
# =============================================================================

def v3a_product_simple_api():
    """
    Simple product= interface: each product is a dict with {name, efficacy, cost}.
    Product efficacy overrides eff_prep parameter.
    """
    hiv = sti.HIV()
    net = sti.StructuredSexual()

    # Inline product specifications
    prep_oral = sti.Prep(
        coverage=0.5,
        product=dict(name='oral', efficacy=0.8, cost=50),
        name='prep_oral'
    )
    
    prep_lnc = sti.Prep(
        coverage=0.2,
        product=dict(name='lenacapavir', efficacy=0.999, cost=40_000),
        name='prep_lnc'
    )

    results = []
    for prep in [prep_oral, prep_lnc]:
        sim = ss.Sim(n_agents=500, diseases=hiv, networks=net, interventions=prep, dt=1/12)
        sim.run()
        results.append(dict(
            product=prep._product.get('name', 'unknown'),
            total_cost=prep._total_cost,
            infections=sim.results.hiv.new_infections[:].sum(),
        ))

    return results


def v3b_product_clark_api():
    """
    NOTIMPLEMENTED — Clark's Product/Supply/Supplies API from PR #432.
    Requires: stisim.interventions.{product,supply,supplies,prep,prep_manager}
    """
    # These imports work only on Clark's branch (PR #432):
    # from stisim.interventions.product import Product
    # from stisim.interventions.supply import Supply
    # from stisim.interventions.supplies import Supplies
    # from stisim.interventions.prep import SuppliedPrep
    # from stisim.interventions.prep_manager import PrepManager

    hiv = sti.HIV()
    net = sti.StructuredSexual()

    # Clark's interface:
    #   product = Product(name='lenacapavir', type='prep', delivery_mode='injectable',
    #                     cost=40_000, eff_by_ti=[0.999])
    #   supply  = Supply(quantity=np.inf, product=product)   # unconstrained
    #   supplies = Supplies(supplies=[supply])
    #
    #   prep_manager = PrepManager()   # must be added as its own intervention
    #   prep = SuppliedPrep(
    #       name='lnc_fsw',
    #       eligibilities=[lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected],
    #       coverages=[0.6],
    #       supplies=supplies,
    #   )
    #   sim = ss.Sim(n_agents=500, diseases=hiv, networks=net,
    #                interventions=[prep_manager, prep], dt=1/12)
    #   sim.run()

    # Questions raised by this API:
    #   - PrepManager is a hidden dependency; easy to forget, silently wrong without it
    #   - PrEP states (on_prep, ti_prep_end, etc.) live on hiv.py — couples disease to intervention
    #   - random.sample breaks CRN; should use ss.bernoulli / ss.choice
    #   - Reuptake is fully disabled (pending design decision)
    pass


# =============================================================================
# V4 — Supply-constrained, multi-group (NOTIMPLEMENTED)
#
# Fixed budget or fixed dose count distributed across groups.
# Allocation proportional to coverage gap (Clark's approach in SuppliedPrep.step()).
# Question: should allocation logic live in the Prep class or in Supplies?
# =============================================================================

def v4_supply_constrained():
    """
    NOTIMPLEMENTED — fixed annual supply of lenacapavir doses, split across FSW + AGYW.
    Annual supply of 10,000 doses; allocated proportional to per-group coverage gap.
    """
    hiv = sti.HIV()
    net = sti.StructuredSexual()

    # Desired API (option A — groups= + supply= on Prep):
    #   prep = sti.Prep(
    #       groups=[
    #           dict(label='fsw',  eligibility=..., coverage=0.6),
    #           dict(label='agyw', eligibility=..., coverage=0.3),
    #       ],
    #       product=dict(name='lenacapavir', efficacy=0.999, cost=40_000),
    #       supply=10_000,   # doses per year; allocated gap-proportionally
    #   )

    # Clark's API (option B — Supplies object):
    #   supply  = Supply(quantity=10_000, product=product)
    #   supplies = Supplies(supplies=[supply])
    #   prep = SuppliedPrep(
    #       name='lnc_combined',
    #       eligibilities=[fsw_elig, agyw_elig],
    #       coverages=[0.6, 0.3],
    #       supplies=supplies,
    #   )
    pass


# =============================================================================
# V5 — Kenya lenacapavir CE: price threshold sweep (aspirational)
#
# Research question: "At what cost per dose would Kenya be willing to pay
# for lenacapavir PrEP?" — CEMA collaborator request.
#
# Also: evaluate current rollout targeting (FSW first, then AGYW).
#
# Requires (not yet in stisim):
#   - Product with cost tracking
#   - Infections-averted analyzer (compared to no-PrEP baseline)
#   - Total cost accumulator (doses used × cost per dose)
#   - Kenya-calibrated HIV parameters (via hiv_kenya localization or pars)
# =============================================================================

def v5_kenya_ce_price_sweep():
    """
    NOTIMPLEMENTED — price threshold analysis for lenacapavir in Kenya.
    Sweep cost per person-year; compute ICER against oral PrEP baseline.
    """
    # Kenya WTP threshold: ~$1,000–2,000/DALY (rough estimate)
    costs_per_year = [50, 200, 500, 1_000, 5_000, 10_000, 40_000]

    # Baseline: no PrEP (or oral PrEP at current coverage)
    # baseline_sim = make_kenya_sim(interventions=[art, vmmc])

    results = []
    for cost in costs_per_year:
        # prep = sti.Prep(
        #     coverage=0.5,
        #     product=dict(name='lenacapavir', efficacy=0.999, cost=cost),
        #     eligibility=lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected,
        # )
        # sim = make_kenya_sim(interventions=[art, vmmc, prep])
        # sim.run()
        # results.append(sc.objdict(
        #     cost=cost,
        #     infections_averted=baseline_infections - sim.results.hiv.new_infections[:].sum(),
        #     total_cost=prep.total_cost,   # NOTIMPLEMENTED: cost accumulator
        #     icer=...,
        # ))
        pass

    # Plot: ICER vs. cost per year; mark Kenya WTP threshold
    # sc.plot(...)
    return results


def v5b_kenya_targeting_comparison():
    """
    NOTIMPLEMENTED — compare targeting strategies for a fixed lenacapavir supply.
    Same total doses; allocated to FSW-only vs. AGYW-only vs. FSW+AGYW vs. PBFW-only.
    """
    # target_strategies = {
    #     'fsw_only':   [dict(label='fsw',  eligibility=fsw_elig,  coverage=0.8)],
    #     'agyw_only':  [dict(label='agyw', eligibility=agyw_elig, coverage=0.4)],
    #     'fsw_agyw':   [dict(label='fsw',  eligibility=fsw_elig,  coverage=0.5),
    #                    dict(label='agyw', eligibility=agyw_elig, coverage=0.3)],
    #     'pbfw_only':  [dict(label='pbfw', eligibility=pbfw_elig, coverage=0.6)],
    # }
    # for label, groups in target_strategies.items():
    #     prep = sti.Prep(groups=groups, product='lenacapavir', supply=10_000)
    #     sim = make_kenya_sim(interventions=[art, vmmc, prep])
    #     sim.run()
    pass


# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    sc.heading('V1: simple FSW PrEP (runs today)')
    v1_simple()
    v1b_simple_with_eligibility()

    sc.heading('V2a: multi-group via multiple Prep instances (runs today)')
    v2a_multigroup_multiple_instances()

    sc.heading('V2b: multi-group via groups= API (runs today)')
    v2b_multigroup_groups_api()

    sc.heading('V3: Product-aware PrEP (simple API)')
    v3a_product_simple_api()

    sc.heading('V4–V5: API sketches — not yet implemented')
    # v3b_product_clark_api()
    # v4_supply_constrained()
    # v5_kenya_ce_price_sweep()
    # v5b_kenya_targeting_comparison()
