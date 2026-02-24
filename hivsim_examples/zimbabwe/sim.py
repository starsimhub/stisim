"""
Zimbabwe HIV model with calibrated parameters and UNAIDS comparison data.

This example demonstrates how to build a data-driven HIV simulation:
- Initial prevalence from CSV by risk group/sex/SW status
- Condom use data by partnership type over time
- ART and VMMC coverage from national program data
- Testing interventions with FSW-targeted, general population, and low-CD4 strategies
- UNAIDS estimates attached to sim.data for plot overlay
"""
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti

datadir = sc.thispath()

# Calibrated parameters
sim_pars = dict(start=1990, stop=2025, total_pop=9_980_999)
sti_pars = dict(
    beta_m2f=0.035,
    eff_condom=0.95,
    rel_init_prev=1.0,
)
nw_pars = dict(
    prop_f0=0.79,
    prop_m0=0.83,
    f1_conc=0.16,
    m1_conc=0.11,
    p_pair_form=0.58,
)


def make_interventions():
    """
    Create Zimbabwe-specific HIV testing interventions.

    Returns FSW-targeted testing, general population testing, low-CD4 testing,
    and PrEP. ART and VMMC are created separately from coverage data.
    """
    scaleup_years = np.arange(1990, 2021)
    years = np.arange(1990, 2041)
    n_scaleup = len(scaleup_years)
    n_future = len(years) - n_scaleup

    fsw_prob = np.concatenate([np.linspace(0, 0.75, n_scaleup), np.linspace(0.75, 0.85, n_future)])
    gp_prob = np.concatenate([np.linspace(0, 0.5, n_scaleup), np.linspace(0.5, 0.6, n_future)])
    low_cd4_prob = np.concatenate([np.linspace(0, 0.85, n_scaleup), np.linspace(0.85, 0.95, n_future)])

    def fsw_eligibility(sim):
        return sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    def other_eligibility(sim):
        return ~sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    def low_cd4_eligibility(sim):
        return (sim.diseases.hiv.cd4 < 200) & ~sim.diseases.hiv.diagnosed

    fsw_testing = sti.HIVTest(years=years, test_prob_data=fsw_prob, name='fsw_testing', eligibility=fsw_eligibility)
    other_testing = sti.HIVTest(years=years, test_prob_data=gp_prob, name='other_testing', eligibility=other_eligibility)
    low_cd4_testing = sti.HIVTest(years=years, test_prob_data=low_cd4_prob, name='low_cd4_testing', eligibility=low_cd4_eligibility)
    prep = sti.Prep()

    return [fsw_testing, other_testing, low_cd4_testing, prep]


def make_sim(**kwargs):
    """
    Create a Zimbabwe HIV simulation with calibrated parameters.

    Loads demographic, behavioral, and epidemiological data from CSV files
    in this directory and creates a fully configured simulation.

    Args:
        **kwargs: Override any parameter passed to sti.Sim (e.g. n_agents, dur).
            If 'interventions' is provided, they are added to the defaults.

    Returns:
        sti.Sim: Configured simulation ready to run.
    """
    # Load data
    init_prev = pd.read_csv(datadir / 'init_prev_hiv.csv')
    condom_data = pd.read_csv(datadir / 'condom_use.csv')
    art_data = pd.read_csv(datadir / 'art_coverage.csv').set_index('year')
    vmmc_data = pd.read_csv(datadir / 'vmmc_coverage.csv').set_index('year')
    hiv_data = pd.read_csv(datadir / 'zimbabwe_hiv_data.csv')

    # Create modules
    hiv = sti.HIV(init_prev_data=init_prev, **sti_pars)
    network = sti.StructuredSexual(condom_data=condom_data, **nw_pars)
    art = sti.ART(coverage_data=art_data)
    vmmc = sti.VMMC(coverage_data=vmmc_data)

    # Combine interventions: testing + ART/VMMC + any user additions
    intvs = make_interventions() + [art, vmmc]
    user_intvs = sc.tolist(kwargs.pop('interventions', []))

    # Merge user sim_pars with defaults (user overrides take precedence)
    user_sim_pars = kwargs.pop('sim_pars', {})
    merged_sim_pars = sc.mergedicts(sim_pars, user_sim_pars)

    sim = sti.Sim(
        demographics='zimbabwe',
        datafolder=datadir,
        diseases=[hiv],
        networks=[network, ss.MaternalNet()],
        interventions=intvs + user_intvs,
        sim_pars=merged_sim_pars,
        data=hiv_data,
        **kwargs,
    )
    return sim
