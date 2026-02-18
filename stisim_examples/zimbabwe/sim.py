"""
Make a Zimbabwe sim for testing and demonstration
"""
import numpy as np
import starsim as ss
import stisim as sti
import sciris as sc


# Zimbabwe-specific parameters
sim_pars = dict(
    start=1990,
    stop=2025,
    total_pop=9_980_999,
)

sti_pars = dict(
    hiv=dict(
        beta_m2f=0.035,
        eff_condom=0.95,
        rel_init_prev=1.0,
    )
)

nw_pars = dict(
    prop_f0=0.79,
    prop_m0=0.83,
    f1_conc=0.16,
    m1_conc=0.11,
    p_pair_form=0.58,
)


def make_custom_interventions(test_years=None):
    """
    Create any custom interventions for Zimbabwe. This example has testing interventions,
    other interventions can be added here. Note that the sim automatically adds
    ART and VMMC using the available data in the folder.

    Args:
        test_years (array): Years for testing coverage. Default: 1990-2040.

    Returns:
        list: List of intervention instances
    """
    if test_years is None:
        test_years = np.arange(1990, 2041)

    scaleup_end = min(2020, test_years[-1])
    scaleup_years = np.arange(test_years[0], scaleup_end + 1)
    n_scaleup = len(scaleup_years)
    n_future = len(test_years) - n_scaleup

    fsw_prob = np.concatenate([np.linspace(0, 0.75, n_scaleup), np.linspace(0.75, 0.85, n_future)])
    gp_prob = np.concatenate([np.linspace(0, 0.5, n_scaleup), np.linspace(0.5, 0.6, n_future)])
    low_cd4_prob = np.concatenate([np.linspace(0, 0.85, n_scaleup), np.linspace(0.85, 0.95, n_future)])

    def fsw_eligibility(sim):
        return sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    def other_eligibility(sim):
        return ~sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    def low_cd4_eligibility(sim):
        return (sim.diseases.hiv.cd4 < 200) & ~sim.diseases.hiv.diagnosed

    return [
        sti.HIVTest(years=test_years, test_prob_data=fsw_prob, name='fsw_testing', eligibility=fsw_eligibility, label='fsw_testing'),
        sti.HIVTest(years=test_years, test_prob_data=gp_prob, name='other_testing', eligibility=other_eligibility, label='other_testing'),
        sti.HIVTest(years=test_years, test_prob_data=low_cd4_prob, name='low_cd4_testing', eligibility=low_cd4_eligibility, label='low_cd4_testing'),
    ]


def make_sim(**kwargs):
    intvs = make_custom_interventions()
    user_intvs = sc.tolist(kwargs.pop('interventions', []))

    sim = sti.Sim(
        location='zimbabwe',
        diseases='hiv',
        data_path=sc.thispath(),
        sim_pars=sim_pars,
        nw_pars=nw_pars,
        sti_pars=sti_pars,
        interventions=intvs + user_intvs,
        **kwargs,
    )
    return sim
