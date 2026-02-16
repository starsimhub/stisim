"""
Zimbabwe-specific simulation configuration helpers.
"""
import numpy as np
import starsim as ss
import stisim as sti

__all__ = ['configure_sim_pars']


def make_hiv_testing(test_years=None):
    """
    Create HIV testing interventions for Zimbabwe.

    Args:
        test_years (array): Years for testing coverage. Default: 1990-2040.

    Returns:
        list: List of HIVTest intervention instances
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


def configure_sim_pars(loc_data, diseases, base_pars):
    """
    Configure simulation parameters from Zimbabwe location data.

    Creates fully configured disease, network, and intervention instances
    from the loaded location data. Uses Zimbabwe-specific network parameters
    and disease parameters from the hiv_zim country repo.

    Args:
        loc_data (sc.objdict): Data loaded from load_location_data()
        diseases (list): List of disease names
        base_pars (dict): Base parameters from user

    Returns:
        dict: Configured parameters including disease instances, networks,
              and interventions to merge with base_pars
    """
    configured = {}

    # Convert diseases to list if needed
    if isinstance(diseases, str):
        diseases = [diseases]

    if 'hiv' in diseases:
        # 1. Create HIV disease with init_prev data and Zimbabwe-specific parameters
        hiv_kwargs = dict(
            beta_m2f=0.035,
            eff_condom=0.95,
            rel_init_prev=1.0,
        )
        if 'hiv' in loc_data.diseases and 'init_prev' in loc_data.diseases.hiv:
            hiv_kwargs['init_prev_data'] = loc_data.diseases.hiv.init_prev
        configured['diseases'] = [sti.HIV(**hiv_kwargs)]

        # 2. Create network with condom data and Zimbabwe-specific parameters
        nw_kwargs = dict(
            prop_f0=0.79,
            prop_m0=0.83,
            f1_conc=0.16,
            m1_conc=0.11,
            p_pair_form=0.58,
        )
        if 'condom_use' in loc_data:
            nw_kwargs['condom_data'] = loc_data.condom_use
        configured['networks'] = [sti.StructuredSexual(**nw_kwargs), ss.MaternalNet()]

        # 3. Create interventions from data
        interventions = []

        # ART from coverage data (Zimbabwe uses absolute counts: n_art)
        if 'art_coverage' in loc_data:
            art_df = loc_data.art_coverage.set_index('year')
            interventions.append(sti.ART(coverage_data=art_df))

        # VMMC from coverage data (Zimbabwe uses absolute counts: n_vmmc)
        if 'vmmc_coverage' in loc_data:
            vmmc_df = loc_data.vmmc_coverage.set_index('year')
            interventions.append(sti.VMMC(coverage_data=vmmc_df))

        # PrEP
        interventions.append(sti.Prep())

        # HIV testing interventions
        interventions.extend(make_hiv_testing())

        configured['interventions'] = interventions

    # Default sim parameters for Zimbabwe
    if 'start' not in base_pars:
        configured['start'] = 1990
    if 'stop' not in base_pars and 'dur' not in base_pars:
        configured['stop'] = 2025
    if 'total_pop' not in base_pars:
        configured['total_pop'] = 9_980_999

    # Pass through calibration/comparison data for plotting overlay
    if 'hiv_data' in loc_data:
        configured['data'] = loc_data.hiv_data

    return configured
