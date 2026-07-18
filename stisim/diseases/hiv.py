"""
Define HIV disease module
"""

import numpy as np
import sciris as sc
import starsim as ss
import stisim as sti
from stisim.diseases.sti import BaseSTI, BaseSTIPars

__all__ = ['HIV', 'HIVPars']


# Default ARTMortalityTable-style dataset: annual HIV mortality hazard for
# agents on ART, stratified by ART duration bin x age bin x CD4 bin, separately
# for males/females and effective/non-suppressive ART. Transcribed from the
# EMOD-HIV campaign example in HIVSim_notes/art.md (art_implementation_notes.md
# section 7 has a discussion of what these numbers mean). Each table has shape
# (n_dur_bins, n_age_bins, n_cd4_bins) = (5, 4, 7), matching the bin edges below.
_ART_MORTALITY_AGE_BINS = np.array([25, 35, 45, 125])  # years; upper edges of <25/25-35/35-45/45+
_ART_MORTALITY_DUR_BINS = np.array([182, 365, 730, 1095, 45625])  # days since ART initiation
_ART_MORTALITY_CD4_BINS = np.array([0, 25, 74.5, 149.5, 274.5, 424.5, 624.5])  # CD4 count

_ART_MORTALITY_EFFECTIVE_MALE = np.array([
    [[0.2015, 0.2015, 0.1128, 0.0625, 0.0312, 0.0206, 0.0162],
     [0.2176, 0.2176, 0.1219, 0.0675, 0.0337, 0.0223, 0.0175],
     [0.2350, 0.2350, 0.1316, 0.0729, 0.0364, 0.0240, 0.0189],
     [0.2538, 0.2538, 0.1421, 0.0787, 0.0393, 0.0260, 0.0205]],
    [[0.0875, 0.0875, 0.0490, 0.0271, 0.0136, 0.0062, 0.0041],
     [0.0945, 0.0945, 0.0529, 0.0293, 0.0146, 0.0067, 0.0044],
     [0.1021, 0.1021, 0.0572, 0.0316, 0.0158, 0.0073, 0.0047],
     [0.1102, 0.1102, 0.0617, 0.0342, 0.0171, 0.0079, 0.0051]],
    [[0.0255, 0.0255, 0.0181, 0.0128, 0.0085, 0.0058, 0.0038],
     [0.0288, 0.0288, 0.0204, 0.0145, 0.0096, 0.0065, 0.0043],
     [0.0326, 0.0326, 0.0231, 0.0164, 0.0108, 0.0074, 0.0049],
     [0.0368, 0.0368, 0.0261, 0.0185, 0.0123, 0.0084, 0.0055]],
    [[0.0164, 0.0164, 0.0116, 0.0083, 0.0055, 0.0037, 0.0025],
     [0.0186, 0.0186, 0.0131, 0.0093, 0.0062, 0.0042, 0.0042],
     [0.0210, 0.0210, 0.0148, 0.0106, 0.0070, 0.0048, 0.0048],
     [0.0237, 0.0237, 0.0168, 0.0119, 0.0079, 0.0054, 0.0054]],
    [[0.0119, 0.0119, 0.0081, 0.0066, 0.0033, 0.0033, 0.0033],
     [0.0135, 0.0135, 0.0092, 0.0074, 0.0037, 0.0037, 0.0037],
     [0.0152, 0.0152, 0.0103, 0.0084, 0.0042, 0.0042, 0.0042],
     [0.0172, 0.0172, 0.0117, 0.0095, 0.0047, 0.0047, 0.0047]],
])

_ART_MORTALITY_NONSUPP_MALE = np.array([
    [[0.2015, 0.2015, 0.1128, 0.0625, 0.0312, 0.0206, 0.0162],
     [0.2176, 0.2176, 0.1219, 0.0675, 0.0337, 0.0223, 0.0175],
     [0.2350, 0.2350, 0.1316, 0.0729, 0.0364, 0.0240, 0.0189],
     [0.2538, 0.2538, 0.1421, 0.0787, 0.0393, 0.0260, 0.0205]],
    [[0.1715, 0.1715, 0.5600, 0.3100, 0.1550, 0.0713, 0.0465],
     [0.1852, 0.1852, 0.6048, 0.3348, 0.1674, 0.0770, 0.0502],
     [0.2000, 0.2000, 0.6532, 0.3616, 0.1808, 0.0832, 0.0542],
     [0.2160, 0.2160, 0.7054, 0.3905, 0.1953, 0.0898, 0.0586]],
    [[0.0532, 0.0532, 0.0362, 0.0293, 0.0171, 0.0116, 0.0095],
     [0.0601, 0.0601, 0.0409, 0.0331, 0.0193, 0.0131, 0.0107],
     [0.0679, 0.0679, 0.0462, 0.0374, 0.0218, 0.0148, 0.0121],
     [0.0768, 0.0768, 0.0522, 0.0422, 0.0246, 0.0168, 0.0137]],
    [[0.0335, 0.0335, 0.0228, 0.0184, 0.0108, 0.0073, 0.0060],
     [0.0379, 0.0379, 0.0258, 0.0208, 0.0122, 0.0083, 0.0068],
     [0.0428, 0.0428, 0.0291, 0.0235, 0.0137, 0.0094, 0.0076],
     [0.0484, 0.0484, 0.0329, 0.0266, 0.0155, 0.0106, 0.0086]],
    [[0.0234, 0.0234, 0.0159, 0.0129, 0.0091, 0.0069, 0.0064],
     [0.0265, 0.0265, 0.0180, 0.0145, 0.0103, 0.0077, 0.0073],
     [0.0299, 0.0299, 0.0203, 0.0164, 0.0116, 0.0088, 0.0082],
     [0.0338, 0.0338, 0.0230, 0.0186, 0.0131, 0.0099, 0.0093]],
])

_ART_MORTALITY_EFFECTIVE_FEMALE = np.array([
    [[0.2015, 0.2015, 0.0993, 0.0518, 0.0259, 0.0171, 0.0135],
     [0.2156, 0.2156, 0.1062, 0.0554, 0.0277, 0.0183, 0.0144],
     [0.2307, 0.2307, 0.1137, 0.0593, 0.0296, 0.0196, 0.0154],
     [0.2468, 0.2468, 0.1216, 0.0634, 0.0317, 0.0209, 0.0165]],
    [[0.0875, 0.0875, 0.0431, 0.0225, 0.0112, 0.0052, 0.0034],
     [0.0936, 0.0936, 0.0461, 0.0241, 0.0120, 0.0055, 0.0036],
     [0.1002, 0.1002, 0.0494, 0.0257, 0.0129, 0.0059, 0.0039],
     [0.1072, 0.1072, 0.0528, 0.0276, 0.0138, 0.0063, 0.0041]],
    [[0.0241, 0.0241, 0.0166, 0.0135, 0.0067, 0.0044, 0.0044],
     [0.0262, 0.0262, 0.0181, 0.0147, 0.0073, 0.0048, 0.0048],
     [0.0286, 0.0286, 0.0197, 0.0160, 0.0080, 0.0052, 0.0052],
     [0.0312, 0.0312, 0.0215, 0.0175, 0.0087, 0.0057, 0.0057]],
    [[0.0149, 0.0149, 0.0103, 0.0084, 0.0042, 0.0042, 0.0042],
     [0.0163, 0.0163, 0.0112, 0.0091, 0.0046, 0.0046, 0.0046],
     [0.0177, 0.0177, 0.0122, 0.0099, 0.0050, 0.0050, 0.0050],
     [0.0193, 0.0193, 0.0133, 0.0108, 0.0054, 0.0054, 0.0054]],
    [[0.0084, 0.0084, 0.0057, 0.0046, 0.0023, 0.0023, 0.0023],
     [0.0092, 0.0092, 0.0062, 0.0051, 0.0025, 0.0025, 0.0025],
     [0.0100, 0.0100, 0.0068, 0.0055, 0.0028, 0.0028, 0.0028],
     [0.0109, 0.0109, 0.0074, 0.0060, 0.0030, 0.0030, 0.0030]],
])

_ART_MORTALITY_NONSUPP_FEMALE = np.array([
    [[0.2015, 0.2015, 0.1128, 0.0625, 0.0312, 0.0206, 0.0162],
     [0.2176, 0.2176, 0.1219, 0.0675, 0.0337, 0.0223, 0.0175],
     [0.2350, 0.2350, 0.1316, 0.0729, 0.0364, 0.0240, 0.0189],
     [0.2538, 0.2538, 0.1421, 0.0787, 0.0393, 0.0260, 0.0205]],
    [[0.1837, 0.1837, 0.0845, 0.0441, 0.0220, 0.0101, 0.0066],
     [0.1965, 0.1965, 0.0904, 0.0472, 0.0236, 0.0108, 0.0071],
     [0.2103, 0.2103, 0.0967, 0.0505, 0.0252, 0.0116, 0.0076],
     [0.2250, 0.2250, 0.1035, 0.0540, 0.0270, 0.0124, 0.0081]],
    [[0.0461, 0.0461, 0.0318, 0.0258, 0.0151, 0.0103, 0.0084],
     [0.0502, 0.0502, 0.0346, 0.0281, 0.0164, 0.0113, 0.0091],
     [0.0547, 0.0547, 0.0378, 0.0306, 0.0179, 0.0123, 0.0100],
     [0.0596, 0.0596, 0.0412, 0.0334, 0.0195, 0.0134, 0.0109]],
    [[0.0286, 0.0286, 0.0197, 0.0160, 0.0113, 0.0086, 0.0080],
     [0.0311, 0.0311, 0.0215, 0.0174, 0.0124, 0.0094, 0.0087],
     [0.0339, 0.0339, 0.0234, 0.0190, 0.0135, 0.0102, 0.0095],
     [0.0370, 0.0370, 0.0255, 0.0207, 0.0147, 0.0111, 0.0104]],
    [[0.0161, 0.0161, 0.0110, 0.0089, 0.0063, 0.0047, 0.0044],
     [0.0176, 0.0176, 0.0119, 0.0097, 0.0068, 0.0052, 0.0048],
     [0.0192, 0.0192, 0.0130, 0.0105, 0.0075, 0.0056, 0.0053],
     [0.0209, 0.0209, 0.0142, 0.0115, 0.0081, 0.0061, 0.0057]],
])

_DEFAULT_ART_MORTALITY_TABLE = dict(
    age_bins=_ART_MORTALITY_AGE_BINS,
    dur_bins=_ART_MORTALITY_DUR_BINS,
    cd4_bins=_ART_MORTALITY_CD4_BINS,
    effective=dict(m=_ART_MORTALITY_EFFECTIVE_MALE, f=_ART_MORTALITY_EFFECTIVE_FEMALE),
    nonsuppressive=dict(m=_ART_MORTALITY_NONSUPP_MALE, f=_ART_MORTALITY_NONSUPP_FEMALE),
)


class HIVPars(BaseSTIPars):
    """
    Parameters for the HIV disease module.

    Holds natural history parameters (CD4 dynamics, acute/latent/falling stage
    durations), transmission rates, ART treatment effects, and care-seeking
    behavior configuration.
    """
    def __init__(self, **kwargs):
        super().__init__()

        # Natural history without treatment
        self.cd4_start = ss.normal(loc=800, scale=50)
        self.cd4_latent = ss.normal(loc=500, scale=50)
        self.dur_acute = ss.lognorm_ex(ss.months(1.7), ss.months(1))  # Duration of acute HIV infection (Bellan 2015: central estimate 1.7 mo)
        self.dur_latent = ss.lognorm_ex(ss.years(10), ss.years(3))  # Duration of latent, untreated HIV infection
        self.dur_falling = ss.lognorm_ex(ss.years(3), ss.years(1))  # Duration of late-stage HIV when CD4 counts fall
        self.p_hiv_death = ss.bernoulli(p=0)  # Death from HIV-related complications; set by make_p_hiv_death()
        self.include_aids_deaths = True

        # Transmission
        self.beta = 0  # Placeholder, replaced by network-specific betas
        self.beta_m2f = 0.05
        self.rel_beta_f2m = 0.5
        self.beta_m2c = ss.permonth(0.025)  # Prenatal MTCT: ~20% over pregnancy
        self.beta_breastfeed = ss.permonth(0.005)  # Postnatal MTCT via breastfeeding (~14% over 12 months without ART)
        self.rel_trans_acute = ss.normal(loc=5.3, scale=0.5)  # Increase transmissibility during acute HIV infection (Bellan 2015: RH=5.3, EHM≈7.3)
        self.rel_trans_falling = ss.normal(loc=8, scale=0.5)  # Increase transmissibility during late HIV infection
        self.rel_sus_age = None  # Age/sex-stratified rel_sus multipliers.
        # List of (age_lo, age_hi, sex, multiplier) tuples; sex is 'f', 'm', or None for both.
        # Multiplied into rel_sus each timestep alongside other contributors (PrEP, connectors, etc.),
        # so values are ceteris-paribus differences, not absolute susceptibilities.
        # Example: [(15, 25, 'f', 1.7), (25, 50, 'f', 1.0), (15, 50, 'm', 1.0)]
        self.eff_condom = 0.9

        # Initialization
        self.init_prev = ss.bernoulli(p=0.05)
        self.rel_init_prev = 1
        self.init_diagnosed = ss.bernoulli(p=0)
        self.dist_ti_init_infected = ss.uniform(low=-10 * 12, high=-5)
        # dist_ti_init_infected=ss.constant(0),  # Experimented with negative values, but safer to use 0

        # Care seeking
        self.care_seeking = ss.normal(loc=1, scale=0.5)  # Distribution of relative care-seeking behavior (wide spread for low-care-seeking extremes)
        self.maternal_care_scale = 2  # Factor for scaling up care-seeking behavior during pregnancy

        # Treatment effects
        self.art_cd4_growth = 0.1  # Unitless parameter defining how CD4 reconstitutes after starting ART - used in a logistic growth function
        self.effective_art_efficacy = 0.96  # Transmission efficacy of effective (virally-suppressive) ART
        self.nonsupp_art_efficacy = 0.35  # Transmission efficacy of non-suppressive ART
        self.time_to_art_efficacy = ss.months(6)  # Time to reach full ART efficacy (in months) - linear increase in efficacy
        self.p_effective_art = ss.bernoulli(p=1.0)  # Probability that a newly-initiated agent achieves viral suppression (vs. non-suppressive ART)
        self.art_cd4_pars = dict(cd4_max=1000, cd4_healthy=500)
        self.dur_on_art = ss.lognorm_ex(ss.years(3), ss.years(1.5))  # Base ART duration (scaled by rel_dur_on_art and trend)
        self.rel_dur_on_art = 1.0  # Calibratable scalar that scales ART duration
        self.dur_on_art_trend = None  # Optional time-varying scale, e.g. sc.objdict(years=np.array([2004, 2015, 2030]), vals=np.array([0.5, 1.0, 1.5]))
        self.dur_post_art = ss.normal()  # Note defined in years!
        self.dur_post_art_scale_factor = 0.1

        # ARTMortalityTable-style differential ART mortality (art.md step 4). Off by
        # default: with use_art_mortality_table=False, on-ART agents keep the original
        # zero-HIV-mortality-while-on-ART behavior (make_p_hiv_death only applies to
        # off-ART agents). Setting it True looks up a nonzero, age/sex/adherence/
        # duration-since-initiation-stratified hazard for on-ART agents instead (see
        # HIV.get_art_mortality_hazard). art_mortality_table ships with the example
        # values from HIVSim_notes/art.md's EMOD campaign excerpt; pass your own dict
        # of the same shape (age_bins/dur_bins/cd4_bins + effective/nonsuppressive
        # tables keyed by 'm'/'f') to use different data.
        self.use_art_mortality_table = False
        self.art_mortality_table = sc.dcp(_DEFAULT_ART_MORTALITY_TABLE)

        self.update(kwargs)
        return


class HIV(BaseSTI):
    """
    HIV disease module.

    Models HIV infection with CD4-driven natural history (acute, latent, and
    falling stages), ART treatment with CD4 reconstitution, diagnosis tracking,
    and AIDS-related mortality. Supports mother-to-child transmission when a
    pregnancy demographic module is present.

    Args:
        pars (dict): Override default parameters from ``HIVPars``.
        init_prev_data: Optional initial prevalence data by age/sex/risk group.
        **kwargs: Additional parameters passed to ``update_pars``.
    """

    def __init__(self, pars=None, init_prev_data=None, **kwargs):
        super().__init__()

        # Parameters
        default_pars = HIVPars()
        self.define_pars(**default_pars)
        self.update_pars(pars, **kwargs)

        # Set initial prevalence
        self.init_prev_data = init_prev_data

        # States
        self.define_states(
            # Natural history
            ss.FloatArr('ti_exposed'),
            ss.FloatArr('ti_acute'),
            ss.BoolState('acute'),
            ss.FloatArr('ti_latent'),
            ss.BoolState('latent'),
            ss.FloatArr('ti_falling'),
            ss.BoolState('falling'),
            ss.BoolState('post_art'),  # After stopping ART, CD4 falls linearly until death
            ss.FloatArr('ti_zero'),  # Time of zero CD4 count - generally corresponds to AIDS death
            ss.FloatArr('ti_dead'),  # Time of HIV/AIDS death

            # Care and treatment states
            ss.FloatArr('baseline_care_seeking'),
            ss.FloatArr('care_seeking'),
            ss.BoolState('art_naive', default=True),        # Never initiated ART
            ss.BoolState('on_art'),                          # Currently on ART (effective or non-suppressive)
            ss.BoolState('on_effective_art'),                # Currently on ART and virally suppressed
            ss.BoolState('on_nonsuppressive_art'),           # Currently on ART but not virally suppressed
            ss.BoolState('art_discontinued'),                # Previously on ART, currently off
            ss.FloatArr('ti_art'),
            ss.FloatArr('ti_stop_art'),

            # CD4 states
            ss.FloatArr('cd4'),             # Current CD4 count
            ss.FloatArr('cd4_start'),       # Initial CD4 count for each agent before an infection
            ss.FloatArr('cd4_preart'),      # CD4 immediately before initiating ART
            ss.FloatArr('cd4_latent'),      # CD4 count during latent infection
            ss.FloatArr('cd4_nadir'),       # Lowest CD4
            ss.FloatArr('cd4_potential'),   # Potential CD4 count if continually treated
            ss.FloatArr('cd4_postart'),     # CD4 after stopping ART

            # Knowledge of HIV status
            ss.BoolState('diagnosed'),
            ss.FloatArr('ti_diagnosed'),
        )

        return

    @property
    def include_mtct(self): return 'pregnancy' in self.sim.demographics

    def plot(self):
        """ Plot key HIV results """
        import matplotlib.pyplot as plt
        res = self.results
        keys = ['new_infections', 'new_deaths', 'prevalence', 'prevalence_15_49', 'new_diagnoses', 'p_on_art']
        keys = [k for k in keys if k in res] # Only plot keys that exist
        with sc.options.with_style('fancy'):
            fig, axs = sc.getrowscols(len(keys), make=True)
            for ax, k in zip(axs.flatten(), keys):
                ax.plot(res.timevec, res[k], lw=2, alpha=0.8)
                ax.set_title(res[k].label)
                ax.set_xlabel('Year')
            sc.figlayout()
        return ss.return_fig(fig)

    def init_results(self):
        """
        Initialize results
        """
        super().init_results()
        results = [
            ss.Result('new_deaths', dtype=int, label='Deaths'),
            ss.Result('cum_deaths', dtype=int, label='Cumulative deaths', auto_plot=False),
            ss.Result('new_diagnoses', dtype=int, label='Diagnoses'),
            ss.Result('cum_diagnoses', dtype=int, label='Cumulative diagnoses', auto_plot=False),
            ss.Result('new_agents_on_art', dtype=int, label='New treated', auto_plot=False),
            ss.Result('prevalence_15_49', dtype=float, label='Prevalence 15-49', scale=False),
            ss.Result('prevalence_sw', dtype=float, label='FSW prevalence', scale=False, auto_plot=False),
            ss.Result('new_infections_sw', dtype=int, label='New infections - FSW', auto_plot=False),
            ss.Result('new_infections_not_sw', dtype=int, label='New infections - Other F', auto_plot=False),
            ss.Result('prevalence_client', dtype=float, label='Client prevalence', scale=False, auto_plot=False),
            ss.Result('new_infections_client', dtype=int, label='New infections - Clients', auto_plot=False),
            ss.Result('new_infections_not_client', dtype=int, label='New infections - Other M', auto_plot=False),
            ss.Result('p_on_art', dtype=float, label='Proportion on ART', scale=False),
        ]

        if self.include_mtct:
            results += [
                ss.Result('n_on_art_pregnant', dtype=int, auto_plot=False),
                ss.Result('p_diagnosed_pregnant', dtype=float, label='Proportion of HIV+ pregnant women diagnosed', scale=False, auto_plot=False),
            ]

        # Add FSW and clients to results:
        if 'structuredsexual' in self.sim.networks.keys():
            for risk_group in range(self.sim.networks.structuredsexual.pars.n_risk_groups):
                for sex in ['female', 'male']:
                    results += [
                        ss.Result('prevalence_risk_group_' + str(risk_group) + '_' + sex, scale=False, auto_plot=False),
                        ss.Result('new_infections_risk_group_' + str(risk_group) + '_' + sex, dtype=int, auto_plot=False),
                    ]

        self.define_results(*results)

        return

    def init_post(self):
        """ Set states """
        ss.Module.init_post(self)  # Skip the disease init_post() since we create infections in a different way
        # Set initial CD4
        self.init_cd4()
        self.init_care_seeking()

        # Make initial cases, some of which may have occured prior to the sim start
        if self.init_prev_data is not None:
            p_init_infection = self.make_init_prev()
            self.pars.init_prev.set(p_init_infection)  # Set the initial prevalence function
        initial_cases = self.pars.init_prev.filter()
        ti_init_cases = self.pars.dist_ti_init_infected.rvs(initial_cases).astype(int)
        self.set_prognoses(initial_cases, ti=ti_init_cases)
        initial_cases_diagnosed = self.pars.init_diagnosed.filter(initial_cases)
        self.diagnosed[initial_cases_diagnosed] = True
        self.ti_diagnosed[initial_cases_diagnosed] = 0
        # Schedule ART start at ti=0 (no delay) so an ART intervention picks them up immediately
        self.ti_art[initial_cases_diagnosed] = 0
        return

    # CD4 functions
    def acute_decline(self, uids):
        """ Acute decline in CD4 """
        acute_start = self.ti_acute[uids]   # Time of infection
        acute_end = self.ti_latent[uids]    # Time to latent infection
        acute_dur = acute_end - acute_start  # Total time in acute phase, rounded to timesteps
        cd4_start = self.cd4_start[uids]
        cd4_end = self.cd4_latent[uids]
        per_timestep_decline = sc.safedivide(cd4_start-cd4_end, acute_dur)
        cd4 = self.cd4[uids] - per_timestep_decline
        return cd4

    def falling_decline(self, uids):
        """ Decline in CD4 during late-stage infection, when counts are falling """
        falling_start = self.ti_falling[uids]
        falling_end = self.ti_zero[uids]
        falling_dur = falling_end - falling_start
        time_falling = self.ti - self.ti_falling[uids]
        cd4_start = self.cd4_latent[uids]
        cd4_end = 1  # To avoid divide by zero problems
        per_timestep_decline = sc.safedivide(cd4_start-cd4_end, falling_dur)
        cd4 = np.maximum(cd4_end, self.cd4[uids] - per_timestep_decline)
        # cd4 = np.maximum(0, cd4_start - per_timestep_decline*time_falling)
        return cd4

    def post_art_decline(self, uids):
        """
        Decline in CD4 after going off treatment
        This implementation has the possibly-undesirable feature that a person
        who goes on ART for a year and then off again might have a slightly shorter
        lifespan than if they'd never started treatment.
        """
        cd4_end = 1  # To avoid divide by zero problems

        # Death immediately
        zero_now_inds = (self.ti_zero[uids] <= self.ti).nonzero()[-1]
        zero_later_inds = (self.ti_zero[uids] > self.ti).nonzero()[-1]
        cd4 = self.cd4[uids]
        cd4[zero_now_inds] = cd4_end

        # Death later
        zero_later_uids = uids.remove(uids[zero_now_inds])
        ti_zero = self.ti_zero[zero_later_uids]
        ti_stop_art = self.ti_stop_art[zero_later_uids]
        post_art_dur = ti_zero - ti_stop_art
        time_post_art = self.ti - ti_stop_art
        cd4_start = self.cd4_postart[zero_later_uids]
        if (post_art_dur <= 0).any():
            post_art_dur[post_art_dur <= 0] = 1
            # error_msg = 'Post-ART duration is negative'
            # raise ValueError(error_msg)
        per_timestep_decline = (cd4_start-cd4_end)/post_art_dur
        cd4[zero_later_inds] = np.maximum(cd4_end, cd4_start - per_timestep_decline*time_post_art)
        return cd4

    def cd4_increase(self, uids):
        """
        Increase CD4 counts for people who are receiving treatment.
        Growth curves are calculated to match EMODs CD4 reconstitution equation for people who initiate treatment
        with a CD4 count of 50 (https://docs.idmod.org/projects/emod-hiv/en/latest/hiv-model-healthcare-systems.html)
        However, here we use a logistic growth function and assume that ART CD4 count depends on CD4 at initiation.
        Sources:

            - https://i-base.info/guides/starting/cd4-increase
            - https://www.sciencedirect.com/science/article/pii/S1876034117302022
            - https://bmcinfectdis.biomedcentral.com/articles/10.1186/1471-2334-8-20
        """
        # Calculate time on ART and CD4 prior to starting
        ti_art = self.ti_art[uids]
        cd4_preart = self.cd4_preart[uids]
        dur_art = self.ti - ti_art

        # Extract growth parameters
        growth_rate = self.pars.art_cd4_growth
        cd4_total_gain = self.cd4_potential[uids] - self.cd4_preart[uids]
        cd4_now = 2*cd4_total_gain/(1+np.exp(-dur_art*growth_rate))-cd4_total_gain+cd4_preart  # Concave logistic

        return cd4_now

    @property
    def symptomatic(self):
        return self.infectious

    @property
    def aids(self):
        return self.cd4 < 200

    def make_p_hiv_death(self, uids=None):
        """
        Calculate per-timestep HIV death probability based on current CD4 count.

        Uses CD4-stratified annual mortality rates, digitized into bins. Rates are
        converted from per-year to per-timestep probabilities.
        """
        cd4_bins = np.array([1000, 500, 350, 200, 50, 0])
        p_hiv_death = ss.peryear(np.array([0.003, 0.003, 0.005, 0.01, 0.05, 0.300])).to_prob(self.dt)
        return p_hiv_death[np.digitize(self.cd4[uids], cd4_bins)]

    def get_art_mortality_hazard(self, uids):
        """
        ARTMortalityTable-style per-timestep HIV death probability for agents
        currently on ART (`pars.use_art_mortality_table=True` only).

        Looks up an annual mortality hazard stratified by age, sex, ART
        adherence category (effective vs. non-suppressive), and duration
        since ART initiation, from `pars.art_mortality_table`. Unlike EMOD's
        ARTMortalityTable, which freezes age/CD4 at the moment of initiation
        (a simplification specific to EMOD's one-time survival-draw design),
        this looks up the agent's CURRENT CD4 every time it's called, since
        STIsim already recomputes CD4 continuously for other purposes (see
        art_implementation_notes.md section 9 for the rationale).
        """
        table = self.pars.art_mortality_table
        age_bins = table['age_bins']
        dur_bins = table['dur_bins']
        cd4_bins = table['cd4_bins']

        age_idx = np.clip(np.digitize(self.sim.people.age[uids], age_bins), 0, len(age_bins) - 1)
        dur_days = (self.ti - self.ti_art[uids]) * self.dt.days
        dur_idx = np.clip(np.digitize(dur_days, dur_bins), 0, len(dur_bins) - 1)
        cd4_idx = np.clip(np.digitize(self.cd4[uids], cd4_bins), 0, len(cd4_bins) - 1)

        male = self.sim.people.male[uids]
        effective = self.on_effective_art[uids]

        annual_rate = np.empty(len(uids))
        for adherence_key, adherence_mask in (('effective', effective), ('nonsuppressive', ~effective)):
            for sex_key, sex_mask in (('m', male), ('f', ~male)):
                mask = adherence_mask & sex_mask
                if not mask.any():
                    continue
                mtable = table[adherence_key][sex_key]  # shape (n_dur, n_age, n_cd4)
                annual_rate[mask] = mtable[dur_idx[mask], age_idx[mask], cd4_idx[mask]]

        return ss.peryear(annual_rate).to_prob(self.dt)

    @staticmethod
    def _interpolate(vals: list, t):
        vals = sorted(vals, key=lambda x: x[0])  # Make sure values are sorted
        assert len({x[0] for x in vals}) == len(vals)  # Make sure time points are unique
        return np.interp(t, [x[0] for x in vals], [x[1] for x in vals], left=vals[0][1], right=vals[-1][1])

    def init_cd4(self):
        """
        Set CD4 counts
        """
        uids = ss.uids(np.isnan(self.cd4_start.raw).nonzero()[-1])
        self.cd4_start[uids] = self.pars.cd4_start.rvs(uids)
        self.cd4_nadir[uids] = sc.dcp(self.cd4_start[uids])
        return

    def init_care_seeking(self):
        """
        Set care seeking behavior
        """
        uids = ss.uids(self.care_seeking.isnan)
        self.care_seeking[uids] = self.pars.care_seeking.rvs(uids)
        self.care_seeking[uids] = np.maximum(self.care_seeking[uids], 0.01)  # Floor to prevent division by zero
        self.baseline_care_seeking[uids] = sc.dcp(self.care_seeking[uids])  # Copy it so pregnancy can modify it
        return

    def step_state(self):
        """
        Carry out autonomous updates at the start of the timestep (prior to transmission)
        """
        ti = self.ti

        # Set initial CD4 counts for new agents:
        self.init_cd4()

        # Handle care seeking behavior. First, initialize, then adjust depending on pregnancy:
        # increase care-seeking for pregnant women and decrease again after postpartum.
        # This makes it much less likely that pregnant women will stop treatment
        self.init_care_seeking()
        if self.include_mtct:
            pregnant = self.sim.demographics.pregnancy.pregnant
            self.care_seeking[pregnant] = self.baseline_care_seeking[pregnant] * self.pars.maternal_care_scale
            self.care_seeking[~pregnant] = self.baseline_care_seeking[~pregnant]

        # Adjust CD4 counts for people receiving treatment - logarithmic increase
        if self.on_art.any():
            art_uids = self.on_art.uids
            self.cd4[art_uids] = self.cd4_increase(art_uids)

        # Adjust CD4 counts for people who have gone off treatment - linear decline
        if self.art_discontinued.any():
            off_art_uids = self.art_discontinued.uids
            self.cd4[off_art_uids] = self.post_art_decline(off_art_uids)

        # Update states for people who have never been on ART (ART removes these)
        # Acute & not on ART
        latent = self.acute & (self.ti_latent <= ti)
        falling = self.latent & (self.ti_falling <= ti)
        self.acute[latent] = False
        self.latent[latent] = True
        self.latent[falling] = False
        self.falling[falling] = True

        # Update CD4 counts
        self.cd4[self.acute.uids] = self.acute_decline(self.acute.uids)
        untreated_latent = self.latent
        self.cd4[untreated_latent.uids] = self.cd4_latent[untreated_latent.uids]
        untreated_falling = self.falling
        if untreated_falling.any():
            self.cd4[untreated_falling.uids] = self.falling_decline(untreated_falling.uids)

        # Update CD4 nadir for anyone not on treatment
        untreated = self.infected & ~self.on_art
        self.cd4_nadir[untreated] = np.minimum(self.cd4_nadir[untreated], self.cd4[untreated])

        # Update transmission
        self.update_transmission()

        # Check CD4
        if np.isnan(self.cd4[self.infected]).any():
            errormsg = 'Invalid entry for CD4'
            raise ValueError(errormsg)

        # Update deaths. We capture deaths from AIDS (i.e., when CD4 count drops to ~0) as well as deaths from
        # serious HIV-related illnesses, which can occur throughout HIV.
        off_art = (self.infected & ~self.on_art).uids
        p_death = self.make_p_hiv_death(uids=off_art)
        self.pars.p_hiv_death.set(0)
        self.pars.p_hiv_death.set(p_death)  # Set the death probability function
        hiv_deaths = self.pars.p_hiv_death.filter(off_art)

        # ARTMortalityTable-style hazard: opt-in (pars.use_art_mortality_table), gives
        # on-ART agents a nonzero, age/sex/adherence/duration-stratified mortality risk
        # instead of the default (off-ART-only hazard, i.e. zero mortality while on ART)
        if self.pars.use_art_mortality_table:
            on_art = (self.infected & self.on_art).uids
            if len(on_art):
                p_death_on_art = self.get_art_mortality_hazard(on_art)
                self.pars.p_hiv_death.set(0)
                self.pars.p_hiv_death.set(p_death_on_art)
                on_art_deaths = self.pars.p_hiv_death.filter(on_art)
                if len(on_art_deaths):
                    hiv_deaths = ss.uids(np.concatenate([hiv_deaths, on_art_deaths]))

        if len(hiv_deaths):
            self.ti_dead[hiv_deaths] = ti
            self.sim.people.request_death(hiv_deaths)
        if self.pars.include_aids_deaths:
            aids_deaths = (self.ti_zero <= ti).uids
            if len(aids_deaths):
                self.ti_dead[aids_deaths] = ti
                self.sim.people.request_death(aids_deaths)
        return

    def step_die(self, uids):
        """ Clear all states for dead agents """

        # Reset boolean states
        self.susceptible[uids] = False
        self.infected[uids] = False
        self.acute[uids] = False
        self.latent[uids] = False
        self.falling[uids] = False
        self.post_art[uids] = False
        self.art_naive[uids] = False
        self.on_art[uids] = False
        self.on_effective_art[uids] = False
        self.on_nonsuppressive_art[uids] = False
        self.art_discontinued[uids] = False
        self.diagnosed[uids] = False

        # Clear time states except for ti_dead
        self.ti_infected[uids] = np.nan
        self.ti_acute[uids] = np.nan
        self.ti_latent[uids] = np.nan
        self.ti_falling[uids] = np.nan
        self.ti_zero[uids] = np.nan
        self.ti_art[uids] = np.nan
        self.ti_stop_art[uids] = np.nan
        self.ti_diagnosed[uids] = np.nan

        # Clear CD4 states
        self.cd4[uids] = np.nan
        self.cd4_start[uids] = np.nan
        self.cd4_preart[uids] = np.nan
        self.cd4_latent[uids] = np.nan
        self.cd4_nadir[uids] = np.nan
        self.cd4_potential[uids] = np.nan
        self.cd4_postart[uids] = np.nan

        # Clear all other states
        self.rel_sus[uids] = np.nan
        self.rel_trans[uids] = np.nan
        self.baseline_care_seeking[uids] = np.nan
        self.care_seeking[uids] = np.nan

        return

    def update_transmission(self):
        """
        Update rel_trans and rel_sus for all agents. These are reset on each timestep then adjusted depending on states.
        Adjustments are made throughout different modules:

           - rel_trans for acute and late-stage untreated infection are adjusted below
           - rel_trans for all people on treatment (including pregnant women) below
           - rel_sus for unborn babies of pregnant WLHIV receiving treatment is adjusted in the ART intervention
        """
        ti = self.ti

        # Reset susceptibility and infectiousness
        self.rel_sus[:] = 1
        self.rel_trans[:] = 1

        # Age/sex-stratified susceptibility. Multiplied in alongside other contributors
        # (PrEP, connectors, etc.) so values are ceteris-paribus, not absolute.
        if self.pars.rel_sus_age is not None:
            ppl = self.sim.people
            for age_lo, age_hi, sex, mult in self.pars.rel_sus_age:
                if sex == 'f':   sex_mask = ppl.female
                elif sex == 'm': sex_mask = ppl.male
                else:            sex_mask = ppl.alive
                in_bin = sex_mask & (ppl.age >= age_lo) & (ppl.age < age_hi)
                self.rel_sus[in_bin] *= mult

        # Update rel_trans to account for acute and late-stage infection
        self.rel_trans[self.acute] *= self.pars.rel_trans_acute.rvs(self.acute.uids)
        self.rel_trans[self.aids] *= self.pars.rel_trans_falling.rvs(self.aids.uids)

        # Update transmission for agents on ART
        # When agents start ART, determine the reduction of transmission (linearly decreasing over 6 months).
        # Efficacy depends on whether the agent is virally suppressed (effective ART) or not (non-suppressive ART).
        if self.on_art.any():
            time_to_full_eff = self.pars.time_to_art_efficacy
            art_uids = self.on_art.uids
            full_eff = np.where(self.on_effective_art[art_uids], self.pars.effective_art_efficacy, self.pars.nonsupp_art_efficacy)
            timesteps_on_art = ti - self.ti_art[art_uids]
            new_on_art = timesteps_on_art < time_to_full_eff/self.dt
            efficacy_to_date = full_eff.copy()
            efficacy_to_date[new_on_art] = timesteps_on_art[new_on_art]*full_eff[new_on_art]/time_to_full_eff.value
            self.rel_trans[art_uids] *= 1 - efficacy_to_date

        return

    def update_results(self):
        """
        Update results at each time step
        """
        super().update_results()
        ti = self.ti

        # Recalculate prevalence so it's for the whole population - the STI module calculates it for adults
        self.results['prevalence'][ti] = sum(self.infected) / len(self.infected)

        self.results['new_deaths'][ti] = np.count_nonzero(self.ti_dead == ti)
        self.results['cum_deaths'][ti] = np.sum(self.results['new_deaths'][:ti + 1])
        self.results['new_diagnoses'][ti] = np.count_nonzero(self.ti_diagnosed == ti)
        self.results['cum_diagnoses'][ti] = np.sum(self.results['new_diagnoses'][:ti + 1])
        self.results['new_agents_on_art'][ti] = sum((self.ti_art == ti) & self.on_art)  # gate on on_art; ti_art may also hold a scheduled future start
        if self.include_mtct:
            pregnant = self.sim.people.pregnancy.pregnant
            self.results['n_on_art_pregnant'][ti] = np.count_nonzero(self.on_art & pregnant)
            n_infected_pregnant = np.count_nonzero(self.infected & pregnant)
            n_diagnosed_pregnant = np.count_nonzero(self.diagnosed & pregnant & self.infected)
            self.results['p_diagnosed_pregnant'][ti] = sc.safedivide(n_diagnosed_pregnant, n_infected_pregnant)
        self.results['p_on_art'][ti] = sc.safedivide(self.results['n_on_art'][ti], self.results['n_infected'][ti])

        # Subset by age group:
        infected_15_19 = self.infected[(self.sim.people.age >= 15) & (self.sim.people.age < 50)]
        if len(infected_15_19): # Avoid divide-by-zero error
            self.results['prevalence_15_49'][ti] = sum(infected_15_19) / len(infected_15_19)

        # Subset by FSW and client:
        if 'structuredsexual' in self.sim.networks.keys():
            fsw_infected = self.infected[self.sim.networks.structuredsexual.fsw]
            client_infected = self.infected[self.sim.networks.structuredsexual.client]

            # Add FSW and clients to results:
            if len(fsw_infected) > 0:
                self.results['prevalence_sw'][ti] = sum(fsw_infected) / len(fsw_infected)
                self.results['new_infections_sw'][ti] = len(((self.ti_infected == ti) & self.sim.networks.structuredsexual.fsw).uids)
                self.results['new_infections_not_sw'][ti] = len(((self.ti_infected == ti) & ~self.sim.networks.structuredsexual.fsw).uids)
            if len(client_infected) > 0:
                self.results['prevalence_client'][ti] = sum(client_infected) / len(client_infected)
                self.results['new_infections_client'][ti] = len(((self.ti_infected == ti) & self.sim.networks.structuredsexual.client).uids)
                self.results['new_infections_not_client'][ti] = len(((self.ti_infected == ti) & ~self.sim.networks.structuredsexual.client).uids)

        return

    def set_prognoses(self, uids, sources=None, ti=None):
        """
        Set prognoses upon infection
        """
        if ti is None:
            ti = self.ti
        else:
            # Check that ti is consistent with uids
            if not (sc.isnumber(ti) or len(ti) == len(uids)):
                errormsg = 'ti for set_prognoses must be int or array of length uids'
                raise ValueError(errormsg)

        self.susceptible[uids] = False
        self.infected[uids] = True
        self.acute[uids] = True

        self.ti_exposed[uids] = ti
        self.ti_infected[uids] = ti
        self.ti_acute[uids] = ti
        self.cd4[uids] = self.cd4_start[uids]

        # Set timing and CD4 count of latent infection
        dur_acute = self.pars.dur_acute.rvs(uids)
        self.ti_latent[uids] = self.ti_acute[uids] + dur_acute.astype(int)
        self.cd4_latent[uids] = self.pars.cd4_latent.rvs(uids)

        # Set time of onset of late-stage CD4 decline
        dur_latent = self.pars.dur_latent.rvs(uids)
        self.ti_falling[uids] = self.ti_latent[uids] + dur_latent.astype(int)
        dur_falling = self.pars.dur_falling.rvs(uids)
        self.ti_zero[uids] = self.ti_falling[uids] + dur_falling.astype(int)

        # Record source/target for ss.infection_log() -- HIV overrides the base
        # Disease.set_prognoses() entirely for its natural-history setup above, so
        # without this call the infection log's append hook never fires for HIV.
        super().set_prognoses(uids, sources)

        return

    def set_congenital(self, uids, sources):
        self.cd4_start[uids] = sc.dcp(self.cd4_start[sources])
        self.cd4_nadir[uids] = sc.dcp(self.cd4_start[uids])
        self.set_prognoses(uids, sources)
        return

    # Treatment-related changes
    def start_art(self, uids, p_effective_art=None):
        """
        Check who is ready to start ART treatment and put them on ART

        Args:
            uids: agents starting ART treatment
            p_effective_art (float/ss.Dist, optional): probability that each agent
                achieves viral suppression (effective ART) rather than non-suppressive
                ART. Pass a float to override `self.pars.p_effective_art`'s probability,
                or an already-initialized `ss.Dist` (e.g. an intervention's own par) to
                use it directly. Defaults to `self.pars.p_effective_art`
                (default: always effective).
        """
        ti = self.ti

        self.on_art[uids] = True
        self.art_discontinued[uids] = False
        newly_treated = uids[self.art_naive[uids]]
        self.art_naive[newly_treated] = False
        self.ti_art[uids] = ti
        self.cd4_preart[uids] = self.cd4[uids]

        # Determine whether each agent achieves viral suppression (effective ART) or
        # remains non-suppressive, e.g. due to poor adherence or drug resistance
        if isinstance(p_effective_art, ss.Dist):
            effective_uids = p_effective_art.filter(uids)
        else:
            if p_effective_art is not None:
                self.pars.p_effective_art.set(p=p_effective_art)
            effective_uids = self.pars.p_effective_art.filter(uids)
        nonsupp_uids = uids.remove(effective_uids)
        self.on_effective_art[effective_uids] = True
        self.on_nonsuppressive_art[effective_uids] = False
        self.on_effective_art[nonsupp_uids] = False
        self.on_nonsuppressive_art[nonsupp_uids] = True

        # Determine when agents goes off ART
        dur_on_art = self.pars.dur_on_art.rvs(uids)
        dur_on_art = dur_on_art * self.pars.rel_dur_on_art
        trend = self.pars.dur_on_art_trend
        if trend is not None:
            year = self.t.now('year')
            time_scale = np.interp(year, trend.years, trend.vals)
            dur_on_art = dur_on_art * time_scale
        self.ti_stop_art[uids] = ti + dur_on_art.astype(int)

        # ART nullifies all states and all future dates in the natural history
        self.acute[uids] = False
        self.latent[uids] = False
        self.falling[uids] = False
        future_latent = uids[self.ti_latent[uids] > ti]
        self.ti_latent[future_latent] = np.nan
        future_falling = uids[self.ti_falling[uids] > ti]
        self.ti_falling[future_falling] = np.nan
        future_zero = uids[self.ti_zero[uids] > ti]  # NB, if they are scheduled to die on this time step, they will
        self.ti_zero[future_zero] = np.nan

        # Set CD4 potential for anyone new to treatment - retreated people have the same potential
        # Extract growth parameters
        if len(newly_treated) > 0:
            cd4_max = self.pars.art_cd4_pars['cd4_max']
            cd4_healthy = self.pars.art_cd4_pars['cd4_healthy']
            cd4_preart = self.cd4_preart[newly_treated]

            # Calculate potential CD4 increase - assuming that growth follows the concave part of a logistic function
            # and that the total gain depends on the CD4 count at initiation
            if (cd4_preart == 0).any():
                raise ValueError('CD4 count is zero at initiation of ART')
            cd4_scale_factor = (cd4_max-cd4_preart)/cd4_healthy*np.log(cd4_max/cd4_preart)
            cd4_total_gain = cd4_preart*cd4_scale_factor
            self.cd4_potential[newly_treated] = self.cd4_preart[newly_treated] + cd4_total_gain

        return

    def stop_art(self, uids=None):
        """
        Check who is stopping ART treatment and put them off ART
        """
        ti = self.ti

        # Remove agents from ART
        if uids is None: uids = self.on_art & (self.ti_stop_art <= ti)
        self.on_art[uids] = False
        self.on_effective_art[uids] = False
        self.on_nonsuppressive_art[uids] = False
        self.art_discontinued[uids] = True
        self.ti_stop_art[uids] = ti
        self.cd4_postart[uids] = sc.dcp(self.cd4[uids])

        # Set decline
        dur_mean = np.log(self.cd4_preart[uids])*self.cd4[uids]/self.cd4_potential[uids]
        dur_scale = dur_mean * self.pars.dur_post_art_scale_factor
        dur_mean = ss.years(dur_mean)
        dur_scale = ss.years(np.maximum(dur_scale, 1e-3))  # Ensure it's not zero
        self.pars.dur_post_art.set(loc=dur_mean, scale=dur_scale)
        dur_post_art = self.pars.dur_post_art.rvs(uids)
        if np.isnan(dur_post_art).any():
            errormsg = 'Invalid entry for post-ART duration'
            raise ValueError(errormsg)
        self.ti_zero[uids] = ti + dur_post_art.astype(int)

        return


