"""
Starsim 3 base class for STIs, kept for the diseases not yet ported to Starsim 4 (HIV, syphilis, BV)
"""

import numpy as np
import starsim as ss
import sciris as sc
import stisim.utils as ut
import pandas as pd

ss_float_ = ss.dtypes.float

__all__ = []  # Internal: see stisim.diseases.sti for the Starsim 4 base classes

# Default stratifications for disease results. Pass age_bins=None or sex_keys=None to a
# disease module to skip creating the corresponding stratified results.
default_age_bins = [0, 15, 20, 25, 30, 35, 50, 65, 100]
default_sex_keys = {'': 'alive', 'f': 'female', 'm': 'male'}


class BaseSTIPars(ss.Pars):
    """
    Base parameters shared by all STI disease modules.

    Holds transmission parameters (beta values, condom efficacy), timestep
    settings, and initial prevalence configuration.
    """
    def __init__(self, **kwargs):
        super().__init__()

        # Settings
        self.include_care = True  # Determines whether testing results are included

        # Time
        self.dt = 'month'

        # Transmission
        self.beta = 0  # Placeholder: no transmission. This will be set in validate_beta
        self.beta_m2f = None
        self.rel_beta_f2m = 0.5
        self.beta_m2c = None
        self.beta_breastfeed = None
        self.beta_m2m = None
        self.eff_condom = 1
        self.rel_init_prev = 1
        self.update(kwargs)
        return


# Main class
class BaseSTI(ss.Infection):
    """
    Base class for sexually transmitted infections.

    Extends ``ss.Infection`` with STI-specific transmission logic, including
    network-specific beta values, condom efficacy, and initial prevalence by
    sex and risk group. All STIsim disease modules inherit from this class.

    Args:
        name (str): Module name used for results and parameter routing.
        pars (dict): Override default parameters from ``BaseSTIPars``.
        init_prev_data: Initial prevalence data (float or DataFrame by age/sex/risk group).
        age_range (list): Age range used for the main (adult) prevalence results. Default: [15, 50].
        age_bins (list): Bin edges for age-stratified results; ``None`` skips them entirely.
        sex_keys (dict): Sex strata for results, mapping result suffix to a people attribute;
            ``None`` skips sex-stratified results, keeping only the whole-population ones.
        **kwargs: Additional parameters passed to ``update_pars``.
    """
    def __init__(self, name=None, pars=None, init_prev_data=None, age_range=None,
                 age_bins=default_age_bins, sex_keys=default_sex_keys, **kwargs):
        super().__init__(name=name)

        # Handle parameters
        default_pars = BaseSTIPars()
        self.define_pars(**default_pars)
        self.update_pars(pars, **kwargs)

        self.define_states(
            # Time of acquisition: S->E, or S->I for diseases with no exposed compartment.
            # Written exactly once per infection and never rewritten, which is what makes
            # `ti_exposed == ti` a correct one-shot test for incidence. In `SEIS` the time of
            # becoming infectious is `ti_infectious`, which `step_state()` snaps forward to the
            # realized transition step, so it cannot be used for that.
            ss.FloatArr('ti_exposed'),
            ss.FloatArr('ti_transmitted_sex'),
            ss.FloatArr('ti_transmitted_mtc'),
            ss.FloatArr('ti_transmitted'),
            ss.FloatArr('new_transmissions'),
            ss.FloatArr('new_transmissions_sex'),
            ss.FloatArr('infection_net'),  # Network index at time of infection (NaN if uninfected)
        )

        # Set initial prevalence
        self.init_prev_data = init_prev_data

        # Results. age_bins=None skips age-stratified results; sex_keys=None skips
        # sex-stratified results, leaving only the whole-population ones.
        self.age_range = age_range if age_range is not None else [15, 50]  # Age range for main results e.g. prevalence
        self.age_bins = np.array(age_bins) if age_bins is not None else None
        self.sex_keys = dict(sex_keys) if sex_keys is not None else {'': 'alive'}

        return

    def make_init_prev(self, uids=None, data=None, active=True):
        """ Initialize prevalence by sex and risk group """
        sim = self.sim
        if data is None: data = self.init_prev_data
        if uids is None: uids = sim.people.auids  # Everyone

        if sc.isnumber(data):
            init_prev = data

        elif isinstance(data, pd.DataFrame):

            init_prev = np.zeros(len(uids), dtype=ss_float_)
            df = data

            nw = sim.networks.structuredsexual
            n_risk_groups = nw.pars.n_risk_groups
            for rg in range(n_risk_groups):
                for sex in ['female', 'male']:
                    for sw in [0, 1]:
                        thisdf = df.loc[(df.risk_group==rg) & (df.sex==sex) & (df.sw==sw)]
                        conditions = sim.people[sex] & (nw.risk_group==rg)
                        if active:
                            conditions = conditions & nw.active(sim.people)
                        if sw:
                            if sex == 'female': conditions = conditions & sim.networks.structuredsexual.fsw
                            if sex == 'male':   conditions = conditions & sim.networks.structuredsexual.client
                        init_prev[conditions[uids]] = thisdf.init_prev.values[0]

        else:
            errormsg = 'Format of init_prev_data must be float or dataframe.'
            raise ValueError(errormsg)

        # Scale and validate
        init_prev = init_prev * self.pars.rel_init_prev
        init_prev = np.clip(init_prev, a_min=0, a_max=1)

        return init_prev

    def validate_beta(self):
        """
        Map disease-specific transmission parameters onto the betamap.

        The parent class (ss.Infection) builds a betamap keyed by standardized
        network name, with [beta_p1→p2, beta_p2→p1] per network. Here we
        override entries for STI-specific networks using disease parameters:

            - structuredsexual: beta_m2f (and rel_beta_f2m for female→male)
            - maternal: beta_m2c (prenatal MTCT)
            - breastfeeding: beta_breastfeed (postnatal MTCT, HIV only)
            - msm: beta_m2m (bidirectional)
        """
        betamap = super().validate_beta()
        if self.pars.beta_m2f is not None:
            for name,net in self.sim.networks.items(): # Every network between men and women (sides m and f; men are p1)
                key = ss.standardize_netkey(name)
                if key == 'structuredsexual' or getattr(net, 'sides', None) == dict(m='male', f='female'):
                    betamap[key][0] = self.pars.beta_m2f
                    betamap[key][1] = self.pars.beta_m2f * self.pars.rel_beta_f2m
        if self.pars.beta_m2c is not None and 'maternal' in betamap:
            betamap['maternal'][0] = ss.permonth(self.pars.beta_m2c)
        if self.pars.beta_breastfeed is not None and 'breastfeeding' in betamap:
            betamap['breastfeeding'][0] = ss.permonth(self.pars.beta_breastfeed)
        if self.pars.beta_m2m is not None and 'msm' in betamap:
            betamap['msm'][0] = self.pars.beta_m2m
            betamap['msm'][1] = self.pars.beta_m2m
        return betamap

    def agehist(self, a):
        """ Return an age histogram """
        aa = self.sim.people.age[a]
        return np.histogram(aa, bins=self.age_bins)[0]

    @property
    def treatable(self):
        """ Assume infected people are treatable, can be overwritten in subclasses """
        return self.infected

    def init_results(self):
        """ Initialize results """
        super().init_results()
        results = sc.autolist()

        # Most results are stored by age and sex
        for sk in self.sex_keys.keys():
            skk = '' if sk == '' else f'_{sk}'
            skl = '' if sk == '' else f' ({sk.upper()})'
            is_subpop = skk != ''  # Sex-stratified results are not auto-plotted
            results += [
                ss.Result(f'incidence{skk}', scale=False, label=f"Incidence{skl}", auto_plot=not is_subpop),
            ]
            if is_subpop:
                results += [
                    ss.Result(f'new_infections{skk}', dtype=int, label=f"New infections{skl}", auto_plot=False),
                    ss.Result(f'prevalence{skk}', scale=False, label=f"Prevalence{skl}", auto_plot=False),
                    ss.Result(f'n_infected{skk}', dtype=int, label=f"Number infected{skl}", auto_plot=False),
                ]

            if self.age_bins is not None:
                for ab1,ab2 in zip(self.age_bins[:-1], self.age_bins[1:]):
                    ask = f'{skk}_{ab1}_{ab2}'
                    asl = f' ({skl}, {ab2}-{ab2})'
                    results += [
                        ss.Result(f'new_infections{ask}', dtype=int, label=f"New infections{asl}", auto_plot=False),
                        ss.Result(f'n_infected{ask}', dtype=int, label=f"Number infected{asl}", auto_plot=False),
                        ss.Result(f'incidence{ask}', scale=False, label=f"Incidence{asl}", auto_plot=False),
                        ss.Result(f'prevalence{ask}', scale=False, label=f"Prevalence{asl}", auto_plot=False),
                    ]

        if self.pars.include_care:
            for sk in self.sex_keys.keys():
                skk = '' if sk == '' else f'_{sk}'
                skl = '' if sk == '' else f' ({sk.upper()})'
                results += [
                    ss.Result('new_care_seekers'+skk, dtype=int, label="New care seekers"+skl, auto_plot=False),
                    ss.Result('new_false_pos'+skk, dtype=int, label="New false positives"+skl, auto_plot=False),
                    ss.Result('new_true_pos'+skk, dtype=int, label="New true positives"+skl, auto_plot=False),
                    ss.Result('new_false_neg'+skk, dtype=int, label="New false negatives"+skl, auto_plot=False),
                    ss.Result('new_true_neg'+skk, dtype=int, label="New true negatives"+skl, auto_plot=False),
                    ss.Result('new_treated_success'+skk, dtype=int, label="Successful treatments"+skl, auto_plot=False),
                    ss.Result('new_treated_failure'+skk, dtype=int, label="Unsuccessful treatments"+skl, auto_plot=False),
                    ss.Result('new_treated_unnecessary'+skk, dtype=int, label="Unnecessary treatments"+skl, auto_plot=False),
                    ss.Result('new_treated'+skk, dtype=int, label="Treatments"+skl, auto_plot=False),
                ]

        # Transmission route results
        results += [
            ss.Result('new_infections_sex', dtype=int, label='New infections (sexual)', auto_plot=False),
            ss.Result('new_infections_mtct', dtype=int, label='New infections (MTCT)', auto_plot=False),
            ss.Result('new_infections_prenatal', dtype=int, label='New infections (prenatal MTCT)', auto_plot=False),
            ss.Result('new_infections_postnatal', dtype=int, label='New infections (postnatal MTCT)', auto_plot=False),
        ]

        self.define_results(*results)

        return

    def step(self):
        """ Override to pass network indices to set_outcomes for MTCT tracking """
        new_cases, sources, networks = self.infect()
        if len(new_cases):
            self.set_outcomes(new_cases, sources, net_ids=networks)
        return new_cases, sources, networks

    def set_outcomes(self, uids, sources, net_ids=None):
        """
        Set prognoses for newly infected agents and track transmission routes.

        Args:
            uids: UIDs of agents newly infected this timestep (not all agents)
            sources: UIDs of the agents who transmitted to each new case
            net_ids: array of network indices indicating which network each
                     infection came from (used to distinguish prenatal vs
                     postnatal MTCT). Passed from step() via infect().
        """
        super().set_outcomes(uids, sources=sources)
        self.new_transmissions[:] = 0
        self.new_transmissions_sex[:] = 0

        # Classify each new infection by its transmission route.
        # net_ids is an integer array (one per new infection) where each value
        # is the index into sim.networks (in insertion order) of the network
        # that caused that infection. We loop over the unique network indices,
        # look up the actual network object, and use isinstance to classify
        # as prenatal MTCT, postnatal MTCT, or horizontal (sexual) transmission.
        is_prenatal = np.zeros(len(uids), dtype=bool)
        is_postnatal = np.zeros(len(uids), dtype=bool)
        if net_ids is not None:
            net_list = list(self.sim.networks.values())
            for idx in np.unique(net_ids):
                net = net_list[idx]
                mask = net_ids == idx
                if isinstance(net, ss.PrenatalNet):
                    is_prenatal |= mask
                elif isinstance(net, ss.PostnatalNet):
                    is_postnatal |= mask

        is_mtct = is_prenatal | is_postnatal

        # Store network index per agent for use in update_results
        if net_ids is not None:
            self.infection_net[uids] = net_ids.astype(float)

        sex_transmission = sources[~is_mtct]
        mtc_transmission = sources[is_mtct]

        # Track per-source transmission counts
        unique_sources, counts = np.unique(sources, return_counts=True)
        self.ti_transmitted[unique_sources] = self.ti
        unique_sources_sex, counts_sex = np.unique(sex_transmission, return_counts=True)
        self.ti_transmitted_sex[sex_transmission] = self.ti
        self.ti_transmitted_mtc[mtc_transmission] = self.ti
        self.new_transmissions[unique_sources] = counts
        self.new_transmissions_sex[unique_sources_sex] = counts_sex

        return

    def update_results(self):
        super().update_results()
        ti = self.ti
        ppl = self.sim.people

        adults = (self.sim.people.age >= self.age_range[0]) & (self.sim.people.age <= self.age_range[1])

        # Incidence is counted at acquisition. The skk == '' pass deliberately overwrites the
        # `new_infections` result that `ss.Infection.set_prognoses()` accumulates (starsim >=
        # 3.6.0), which excludes seed infections as prevalent rather than incident; STIsim counts
        # them, so that `new_infections` agrees with the stratified `new_infections_*` results.
        newly_infected = self.ti_exposed == (ti - 1)

        # Main results, looping over people keys and attributes.
        for pkey, pattr in self.sex_keys.items():
            skk = '' if pkey == '' else f'_{pkey}'

            # Collate results
            new_inf = newly_infected & ppl[pattr]
            n_sus = self.susceptible & ppl[pattr]
            n_inf = self.infected & ppl[pattr]

            # Store main results
            self.results[f'prevalence{skk}'][ti] = ut.cond_prob(self.infected, adults & ppl[pattr])
            self.results[f'new_infections{skk}'][ti] = ut.count(new_inf)
            self.results[f'n_infected{skk}'][ti] = ut.count(n_inf)
            self.results[f'incidence{skk}'][ti] = ut.countdiv(new_inf, n_sus)

            # Compute and store age results
            if self.age_bins is not None:
                age_results = dict(
                    new_infections  = self.agehist(new_inf),
                    n_infected      = self.agehist(n_inf),
                    incidence       = ut.div(self.agehist(new_inf), self.agehist(n_sus)),
                    prevalence      = ut.div(self.agehist(n_inf), self.agehist(ppl[pattr])),
                )

                for akey, ares in age_results.items():
                    ai = 0
                    for ab1, ab2 in zip(self.age_bins[:-1], self.age_bins[1:]):
                        ask = f'{skk}_{ab1}_{ab2}'
                        self.results[f'{akey}{ask}'][ti] = ares[ai]
                        ai += 1

        # Transmission route results: classify new infections by the network
        # they came from (stored in infection_net by set_outcomes).
        # TODO: refactor once starsim tracks transmission route natively
        #   (see https://github.com/starsimhub/starsim/issues/1232)
        new_this_ti = newly_infected & ppl.alive # As new_infections counts them
        n_prenatal = 0
        n_postnatal = 0
        n_horizontal = 0
        if ut.count(new_this_ti):
            for idx, net in enumerate(self.sim.networks.values()):
                infected_via = new_this_ti & (self.infection_net == float(idx))
                n = ut.count(infected_via)
                if isinstance(net, ss.PrenatalNet):
                    n_prenatal += n
                elif isinstance(net, ss.PostnatalNet):
                    n_postnatal += n
                else:
                    n_horizontal += n
        self.results['new_infections_prenatal'][ti] = n_prenatal
        self.results['new_infections_postnatal'][ti] = n_postnatal
        self.results['new_infections_mtct'][ti] = n_prenatal + n_postnatal
        self.results['new_infections_sex'][ti] = n_horizontal

        return
