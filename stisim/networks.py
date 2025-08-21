"""
Define sexual network for STI transmission.

Overview:

- Risk groups: agents are randomly assigned into one of 3 main risk groups:

    - 0 = marry and remain married to a single partner throughout their lifetime
    - 1 = marry and then divorce or who have concurrent partner(s) during their marriage
    - 2 = never marry

- In addition, a proportion of each of the groups above engages in sex work.
"""

import starsim as ss
import sciris as sc
import numpy as np
import pandas as pd
from collections import defaultdict
from bisect import bisect_left

ss_float_ = ss.dtypes.float
ss_int_ = ss.dtypes.int

# Specify all externally visible functions this file defines; see also more definitions below
__all__ = ['NetworkPars', 'StructuredSexual', 'PriorPartners', 'AgeMatchedMSM', 'AgeApproxMSM']


class NoPartnersFound(Exception):
    # Raise this exception if the matching algorithm wasn't able to match any partners
    pass


class NetworkPars(ss.Pars):
    def __init__(self, **kwargs):
        super().__init__()

        # Settings
        self.recall_prior = False  # Whether to recall prior partners
        self.n_risk_groups = 3  # Number of risk groups

        self.f_age_group_bins = dict(  # For separating women into age groups: teens, young women, adult women
            teens=(0, 20),
            young=(20, 25),
            adult=(25, np.inf),
        )

        # Age of sexual debut
        self.debut = ss.lognorm_ex(20, 3)
        self.debut_pars_f = [20, 3]
        self.debut_pars_m = [21, 3]

        # Risk groups
        self.p_lo_risk = ss.bernoulli(p=0)
        self.p_hi_risk = ss.bernoulli(p=0)
        self.prop_f0 = 0.85
        self.prop_m0 = 0.8
        self.prop_f2 = 0.01
        self.prop_m2 = 0.02

        # Age difference preferences
        self.age_diff_pars = dict(
            teens=[(7, 3), (6, 3), (5, 1)],  # (mu,stdev) for levels 0, 1, 2
            young=[(8, 3), (7, 3), (5, 2)],
            adult=[(8, 3), (7, 3), (5, 2)],
        )

        # Concurrency preferences
        self.concurrency_dist = ss.poisson(lam=1)
        self.f0_conc = 0.0001
        self.f1_conc = 0.01
        self.f2_conc = 0.1
        self.m0_conc = 0.0001
        self.m1_conc = 0.2
        self.m2_conc = 0.5

        # Relationship initiation, stability, and duration
        self.p_pair_form = ss.bernoulli(p=0.5)  # Probability of a (stable) pair forming between two matched people
        self.match_dist = ss.bernoulli(p=0)  # Placeholder value replaced by risk-group stratified values below
        self.p_matched_stable = [0.9, 0.5, 0]  # Probability of a stable pair forming between matched people (otherwise casual)
        self.p_mismatched_casual = [0.5, 0.5, 0.5]  # Probability of a casual pair forming between mismatched people (otherwise instantanous)

        # Durations of stable and casual relationships
        self.stable_dur_pars = dict(
            teens=[
                # (mu,stdev) for levels 0, 1, 2
                [ss.years(100),  ss.years(1)],
                [ss.years(8),  ss.years(2)],
                [ss.months(1e-4), ss.months(1e-4)]
            ],
            young=[
                [ss.years(100),  ss.years(1)],
                [ss.years(10),  ss.years(3)],
                [ss.months(1e-4), ss.months(1e-4)]
            ],
            adult=[
                [ss.years(100),  ss.years(1)],
                [ss.years(12),  ss.years(3)],
                [ss.months(1e-4), ss.months(1e-4)]
            ],
        )
        self.casual_dur_pars = dict(
            teens=[[ss.years(1), ss.years(3)]]*3,
            young=[[ss.years(1), ss.years(3)]]*3,
            adult=[[ss.years(1), ss.years(3)]]*3,
        )

        # Acts
        self.acts = ss.lognorm_ex(ss.freqperyear(80), ss.freqperyear(30))  # Coital acts/year

        # Condoms
        self.condom_data = None

        # Sex work parameters
        self.fsw_shares = ss.bernoulli(p=0.05)
        self.client_shares = ss.bernoulli(p=0.12)
        self.sw_seeking_rate = ss.probpermonth(1.0)  # Monthly rate at which clients seek FSWs (1 new SW partner / month)
        self.sw_seeking_dist = ss.bernoulli(p=0.5)  # Placeholder value replaced by dt-adjusted sw_seeking_rate
        self.sw_beta = 1
        self.sw_intensity = ss.random()  # At each time step, FSW may work with varying intensity

        # Distributions derived from parameters above - don't adjust
        self.age_diffs = ss.normal()
        self.dur_dist = ss.lognorm_ex()
        self.update(kwargs)
        return


class StructuredSexual(ss.SexualNetwork):
    """
    Structured sexual network
    """

    def __init__(self, pars=None, condom_data=None, name=None, **kwargs):

        super().__init__(name=name)

        # Set edge attributes
        self.meta.sw = bool
        self.meta.condoms = ss_float_
        self.meta.age_p1 = ss_float_
        self.meta.age_p2 = ss_float_
        self.meta.edge_type = ss_float_  # edge type tracks stable/casual/onetime

        # Set parameters
        default_pars = NetworkPars()
        self.define_pars(**default_pars)
        self.update_pars(pars, **kwargs)

        # Set condom use
        if self.pars.condom_data is not None:
            self.pars.condom_data = self.process_condom_data(self.pars.condom_data)

        self.edge_types = {'stable': 0, 'casual': 1, 'onetime': 2, 'sw': 3}

        # Add states
        self.define_states(
            ss.BoolArr('participant', default=True),
            ss.FloatArr('debut', default=0),
            ss.FloatArr('risk_group'),  # Which risk group an agent belongs to
            ss.BoolArr('fsw'),  # Whether an agent is a female sex worker
            ss.BoolArr('client'),  # Whether an agent is a client of sex workers
            ss.FloatArr('concurrency'),  # Preferred number of concurrent partners
            ss.FloatArr('partners', default=0),  # Actual number of concurrent partners
            ss.FloatArr('partners_12', default=0),  # Number of partners over the past 12m
            ss.FloatArr('lifetime_partners', default=0),  # Lifetime total number of partners
            ss.FloatArr('casual_partners', default=0),
            ss.FloatArr('stable_partners', default=0),
            ss.FloatArr('onetime_partners', default=0),
            ss.FloatArr('sw_partners', default=0),
            ss.FloatArr('lifetime_casual_partners', default=0),
            ss.FloatArr('lifetime_stable_partners', default=0),
            ss.FloatArr('lifetime_onetime_partners', default=0),
            ss.FloatArr('lifetime_sw_partners', default=0),
            ss.FloatArr('sw_intensity'),  # Intensity of sex work
            reset = True, # To allow redefining participant
        )

        self.relationship_durs = defaultdict(list)

        return

    @staticmethod
    def process_condom_data(condom_data):
        if sc.isnumber(condom_data):
            return condom_data
        elif isinstance(condom_data, pd.DataFrame):
            df = condom_data.melt(id_vars=['partnership'])
            dd = dict()
            for pcombo in df.partnership.unique():
                key = tuple(map(int, pcombo[1:-1].split(','))) if pcombo != '(fsw,client)' else ('fsw','client')
                thisdf = df.loc[df.partnership == pcombo]
                dd[key] = dict()
                dd[key]['year'] = thisdf.variable.values.astype(int)
                dd[key]['val'] = thisdf.value.values
        return dd

    def get_age_risk_pars(self, uids, par):
        loc = np.full(uids.shape, fill_value=np.nan)
        scale = np.full(uids.shape, fill_value=np.nan)
        for a_label, (age_lower, age_upper) in self.pars.f_age_group_bins.items():
            for rg in range(self.pars.n_risk_groups):
                in_risk_group = (self.sim.people.age[uids] >= age_lower) & (self.sim.people.age[uids] < age_upper) & (self.risk_group[uids] == rg)
                p0 = par[a_label][rg][0]
                p1 = par[a_label][rg][1]
                # Scale the parameters by the time step if specified
                # TODO: fix this
                if isinstance(p0, ss.dur):
                    p0 = p0.months
                    p1 = p1.months
                loc[in_risk_group] = p0
                scale[in_risk_group] = p1
        if np.isnan(scale).any() or np.isnan(loc).any():
            errormsg = 'Invalid entries for age difference preferences.'
            raise ValueError(errormsg)
        return loc, scale

    def init_pre(self, sim):
        super().init_pre(sim)

        # Checks
        if self.pars.recall_prior:
            isprior = [isinstance(nw, PriorPartners) for nw in self.sim.networks.values()]
            if not any(isprior):
                errormsg = 'PriorPartners network is required if recall_prior is True.'
                raise ValueError(errormsg)

        # Process condom data
        if self.pars.condom_data is not None:
            if isinstance(self.pars.condom_data, dict):
                for rgtuple, valdict in self.pars.condom_data.items():
                    yearvec = self.t.yearvec
                    self.pars.condom_data[rgtuple]['simvals'] = sc.smoothinterp(yearvec, valdict['year'], valdict['val'])
        return

    def init_post(self):
        super().init_post(add_pairs=False)
        self.set_network_states()
        return

    def set_network_states(self, upper_age=None):
        self.set_risk_groups(upper_age=upper_age)
        self.set_concurrency(upper_age=upper_age)
        self.set_sex_work(upper_age=upper_age)
        self.set_debut(upper_age=upper_age)
        return

    @property
    def over_debut(self):
        return self.sim.people.age > self.debut

    def _get_uids(self, upper_age=None, by_sex=True):
        people = self.sim.people
        if upper_age is None: upper_age = 1000
        within_age = people.age <= upper_age
        if by_sex:
            f_uids = (within_age & people.female).uids
            m_uids = (within_age & people.male).uids
            return f_uids, m_uids
        else:
            uids = within_age.uids
            return uids

    def set_risk_groups(self, upper_age=None):
        """ Assign each person to a risk group """
        ppl = self.sim.people
        uids = self._get_uids(upper_age=upper_age, by_sex=False)

        p_lo = np.full(len(uids), fill_value=np.nan, dtype=ss_float_)
        p_lo[ppl.female[uids]] = self.pars.prop_f0
        p_lo[ppl.male[uids]] = self.pars.prop_m0
        self.pars.p_lo_risk.set(p=p_lo)
        lo_risk, hi_med_risk = self.pars.p_lo_risk.split(uids)

        p_hi = np.full(len(hi_med_risk), fill_value=np.nan, dtype=ss_float_)
        p_hi[ppl.female[hi_med_risk]] = self.pars.prop_f2/(1-self.pars.prop_f0)
        p_hi[ppl.male[hi_med_risk]] = self.pars.prop_m2/(1-self.pars.prop_m0)
        self.pars.p_hi_risk.set(p=p_hi)
        hi_risk, med_risk = self.pars.p_hi_risk.split(hi_med_risk)

        self.risk_group[lo_risk] = 0
        self.risk_group[med_risk] = 1
        self.risk_group[hi_risk] = 2
        return

    def set_concurrency(self, upper_age=None):
        """ Assign each person a preferred number of simultaneous partners """
        people = self.sim.people
        if upper_age is None: upper_age = 1000
        in_age_lim = (people.age < upper_age)
        uids = in_age_lim.uids

        lam = np.full(uids.shape, fill_value=np.nan, dtype=ss_float_)
        for rg in range(self.pars.n_risk_groups):
            f_conc = self.pars[f'f{rg}_conc']
            m_conc = self.pars[f'm{rg}_conc']
            in_risk_group = self.risk_group == rg
            in_group = in_risk_group & in_age_lim
            f_in = (people.female & in_group)[uids]
            m_in = (people.male   & in_group)[uids]
            if f_in.any(): lam[f_in] = f_conc
            if m_in.any(): lam[m_in] = m_conc

        self.pars.concurrency_dist.set(lam=lam)
        self.concurrency[uids] = self.pars.concurrency_dist.rvs(uids) + 1

        return

    def set_sex_work(self, upper_age=None):
        f_uids, m_uids = self._get_uids(upper_age=upper_age)
        self.fsw[f_uids] = self.pars.fsw_shares.rvs(f_uids)
        self.client[m_uids] = self.pars.client_shares.rvs(m_uids)
        return

    def set_debut(self, upper_age=None):
        uids = self._get_uids(upper_age=upper_age, by_sex=False)
        par1 = np.full(len(uids), fill_value=np.nan, dtype=ss_float_)
        par2 = np.full(len(uids), fill_value=np.nan, dtype=ss_float_)
        par1[self.sim.people.female[uids]] = self.pars.debut_pars_f[0]
        par2[self.sim.people.female[uids]] = self.pars.debut_pars_f[1]
        par1[self.sim.people.male[uids]] = self.pars.debut_pars_m[0]
        par2[self.sim.people.male[uids]] = self.pars.debut_pars_m[1]
        self.pars.debut.set(mean=par1, std=par2)
        self.debut[uids] = self.pars.debut.rvs(uids)
        return

    def match_pairs(self):
        """
        Match pairs by age, using sorting rather than the linear sum assignment
        """
        ppl = self.sim.people

        # Find people eligible for a relationship
        active = self.over_debut
        underpartnered = self.partners < self.concurrency
        f_eligible = active & ppl.female & underpartnered
        m_eligible = active & ppl.male & underpartnered
        f_looking = self.pars.p_pair_form.filter(f_eligible.uids)  # ss.uids of women looking for partners

        if len(f_looking) == 0 or m_eligible.count() == 0:
            raise NoPartnersFound()

        # Get mean age differences and desired ages
        loc, scale = self.get_age_risk_pars(f_looking, self.pars.age_diff_pars)
        self.pars.age_diffs.set(loc=loc, scale=scale)
        age_gaps = self.pars.age_diffs.rvs(f_looking)   # Sample the age differences
        desired_ages = ppl.age[f_looking] + age_gaps    # Desired ages of the male partners
        m_ages = ppl.age[m_eligible]            # Ages of eligible males
        ind_m = np.argsort(m_ages, stable=True)
        ind_f = np.argsort(desired_ages, stable=True)

        # If there are no agents in either group, return empty arrays
        if len(ind_m) == 0 or len(ind_f) == 0:
            raise NoPartnersFound()

        # drop all males that are younger than one standard deviation below the lowest desired age.
        youngest_preferred_male_age = desired_ages[ind_f[0]]
        youngest_male_age = m_ages[ind_m[0]]

        if youngest_male_age < youngest_preferred_male_age:
            # remove the youngest males until the youngest is at least as old as the youngest preferred
            cutoff_index = bisect_left(m_ages[ind_m], youngest_preferred_male_age)
            ind_m = ind_m[cutoff_index:]

        elif youngest_preferred_male_age < youngest_male_age:
            # remove the youngest females until the youngest preferred is at least as old as the youngest male
            cutoff_index = bisect_left(desired_ages[ind_f], youngest_male_age)
            ind_f = ind_f[cutoff_index:]

        # Check again for empty arrays after filtering
        if len(ind_m) == 0 or len(ind_f) == 0:
            raise NoPartnersFound()

        # Check the upper limit of the age spectrum
        oldest_preferred_male_age = desired_ages[ind_f[-1]]
        oldest_male_age = m_ages[ind_m[-1]]

        if oldest_male_age > oldest_preferred_male_age:
            cutoff_index = bisect_left(m_ages[ind_m], oldest_preferred_male_age)
            ind_m = ind_m[:cutoff_index]

        elif oldest_preferred_male_age > oldest_male_age:
            cutoff_index = bisect_left(desired_ages[ind_f], oldest_male_age)
            ind_f = ind_f[:cutoff_index]

        # draw n samples from the larger of the two groups, where n is the number of samples in the smaller group
        if len(ind_m) < len(ind_f):
            ind_f_subset = np.random.choice(len(ind_f), size=len(ind_m), replace=False)
            ind_f_subset.sort()
            ind_f = ind_f[ind_f_subset]
        elif len(ind_f) < len(ind_m):
            ind_m_subset = np.random.choice(len(ind_m), size=len(ind_f), replace=False)
            ind_m_subset.sort()
            ind_m = ind_m[ind_m_subset]

        if len(ind_m) == 0 or len(ind_f) == 0:
            raise NoPartnersFound()

        p1 = m_eligible.uids[ind_m]
        p2 = f_looking[ind_f]

        return p1, p2

    def add_pairs_sw(self):
        ppl = self.sim.people

        try:
            p1, p2 = self.match_sex_workers()
        except NoPartnersFound:
            return

        match_count = len(p1)
        beta = np.ones(match_count, dtype=ss_float_)
        condoms = np.zeros(match_count, dtype=ss_float_)
        acts = (self.pars.acts.rvs(p2)).astype(int)
        dur = np.full(match_count, fill_value=1)  # Measured in timesteps
        age_p1 = ppl.age[p1]
        age_p2 = ppl.age[p2]
        edge_types = np.full(match_count, dtype=ss_float_, fill_value=self.edge_types['sw'])

        self.append(p1=p1, p2=p2, beta=beta, condoms=condoms, dur=dur, acts=acts, sw=[True]*match_count, age_p1=age_p1, age_p2=age_p2, edge_type=edge_types)

        p1_edges, p1_counts = np.unique(p1, return_counts=True)
        p2_edges, p2_counts = np.unique(p2, return_counts=True)

        # update partner counts
        self.lifetime_sw_partners[p1_edges] += p1_counts
        self.lifetime_sw_partners[p2_edges] += p2_counts

        return

    def add_pairs_nonsw(self):
        ppl = self.sim.people

        try:
            p1, p2 = self.match_pairs()
        except NoPartnersFound:
            return

        matched_risk = (self.risk_group[p1] == self.risk_group[p2])
        mismatched_risk = (self.risk_group[p1] != self.risk_group[p2])

        # Set the probability of forming a partnership
        p_match = np.full(len(p1), fill_value=np.nan, dtype=ss_float_)
        for rg in range(self.pars.n_risk_groups):
            p_match[matched_risk & (self.risk_group[p1] == rg)] = self.pars.p_matched_stable[rg]
            p_match[mismatched_risk & (self.risk_group[p2] == rg)] = self.pars.p_mismatched_casual[rg]
        self.pars.match_dist.set(p=p_match)
        matches = self.pars.match_dist.rvs(p2)

        stable = matches & matched_risk
        casual = matches & mismatched_risk
        any_match = stable | casual

        match_count = len(p1)

        beta = np.ones(match_count, dtype=ss_float_)
        condoms = np.zeros(match_count, dtype=ss_float_)
        acts = (self.pars.acts.rvs(p2)).astype(int)
        dur = np.full(match_count, fill_value=1)  # Measured in timesteps
        age_p1 = ppl.age[p1]
        age_p2 = ppl.age[p2]
        edge_types = np.full(match_count, dtype=ss_float_, fill_value=np.nan)
        edge_types[stable] = self.edge_types['stable']
        edge_types[casual] = self.edge_types['casual']

        # Set duration
        dur_mean = np.full(match_count, fill_value=np.nan, dtype=ss_float_)
        dur_std = np.full(match_count, fill_value=np.nan, dtype=ss_float_)
        for which, bools in {'stable': stable, 'casual': casual}.items():
            if bools.any():
                uids = p2[bools]
                thesepars = self.pars[f'{which}_dur_pars']
                mean, std = self.get_age_risk_pars(uids, thesepars)
                dur_mean[bools] = mean
                dur_std[bools] = std
        self.pars.dur_dist.set(mean=dur_mean[any_match], std=dur_std[any_match])
        dur[any_match] = self.pars.dur_dist.rvs(p2[any_match])
        # dur[any_match] = self.pars.dur_dist.rvs(sum(any_match))

        edge_types[(dur == 1)] = self.edge_types['onetime']

        # track the duration of all new relationships
        relationships = (dur > 1)

        for (a, b, reldur) in zip(p1[relationships], p2[relationships], dur[relationships]):
            pair = (min(a,b), max(a,b))
            self.relationship_durs[pair].append({'start': self.ti, 'dur': reldur}) # set dur to intended duration. When the relationship actually ends, this will be updated

        self.append(p1=p1, p2=p2, beta=beta, condoms=condoms, dur=dur, acts=acts, sw=[False]*match_count, age_p1=age_p1, age_p2=age_p2, edge_type=edge_types)

        # Checks
        if (self.sim.people.female[p1].any() or self.sim.people.male[p2].any()) and (self.name == 'structuredsexual'):
            errormsg = 'Same-sex pairings should not be possible in this network'
            raise ValueError(errormsg)
        if len(p1) != len(p2):
            errormsg = 'Unequal lengths in edge list'
            raise ValueError(errormsg)

        # update partner counts
        for key, edge_type in self.edge_types.items():
            p1_edges = p1[edge_types==edge_type]
            p2_edges = p2[edge_types==edge_type]
            self.partners[p1_edges] += 1
            self.partners[p2_edges] += 1
            self.lifetime_partners[p1_edges] += 1
            self.lifetime_partners[p2_edges] += 1
            getattr(self, f'{key}_partners')[p1_edges] += 1
            getattr(self, f'{key}_partners')[p2_edges] += 1
            getattr(self, f'lifetime_{key}_partners')[p1_edges] += 1
            getattr(self, f'lifetime_{key}_partners')[p2_edges] += 1

        return

    def add_pairs(self):
        self.add_pairs_nonsw()
        self.add_pairs_sw()
        return

    def match_sex_workers(self):
        """ Match sex workers to clients """

        # Find people eligible for a relationship
        active = self.over_debut
        active_fsw = active & self.fsw
        active_clients = active & self.client
        self.sw_intensity[active_fsw.uids] = self.pars.sw_intensity.rvs(active_fsw.uids)

        # Find clients who will seek FSW
        self.pars.sw_seeking_dist.pars.p = self.pars.sw_seeking_rate.to_prob()
        m_looking = self.pars.sw_seeking_dist.filter(active_clients.uids)

        if len(m_looking) == 0 or len(active_fsw.uids) == 0:
            raise NoPartnersFound()

        # Attempt to assign a sex worker to every client by repeat sampling the sex workers.
        # FSW with higher work intensity will be sampled more frequently
        if len(m_looking) > len(active_fsw.uids):
            n_repeats = (self.sw_intensity[active_fsw]*10).astype(int)+1
            fsw_repeats = np.repeat(active_fsw.uids, n_repeats)
            if len(fsw_repeats) < len(m_looking):
                fsw_repeats = np.repeat(fsw_repeats, 10)  # 10x the number of clients each sex worker can have

            # Might still not have enough FSW, so form as many pairs as possible
            n_pairs = min(len(fsw_repeats), len(m_looking))
            if len(fsw_repeats) < len(m_looking):
                p1 = m_looking[:n_pairs]
                p2 = fsw_repeats
            else:
                unique_sw, counts_sw = np.unique(fsw_repeats, return_counts=True)
                count_repeats = np.repeat(counts_sw, counts_sw)
                weights = self.sw_intensity[fsw_repeats] / count_repeats
                choices = np.argsort(-weights)[:n_pairs]
                p2 = fsw_repeats[choices]
                p1 = m_looking

        else:
            n_pairs = len(m_looking)
            weights = self.sw_intensity[active_fsw]
            choices = np.argsort(-weights)[:n_pairs]
            p2 = active_fsw.uids[choices]
            p1 = m_looking

        return p1, p2

    def end_pairs(self):
        people = self.sim.people

        self.edges.dur = self.edges.dur - 1  # Decrement the duration of each partnership, noting that dur is timesteps

        # Non-alive agents are removed
        alive_bools = people.alive[ss.uids(self.edges.p1)] & people.alive[ss.uids(self.edges.p2)]
        active = (self.edges.dur > 0) & alive_bools

        # If there's a prior partner network, add the newly dissolved partnerships to the prior partners
        if self.pars.recall_prior:
            prior_network = self.sim.networks.get('priorpartners')
            if prior_network is not None:
                # Get the uids of the partners that just ended
                ended_p1 = self.edges.p1[~active]
                ended_p2 = self.edges.p2[~active]
                durs = np.zeros_like(ended_p1, dtype=ss_float_)
                betas = np.zeros_like(ended_p1, dtype=ss_float_)

                # Add these to the prior partners network
                prior_network.append(p1=ended_p1, p2=ended_p2, dur=durs, beta=betas)

        # For gen pop contacts that are due to expire, decrement the partner count
        inactive_gp = ~active & (~self.edges.sw)

        p1_edges = self.edges.p1[inactive_gp]
        p2_edges = self.edges.p2[inactive_gp]
        edge_types = self.edges.edge_type[inactive_gp]
        onetimes = edge_types == self.edge_types['onetime']
        casuals  = edge_types == self.edge_types['casual']
        stables  = edge_types == self.edge_types['stable']
        sw = edge_types == self.edge_types['sw']

        self.partners[p1_edges] -= 1
        self.partners[p2_edges] -= 1

        self.onetime_partners[p1_edges[onetimes]] -= 1
        self.onetime_partners[p2_edges[onetimes]] -= 1
        self.casual_partners[p1_edges[casuals]] -= 1
        self.casual_partners[p2_edges[casuals]] -= 1
        self.stable_partners[p1_edges[stables]] -= 1
        self.stable_partners[p2_edges[stables]] -= 1
        self.sw_partners[p1_edges[sw]] -= 1
        self.sw_partners[p2_edges[sw]] -= 1
        # for a, b in zip(p1_edges[(casuals + stables)], p2_edges[(casuals + stables)]):
        #     pair = (min(a,b), max(a,b))
        #     self.relationship_durs[pair][-1]['dur'] = self.ti - self.relationship_durs[pair][-1]['start']

        # For all contacts that are due to expire, remove them from the contacts list
        if len(active) > 0:
            for k in self.meta_keys():
                self.edges[k] = (self.edges[k][active])

        return

    def net_beta(self, disease_beta=None, uids=None, disease=None):
        if uids is None: uids = Ellipsis
        p_condom = self.edges.condoms[uids]
        eff_condom = disease.pars.eff_condom
        p_trans_condom = (1 - disease_beta*(1-eff_condom))**(self.edges.acts[uids]*p_condom)
        p_trans_no_condom = (1 - disease_beta)**(self.edges.acts[uids]*(1-p_condom))
        p_trans = 1 - p_trans_condom * p_trans_no_condom
        result = p_trans * self.edges.beta[uids]
        return result

    def set_condom_use(self):
        """ Set condom use """
        if self.pars.condom_data is not None:
            if isinstance(self.pars.condom_data, dict):
                for rgm in range(self.pars.n_risk_groups):
                    for rgf in range(self.pars.n_risk_groups):
                        risk_pairing = (self.risk_group[self.p1] == rgm) & (self.risk_group[self.p2] == rgf)
                        self.edges.condoms[risk_pairing] = self.pars.condom_data[(rgm, rgf)]['simvals'][self.ti]
                self.edges.condoms[self.edges.sw] = self.pars.condom_data[('fsw','client')]['simvals'][self.ti]

            elif sc.isnumber(self.pars.condom_data):
                self.edges.condoms[:] = self.pars.condom_data

            else:
                raise Exception("Unknown condom data input type")

        return

    def count_partners(self):
        """ Count the number of partners each person has had over the past 3/12 months """
        self.lifetime_partners
        return

    def step(self):
        self.end_pairs()
        self.set_network_states(upper_age=self.t.dt_year)
        self.add_pairs()
        self.set_condom_use()
        self.count_partners()

        return


class PriorPartners(ss.DynamicNetwork):
    """
    Lightweight network for storing prior partners, for use in partner notification
    In this network, 'dur' refers to the duration of time since the relationship ended
    """
    def __init__(self, pars=None, name='priorpartners', **kwargs):
        super().__init__(name=name)
        self.define_pars(
            dur_recall=ss.years(1),  # How long to remember prior relationships
        )
        self.update_pars(pars=pars, **kwargs)
        return

    def step(self):
        self.end_pairs()
        self.edges.dur += 1  # Increment the duration since relationship ended
        return

    def end_pairs(self):
        people = self.sim.people
        max_dur = int(self.pars.dur_recall.value)
        active = (self.edges.dur < max_dur) & people.alive[self.edges.p1] & people.alive[self.edges.p2]
        for k in self.meta_keys():
            self.edges[k] = self.edges[k][active]
        return len(active)


class AgeMatchedMSM(StructuredSexual):

    def __init__(self, pars=None, **kwargs):
        super().__init__(name='msm')
        self.define_pars(
            msm_share=ss.bernoulli(p=0.015),
        )
        self.update_pars(pars=pars, **kwargs)

        return

    def set_network_states(self, upper_age=None):
        self.set_msm(upper_age=upper_age)
        return

    def set_msm(self, upper_age=None):
        _, m_uids = self._get_uids(upper_age=upper_age)
        self.participant[m_uids] = self.pars.msm_share.rvs(m_uids)
        return

    def match_pairs(self):
        """ Match males by age using sorting """
        ppl = self.sim.people

        # Find people eligible for a relationship
        active = self.over_debut
        underpartnered = self.partners < self.concurrency
        m_eligible = active & ppl.male & underpartnered
        m_looking = self.pars.p_pair_form.filter(m_eligible.uids)

        if len(m_looking) == 0:
            raise NoPartnersFound()

        # Match mairs by sorting the men looking for partners by age, then matching pairs by taking
        # 2 people at a time from the sorted list
        m_ages = ppl.age[m_looking]
        ind_m = np.argsort(m_ages)
        p1 = m_looking[ind_m][::2]
        p2 = m_looking[ind_m][1::2]
        maxlen = min(len(p1), len(p2))
        p1 = p1[:maxlen]
        p2 = p2[:maxlen]

        # Make sure everyone only appears once (?)
        if len(np.intersect1d(p1, p2)):
            errormsg = 'Some people appear in both p1 and p2'
            raise ValueError(errormsg)

        return p1, p2


class AgeApproxMSM(StructuredSexual):

    def __init__(self, **kwargs):
        super().__init__(name='msm', **kwargs)

    def match_pairs(self, ppl):
        """ Match pairs using age preferences """

        # Find people eligible for a relationship
        active = self.over_debut()
        underpartnered = self.partners < self.concurrency
        m_eligible = active & ppl.male & underpartnered
        m_looking = self.pars.p_pair_form.filter(m_eligible.uids)

        # Split the total number of males looking for partners into 2 groups
        # The first group will be matched with the second group
        group1 = m_looking[::2]
        group2 = m_looking[1::2]
        loc, scale = self.get_age_risk_pars(group1, self.pars.age_diff_pars)
        self.pars.age_diffs.set(loc=loc, scale=scale)
        age_gaps = self.pars.age_diffs.rvs(group1)
        desired_ages = ppl.age[group1] + age_gaps
        g2_ages = ppl.age[group2]
        ind_p1 = np.argsort(g2_ages)
        ind_p2 = np.argsort(desired_ages)
        p1 = m_eligible.uids[ind_p1]
        p2 = group2[ind_p2]
        maxlen = min(len(p1), len(p2))
        p1 = p1[:maxlen]
        p2 = p2[:maxlen]

        return p1, p2