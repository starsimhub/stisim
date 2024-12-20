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
import scipy.optimize as spo
import scipy.spatial as spsp

ss_float_ = ss.dtypes.float

# Specify all externally visible functions this file defines; see also more definitions below
__all__ = ['StructuredSexual', 'FastStructuredSexual', 'AgeMatchedMSM', 'AgeApproxMSM']


class NoPartnersFound(Exception):
    # Raise this exception if the matching algorithm wasn't able to match any partners
    pass


class StructuredSexual(ss.SexualNetwork):
    """
    Structured sexual network
    """

    def __init__(self, pars=None, key_dict=None, condom_data=None, name=None, **kwargs):

        key_dict = sc.mergedicts({
            'sw': bool,
            'condoms': ss_float_,
            'age_p1': ss_float_,
            'age_p2': ss_float_,
        }, key_dict)

        super().__init__(key_dict=key_dict, name=name)

        self.define_pars(
            # Settings - generally shouldn't be adjusted
            unit='month',
            store_register=False,
            n_risk_groups=3,
            f_age_group_bins=dict(  # For separating women into age groups: teens, young women, adult women
                teens=(0, 20),
                young=(20, 25),
                adult=(25, np.inf),
            ),

            # Age of sexual debut
            debut=ss.lognorm_ex(20, 3),
            debut_pars_f=[20, 3],
            debut_pars_m=[21, 3],

            # Risk groups
            p_lo_risk=ss.bernoulli(p=0),
            p_hi_risk=ss.bernoulli(p=0),
            prop_f0=0.85,
            prop_m0=0.8,
            prop_f2=0.01,
            prop_m2=0.02,

            # Age difference preferences
            age_diff_pars=dict(
                teens=[(7, 3), (6, 3), (5, 1)],  # (mu,stdev) for levels 0, 1, 2
                young=[(8, 3), (7, 3), (5, 2)],
                adult=[(8, 3), (7, 3), (5, 2)],
            ),

            # Concurrency preferences
            concurrency_dist=ss.poisson(lam=1),
            f0_conc=0.0001,
            f1_conc=0.01,
            f2_conc=0.1,
            m0_conc=0.0001,
            m1_conc=0.2,
            m2_conc=0.5,

            # Relationship initiation, stability, and duration
            p_pair_form=ss.bernoulli(p=0.5),  # Probability of a (stable) pair forming between two matched people
            match_dist=ss.bernoulli(p=0),  # Placeholder value replaced by risk-group stratified values below
            p_matched_stable=[0.9, 0.5, 0],  # Probability of a stable pair forming between matched people (otherwise casual)
            p_mismatched_casual=[0.5, 0.5, 0.5],  # Probability of a casual pair forming between mismatched people (otherwise instantanous)

            # Durations of stable and casual relationships
            stable_dur_pars=dict(
                teens=[
                    # (mu,stdev) for levels 0, 1, 2
                    [ss.dur(100, 'year'),  ss.dur(1, 'year')],
                    [ss.dur(8, 'year'),  ss.dur(2, 'year')],
                    [ss.dur(1e-4, 'month'), ss.dur(1e-4, 'month')]
                ],
                young=[
                    [ss.dur(100, 'year'),  ss.dur(1, 'year')],
                    [ss.dur(10, 'year'),  ss.dur(3, 'year')],
                    [ss.dur(1e-4, 'month'), ss.dur(1e-4, 'month')]
                ],
                adult=[
                    [ss.dur(100, 'year'),  ss.dur(1, 'year')],
                    [ss.dur(12, 'year'),  ss.dur(3, 'year')],
                    [ss.dur(1e-4, 'month'), ss.dur(1e-4, 'month')]
                ],
            ),
            casual_dur_pars=dict(
                teens=[[ss.dur(1, 'year'), ss.dur(3, 'year')]]*3,
                young=[[ss.dur(1, 'year'), ss.dur(3, 'year')]]*3,
                adult=[[ss.dur(1, 'year'), ss.dur(3, 'year')]]*3,
            ),

            # Acts
            acts=ss.lognorm_ex(ss.peryear(80), ss.peryear(30)),  # Coital acts/year

            # Sex work parameters
            fsw_shares=ss.bernoulli(p=0.05),
            client_shares=ss.bernoulli(p=0.12),
            sw_seeking_rate=ss.rate(1, 'month'),  # Monthly rate at which clients seek FSWs (1 new SW partner / month)
            sw_seeking_dist=ss.bernoulli(p=0.5),  # Placeholder value replaced by dt-adjusted sw_seeking_rate
            sw_beta=1,  
            sw_intensity=ss.random(),  # At each time step, FSW may work with varying intensity

            # Distributions derived from parameters above - don't adjust
            age_diffs=ss.normal(),
            dur_dist=ss.lognorm_ex(),
        )

        self.update_pars(pars=pars, **kwargs)

        # Set condom use
        self.condom_data = None
        if condom_data is not None:
            self.condom_data = self.process_condom_data(condom_data)

        # Store register
        if self.pars.store_register:
            self.breakup_register = [[]]*12

        # Add states
        self.define_states(
            ss.BoolArr('participant', default=True),
            ss.FloatArr('risk_group'),  # Which risk group an agent belongs to
            ss.BoolArr('fsw'),  # Whether an agent is a female sex worker
            ss.BoolArr('client'),  # Whether an agent is a client of sex workers
            ss.FloatArr('concurrency'),  # Preferred number of concurrent partners
            ss.FloatArr('partners', default=0),  # Actual number of concurrent partners
            ss.FloatArr('partners_12', default=0),  # Number of partners over the past 12m
            ss.FloatArr('lifetime_partners', default=0),  # Lifetime total number of partners
            ss.FloatArr('sw_intensity'),  # Intensity of sex work
        )

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
        loc = np.full(uids.shape, fill_value=np.nan, dtype=ss_float_)
        scale = np.full(uids.shape, fill_value=np.nan, dtype=ss_float_)
        for a_label, (age_lower, age_upper) in self.pars.f_age_group_bins.items():
            for rg in range(self.pars.n_risk_groups):
                in_risk_group = (self.sim.people.age[uids] >= age_lower) & (self.sim.people.age[uids] < age_upper) & (self.risk_group[uids] == rg)
                loc[in_risk_group] = par[a_label][rg][0]
                scale[in_risk_group] = par[a_label][rg][1]
        if np.isnan(scale).any() or np.isnan(loc).any():
            errormsg = 'Invalid entries for age difference preferences.'
            raise ValueError(errormsg)
        return loc, scale

    def init_pre(self, sim):
        super().init_pre(sim)
        if self.condom_data is not None:
            if isinstance(self.condom_data, dict):
                for rgtuple, valdict in self.condom_data.items():
                    self.condom_data[rgtuple]['simvals'] = sc.smoothinterp(sim.timevec, valdict['year'], valdict['val'])
        # self.init_results()
        return

    def init_results(self):
        self.define_results(
            ss.Result('share_active', dtype=float, scale=False),
            ss.Result('partners_f_mean', dtype=float, scale=False),
            ss.Result('partners_m_mean', dtype=float, scale=False),
        )
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

    def over_debut(self):
        return self.sim.people.age > self.debut

    def _get_uids(self, upper_age=None, by_sex=True):
        people = self.sim.people
        if upper_age is None: upper_age = 1000
        within_age = people.age < upper_age
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

    def match_pairs(self, ppl):
        """
        Match pairs by age
        """

        # Find people eligible for a relationship
        active = self.over_debut()
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
        dist_mat = spsp.distance_matrix(m_ages[:, np.newaxis], desired_ages[:, np.newaxis])
        ind_m, ind_f = spo.linear_sum_assignment(dist_mat)
        p1 = m_eligible.uids[ind_m]
        p2 = f_looking[ind_f]

        return p1, p2

    def add_pairs(self, ti=None):
        """ Add pairs """
        ppl = self.sim.people
        dt = self.t.dt

        # Obtain new pairs
        try:
            p1_gp, p2_gp = self.match_pairs(ppl)
            p1_sw, p2_sw = self.match_sex_workers(ppl)
        except NoPartnersFound:
            return

        p1 = p1_gp.concat(p1_sw)
        p2 = p2_gp.concat(p2_sw)
        sw = np.array([False]*len(p1_gp) + [True]*len(p1_sw))

        # Initialize beta, acts, duration
        beta = np.ones(len(p2), dtype=ss_float_)
        condoms = np.zeros(len(p2), dtype=ss_float_)  # FILLED IN LATER
        acts = (self.pars.acts.rvs(p2)).astype(int)  # Number of acts per timestep - does not depend on commitment/risk group
        dur = np.full(len(p2), dtype=ss_float_, fill_value=dt)  # Default duration is dt, replaced for stable matches
        age_p1 = ppl.age[p1]
        age_p2 = ppl.age[p2]

        # Determine whether the pair have matched risk profiles
        # Partners with mismatched risk profiles may still form a casual partnership
        matched_risk = (self.risk_group[p1] == self.risk_group[p2]) & ~sw
        mismatched_risk = (self.risk_group[p1] != self.risk_group[p2]) & ~sw

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

        # Set duration
        dur_mean = np.full(sum(any_match), fill_value=np.nan, dtype=ss_float_)
        dur_std = np.full(sum(any_match), fill_value=np.nan, dtype=ss_float_)
        for which, bools in {'stable': stable, 'casual': casual}.items():
            if bools.any():
                uids = p2[bools]
                mean, std = self.get_age_risk_pars(uids, self.pars[f'{which}_dur_pars'])
                inds = bools[any_match].nonzero()[-1]
                dur_mean[inds] = mean
                dur_std[inds] = std
        self.pars.dur_dist.set(mean=dur_mean, std=dur_std)
        dur[any_match] = self.pars.dur_dist.rvs(p2[any_match])

        self.append(p1=p1, p2=p2, beta=beta, condoms=condoms, dur=dur, acts=acts, sw=sw, age_p1=age_p1, age_p2=age_p2)

        # Checks
        if (self.sim.people.female[p1].any() or self.sim.people.male[p2].any()) and (self.name == 'structuredsexual'):
            errormsg = 'Same-sex pairings should not be possible in this network'
            raise ValueError(errormsg)
        if len(p1) != len(p2):
            errormsg = 'Unequal lengths in edge list'
            raise ValueError(errormsg)

        # Add partner counts, not including SW partners
        unique_p1, counts_p1 = np.unique(p1_gp, return_counts=True)
        unique_p2, counts_p2 = np.unique(p2_gp, return_counts=True)
        self.partners[unique_p1] += counts_p1
        self.partners[unique_p2] += counts_p2
        self.lifetime_partners[unique_p1] += counts_p1
        self.lifetime_partners[unique_p2] += counts_p2

        return

    def match_sex_workers(self, ppl):
        """ Match sex workers to clients """

        # Find people eligible for a relationship
        active = self.over_debut()
        active_fsw = active & self.fsw
        active_clients = active & self.client
        self.sw_intensity[active_fsw.uids] = self.pars.sw_intensity.rvs(active_fsw.uids)

        # Find clients who will seek FSW
        self.pars.sw_seeking_dist.pars.p = np.clip(self.pars.sw_seeking_rate, 0, 1)
        m_looking = self.pars.sw_seeking_dist.filter(active_clients.uids)

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

        # Update the breakup register
        if self.pars.store_register:
            over_12m = self.breakup_register[11]
            if len(over_12m):
                u, c = np.unique(over_12m, return_counts=True)
                self.partners_12[u] -= c
            self.breakup_register = self.breakup_register[:11]  # Forget partners from >12m ago
            just_ended = (self.edges.dur == 0) & alive_bools
            je1 = ss.uids(self.edges.p1[just_ended])
            je2 = ss.uids(self.edges.p2[just_ended])
            je_uids = je1.concat(je2)
            self.breakup_register.insert(0, je1.concat(je2))
            if len(je_uids):
                u, c = np.unique(je_uids, return_counts=True)
                self.partners_12[u] += c

        # For gen pop contacts that are due to expire, decrement the partner count
        inactive_gp = ~active & (~self.edges.sw)
        self.partners[ss.uids(self.edges.p1[inactive_gp])] -= 1
        self.partners[ss.uids(self.edges.p2[inactive_gp])] -= 1

        # For all contacts that are due to expire, remove them from the contacts list
        if len(active) > 0:
            for k in self.meta_keys():
                self.edges[k] = (self.edges[k][active])

        return

    def update_results(self):
        ti = self.ti
        self.results.share_active[ti] = len(self.active(self.sim.people).uids)/len(self.sim.people)

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
        if self.condom_data is not None:
            if isinstance(self.condom_data, dict):
                for rgm in range(self.pars.n_risk_groups):
                    for rgf in range(self.pars.n_risk_groups):
                        risk_pairing = (self.risk_group[self.p1] == rgm) & (self.risk_group[self.p2] == rgf)
                        self.edges.condoms[risk_pairing] = self.condom_data[(rgm, rgf)]['simvals'][self.ti]
                self.edges.condoms[self.edges.sw] = self.condom_data[('fsw','client')]['simvals'][self.ti]

            elif sc.isnumber(self.condom_data):
                self.edges.condoms[:] = self.condom_data

            else:
                raise Exception("Unknown condom data input type")

        return

    def count_partners(self):
        """ Count the number of partners each person has had over the past 3/12 months """
        self.lifetime_partners
        return

    def step(self):
        self.end_pairs()
        self.set_network_states(upper_age=self.t.dt)
        self.add_pairs()
        self.set_condom_use()

        self.count_partners()

        return


class FastStructuredSexual(StructuredSexual):

    def __init__(self, **kwargs):
        super().__init__(name='structuredsexual', **kwargs)

    def match_pairs(self, ppl):
        """
        Match pairs by age, using sorting rather than the linear sum assignment
        """

        # Find people eligible for a relationship
        active = self.over_debut()
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
        ind_m = np.argsort(m_ages)  # Use sort instead of linear_sum_agreement
        ind_f = np.argsort(desired_ages)
        p1 = m_eligible.uids[ind_m]
        p2 = f_looking[ind_f]
        maxlen = min(len(p1), len(p2))
        p1 = p1[:maxlen]
        p2 = p2[:maxlen]

        return p1, p2


class AgeMatchedMSM(StructuredSexual):

    def __init__(self, **kwargs):
        super().__init__(name='msm', **kwargs)

    def match_pairs(self, ppl):
        """ Match males by age using sorting """

        # Find people eligible for a relationship
        active = self.over_debut()
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
        """ Match"""

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