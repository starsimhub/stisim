"""Heterosexual (MF) sexual network and its parameters."""
import numpy as np
import starsim as ss
from collections import defaultdict
from bisect import bisect_left

from .base import BaseNetwork, NoPartnersFound, ss_float

__all__ = ['MFPars', 'MFNetwork']


def _mf_states():
    """States specific to the heterosexual (MF) network (excludes shared)."""
    return [
        ss.FloatArr('risk_group'),
        ss.FloatArr('concurrency'),
        ss.FloatArr('partners', default=0),
        ss.FloatArr('partners_12', default=0),
        ss.FloatArr('lifetime_partners', default=0),
        ss.FloatArr('casual_partners', default=0),
        ss.FloatArr('stable_partners', default=0),
        ss.FloatArr('onetime_partners', default=0),
        ss.FloatArr('lifetime_casual_partners', default=0),
        ss.FloatArr('lifetime_stable_partners', default=0),
        ss.FloatArr('lifetime_onetime_partners', default=0),
    ]


class MFPars(ss.Pars):
    """MF-specific parameters: risk groups, concurrency, age-diff preferences, durations."""
    def __init__(self, **kwargs):
        super().__init__()

        # Risk groups
        self.n_risk_groups = 3
        self.f_age_group_bins = dict(
            teens=(0, 20),
            young=(20, 25),
            adult=(25, np.inf),
        )
        self.p_lo_risk = ss.bernoulli(p=0)
        self.p_hi_risk = ss.bernoulli(p=0)
        self.prop_f0 = 0.85
        self.prop_m0 = 0.8
        self.prop_f2 = 0.01
        self.prop_m2 = 0.02

        # Age difference preferences
        self.age_diff_pars = dict(
            teens=[(7, 3), (6, 3), (5, 1)],
            young=[(8, 3), (7, 3), (5, 2)],
            adult=[(8, 3), (7, 3), (5, 2)],
        )

        # Concurrency
        self.concurrency_dist = ss.poisson(lam=1)
        self.f0_conc = 0.0001
        self.f1_conc = 0.01
        self.f2_conc = 0.1
        self.m0_conc = 0.0001
        self.m1_conc = 0.2
        self.m2_conc = 0.5
        # Multiplier on FSW concurrency: values <1 mean that active sex workers
        # have fewer non-sex-work partners. Default 1.0 = no effect.
        self.fsw_mf_conc_mult = 1.0

        # Relationship initiation, stability, and duration
        self.p_pair_form = ss.bernoulli(p=0.5)              # Probability of a (stable) pair forming between two matched people
        self.match_dist = ss.bernoulli(p=0)                 # Placeholder value replaced by risk-group stratified values below
        self.p_matched_stable = [0.9, 0.5, 0]               # Probability of a stable pair forming between matched people (otherwise casual)
        self.p_mismatched_casual = [0.5, 0.5, 0.5]          # Probability of a casual pair forming between mismatched people (otherwise instantaneous)

        # Durations of stable and casual relationships
        self.stable_dur_pars = dict(
            teens=[
                [ss.years(100), ss.years(1)],
                [ss.years(8),   ss.years(2)],
                [ss.months(1e-4), ss.months(1e-4)],
            ],
            young=[
                [ss.years(100), ss.years(1)],
                [ss.years(10),  ss.years(3)],
                [ss.months(1e-4), ss.months(1e-4)],
            ],
            adult=[
                [ss.years(100), ss.years(1)],
                [ss.years(12),  ss.years(3)],
                [ss.months(1e-4), ss.months(1e-4)],
            ],
        )
        self.casual_dur_pars = dict(
            teens=[[ss.years(1), ss.years(3)]]*3,
            young=[[ss.years(1), ss.years(3)]]*3,
            adult=[[ss.years(1), ss.years(3)]]*3,
        )

        # Distributions populated by other methods at runtime
        self.age_diffs = ss.normal()
        self.dur_dist = ss.lognorm_ex()

        # relationship search taper for older women -- used only for match method closest_age_tapered_seeking
        self.p_actually_looking = ss.bernoulli(p=0)  # Placeholder to be replaced by agent-based calculation per-timestep
        self._f_partnership_taper_offset = 0  # Placeholder: mean_age_gap_target + 3*sd, maximum of all age groupings. Set in init_post()
        self.f_partnership_taper_cut = 55  # max age, over which females no longer search for new relationships

        # Pair-formation algorithm: string key in matchers.MATCHERS, or a callable.
        self.match_method = 'closest_age_tapered_seeking'

        self.update(kwargs)
        return


class MFNetwork(BaseNetwork):
    """Heterosexual contact network with risk groups, concurrency, and age preferences.

    Agents are assigned to one of three risk groups (low, medium, high) that
    govern partnership formation, concurrency, and relationship duration.
    Partnerships are formed each timestep by matching under-partnered agents
    using age preferences. Edges are typed as ``stable``, ``casual``, or
    ``onetime``. Does not model sex work — combine with :class:`SWNetwork`
    or use :class:`StructuredSexual` for an MF + SW bundle.

    Args:
        pars (dict): Parameter overrides (see :class:`MFPars` for defaults).
        condom_data: Optional condom-use data (DataFrame or dict).
        name (str): Network name (default: auto-assigned).
        **kwargs: Additional parameter overrides forwarded to ``update_pars``.
    """

    def __init__(self, pars=None, condom_data=None, name=None, **kwargs):
        super().__init__(name=name)
        self.define_pars(**MFPars())
        self.update_pars(pars, **kwargs)
        self.edge_types = {'stable': 0, 'casual': 1, 'onetime': 2}
        self.define_states(*_mf_states())
        self.relationship_durs = defaultdict(list)
        return

    def init_post(self):
        """
        This is intended to auto-set the female partnership looking rate for pairing based on input age gap data.
        It is needed if using the closest_age_tapered_seeking matching algorithm (otherwise ignored).
        """
        super().init_post()
        uids_born = (self.sim.people.age > 0).uids
        loc, scale = self.get_age_risk_pars(uids_born, self.pars.age_diff_pars)
        self.pars._f_partnership_taper_offset = max(loc + 3 * scale)  # mean age + 3*sd
        return

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

    def set_network_states(self, upper_age=None):
        super().set_network_states(upper_age=upper_age)  # set_debut
        self.set_risk_groups(upper_age=upper_age)
        self.set_concurrency(upper_age=upper_age)
        return

    def set_risk_groups(self, upper_age=None):
        """ Assign each person to a risk group """
        ppl = self.sim.people
        uids = self._get_uids(upper_age=upper_age, by_sex=False)

        p_lo = np.full(len(uids), fill_value=np.nan, dtype=ss_float)
        p_lo[ppl.female[uids]] = self.pars.prop_f0
        p_lo[ppl.male[uids]] = self.pars.prop_m0
        self.pars.p_lo_risk.set(p=p_lo)
        lo_risk, hi_med_risk = self.pars.p_lo_risk.split(uids)

        p_hi = np.full(len(hi_med_risk), fill_value=np.nan, dtype=ss_float)
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

        mu = np.full(uids.shape, fill_value=np.nan, dtype=ss_float)
        for rg in range(self.pars.n_risk_groups):
            f_conc = self.pars[f'f{rg}_conc']
            m_conc = self.pars[f'm{rg}_conc']
            in_risk_group = self.risk_group == rg
            in_group = in_risk_group & in_age_lim
            f_in = (people.female & in_group)[uids]
            m_in = (people.male   & in_group)[uids]
            if f_in.any(): mu[f_in] = f_conc
            if m_in.any(): mu[m_in] = m_conc

        # Note: the FSW-specific multiplier (fsw_mf_conc_mult) is applied in
        # StructuredSexual.set_network_states post-hoc, because set_sex_work
        # populates self.fsw AFTER set_concurrency runs in the MF flow.

        dist = self.pars.concurrency_dist
        if isinstance(dist, ss.nbinom):
            n = dist.pars['n']
            p = n / (n + mu)
            dist.set(n=n, p=p)
        else:
            dist.set(lam=mu)
        self.concurrency[uids] = dist.rvs(uids) + 1

        return

    def _get_eligible(self):
        """Return (f_looking, m_eligible) ss.uids. Raises NoPartnersFound if either is empty."""
        ppl = self.sim.people
        active = self.over_debut
        underpartnered = self.partners < self.concurrency
        f_eligible = active & ppl.female & underpartnered
        m_eligible = active & ppl.male & underpartnered
        f_looking = self.pars.p_pair_form.filter(f_eligible.uids)
        if len(f_looking) == 0 or m_eligible.count() == 0:
            raise NoPartnersFound()
        return f_looking, m_eligible

    def _sample_desired_ages(self, f_looking):
        """Sample desired male partner ages for the given f_looking uids."""
        loc, scale = self.get_age_risk_pars(f_looking, self.pars.age_diff_pars)
        self.pars.age_diffs.set(loc=loc, scale=scale)
        age_gaps = self.pars.age_diffs.rvs(f_looking)
        return self.sim.people.age[f_looking] + age_gaps

    def match_pairs(self):
        """Dispatch to the matcher named by ``pars.match_method``.

        ``match_method`` is either a string key in ``matchers.MATCHERS``
        or a callable ``f(net) -> (p1, p2)``. See ``matchers.py``.
        """
        from .matchers import MATCHERS
        m = self.pars.match_method
        if callable(m):
            return m(self)
        return MATCHERS[m](self)

    def add_pairs(self):
        """ Match and add stable/casual/onetime partnerships for this timestep. Assigns relationship type, duration, and acts based on risk group and age, and updates partner counts. """
        ppl = self.sim.people

        try:
            p1, p2 = self.match_pairs()
        except NoPartnersFound:
            return

        matched_risk = (self.risk_group[p1] == self.risk_group[p2])
        mismatched_risk = (self.risk_group[p1] != self.risk_group[p2])

        # Set the probability of forming a partnership
        p_match = np.full(len(p1), fill_value=np.nan, dtype=ss_float)
        for rg in range(self.pars.n_risk_groups):
            p_match[matched_risk & (self.risk_group[p1] == rg)] = self.pars.p_matched_stable[rg]
            p_match[mismatched_risk & (self.risk_group[p2] == rg)] = self.pars.p_mismatched_casual[rg]
        self.pars.match_dist.set(p=p_match)
        matches = self.pars.match_dist.rvs(p2)

        stable = matches & matched_risk
        casual = matches & mismatched_risk
        any_match = stable | casual

        match_count = len(p1)

        beta = np.ones(match_count, dtype=ss_float)
        condoms = np.zeros(match_count, dtype=ss_float)
        acts = (self.pars.acts.rvs(p2)).astype(int)
        dur = np.full(match_count, fill_value=1)  # Measured in timesteps
        age_p1 = ppl.age[p1]
        age_p2 = ppl.age[p2]
        edge_types = np.full(match_count, dtype=ss_float, fill_value=np.nan)
        edge_types[stable] = self.edge_types['stable']
        edge_types[casual] = self.edge_types['casual']

        # Set duration
        dur_mean = np.full(match_count, fill_value=np.nan, dtype=ss_float)
        dur_std = np.full(match_count, fill_value=np.nan, dtype=ss_float)
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

        for (a, b, reldur, etype) in zip(p1[relationships], p2[relationships], dur[relationships], edge_types[relationships]):
            pair = (min(a,b), max(a,b))
            self.relationship_durs[pair].append({'start': self.ti, 'dur': reldur, 'edge_type': int(etype)}) # set dur to intended duration. When the relationship actually ends, this will be updated

        self.append(p1=p1, p2=p2, beta=beta, condoms=condoms, dur=dur, acts=acts, age_p1=age_p1, age_p2=age_p2, edge_type=edge_types)

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

    def _on_edge_dissolution(self, active):
        """Record dissolved partnerships and decrement partner counts."""
        self._record_prior_partners(active)
        self._decrement_partners(active)

    def _record_prior_partners(self, active):
        """If ``recall_prior``, push dissolved partnerships into the PriorPartners network."""
        if not self.pars.recall_prior:
            return
        prior_network = self.sim.networks.get('priorpartners')
        if prior_network is None:
            return
        ended_p1 = self.edges.p1[~active]
        ended_p2 = self.edges.p2[~active]
        durs = np.zeros_like(ended_p1, dtype=ss_float)
        betas = np.zeros_like(ended_p1, dtype=ss_float)
        prior_network.append(p1=ended_p1, p2=ended_p2, dur=durs, beta=betas)

    def _decrement_partners(self, active):
        """Decrement partner counters for expired non-SW edges.

        SW edges (when present) are skipped because ``add_pairs`` doesn't
        increment ``partners`` for them — only ``lifetime_sw_partners``.
        """
        inactive = ~active
        if 'sw' in self.edge_types:
            inactive = inactive & (self.edges.edge_type != self.edge_types['sw'])
        p1e = self.edges.p1[inactive]
        p2e = self.edges.p2[inactive]
        edge_types = self.edges.edge_type[inactive]
        self.partners[p1e] -= 1
        self.partners[p2e] -= 1
        for key in ('stable', 'casual', 'onetime'):
            if key not in self.edge_types:
                continue
            mask = edge_types == self.edge_types[key]
            if not mask.any():
                continue
            getattr(self, f'{key}_partners')[p1e[mask]] -= 1
            getattr(self, f'{key}_partners')[p2e[mask]] -= 1
