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

from starsim import uids

ss_float = ss.dtypes.float
ss_int = ss.dtypes.int

# Specify all externally visible functions this file defines; see also more definitions below
__all__ = ['NetworkPars', 'BasePars', 'MFPars', 'SWPars',
           'StructuredSexual', 'BaseNetwork', 'MFNetwork', 'SWNetwork',
           'PriorPartners', 'AgeMatchedMSM', 'AgeApproxMSM']


class NoPartnersFound(Exception):
    # Raise this exception if the matching algorithm wasn't able to match any partners
    pass


# %% State helpers — lists of agent-level states fed to ``define_states``

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


def _sw_states():
    """States specific to the sex-work (SW) network (excludes shared)."""
    return [
        ss.BoolArr('ever_fsw'),                # Lifetime SW fate (female)
        ss.BoolArr('ever_client'),             # Lifetime SW fate (male)
        ss.FloatArr('age_sw_start'),           # Age of SW entry (NaN unless ever_fsw)
        ss.FloatArr('dur_sw'),                 # SW duration in years (NaN unless ever_fsw)
        ss.FloatArr('age_client_start'),
        ss.FloatArr('dur_client'),
        ss.FloatArr('sw_intensity'),
        ss.FloatArr('sw_partners', default=0),
        ss.FloatArr('lifetime_sw_partners', default=0),
    ]


# %% Parameter classes
#
# ``BasePars`` defaults are layered into ``BaseNetwork`` via ``define_pars``.
# ``MFPars``/``SWPars`` are layered on top by their respective subclasses.
# ``NetworkPars`` is the union — used by ``sti.Sim`` to route user kwargs
# into the right network module.

class BasePars(ss.Pars):
    """Defaults shared by every sexual network (debut, acts, condom data, recall_prior)."""
    def __init__(self, **kwargs):
        super().__init__()
        self.recall_prior = False
        self.debut_f = ss.lognorm_ex(20, 3)
        self.debut_m = ss.lognorm_ex(21, 3)
        self.acts = ss.lognorm_ex(ss.freqperyear(80), ss.freqperyear(30))
        self.condom_data = None
        self.condom_smoothness = None  # Smoothness for sc.smoothinterp on time-varying condom_data
        self.update(kwargs)
        return


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

        self.update(kwargs)
        return


class SWPars(ss.Pars):
    """SW-specific parameters: shares, seeking rates, intensity, age-window distributions.

    The ``fsw``/``client`` properties are derived from the ``ever_*`` flag
    plus a per-agent age window. For each agent flagged ``ever_*``,
    ``set_sex_work`` draws ``age_*_start`` and ``dur_*`` from these
    distributions; the agent counts as currently SW iff
    ``age_*_start <= age < age_*_start + dur_*``. Setting any ``*_dist`` to
    ``None`` makes ``set_sex_work`` fall back to ``debut`` (for missing
    entry-age) or ``inf`` (for missing duration).

    Distribution defaults are placeholders — tune to the modeling context.
    """
    def __init__(self, **kwargs):
        super().__init__()

        self.fsw_shares = ss.bernoulli(p=0.05)         # Lifetime proportion of females ever FSW
        self.client_shares = ss.bernoulli(p=0.12)      # Lifetime proportion of males ever client
        self.sw_seeking_rate = ss.probpermonth(1.0)    # Monthly rate at which clients seek FSWs
        self.sw_seeking_dist = ss.bernoulli(p=0.5)     # Placeholder; replaced by dt-adjusted sw_seeking_rate
        self.sw_beta = 1
        self.sw_intensity = ss.random()                # FSW may work with varying intensity each step

        # Per-agent age window — placeholders; tune to context
        self.age_sw_start = ss.normal(loc=20, scale=3)
        self.dur_sw = ss.lognorm_ex(mean=5, std=3)
        self.age_client_start = ss.normal(loc=25, scale=5)
        self.dur_client = ss.lognorm_ex(mean=10, std=5)

        self.update(kwargs)
        return


class NetworkPars(ss.Pars):
    """Combined defaults — union of base + MF + SW pars.

    Used by :class:`sti.Sim` to filter network kwargs into the right bucket.

    Args:
        **kwargs: Any parameter name-value pairs to override defaults.
    """
    def __init__(self, **kwargs):
        super().__init__()
        for src in (BasePars(), MFPars(), SWPars()):
            self.update(dict(src), create=True)
        self.update(kwargs)
        return


# %% Network base class
#
# ``BaseNetwork`` provides edge-meta layout, shared pars/states (debut,
# participant), condom-data processing, end-of-step accounting, and the
# ``step`` orchestrator. Concrete networks (``MFNetwork``, ``SWNetwork``)
# extend it. ``StructuredSexual`` layers SW state on top of ``MFNetwork``.

class BaseNetwork(ss.SexualNetwork):
    """Shared infrastructure for the heterosexual and sex-work networks.

    Provides the common edge meta layout, condom-data processing, debut
    setting, end-of-step edge accounting, and the orchestrating ``step``
    method. Subclasses (``MFNetwork``, ``SWNetwork``) define their own
    parameters, states, edge types, ``set_network_states``, ``add_pairs``,
    and ``set_condom_use``. The MSM subclasses inherit this via
    ``MFNetwork``.
    """

    def __init__(self, name=None, **kwargs):
        super().__init__(name=name)
        self.meta.condoms = ss_float
        self.meta.age_p1 = ss_float
        self.meta.age_p2 = ss_float
        self.meta.edge_type = ss_float
        self.define_pars(**BasePars())
        self.define_states(
            ss.BoolArr('participant', default=True),
            ss.FloatArr('debut', default=0),
            reset=True,
        )

    @staticmethod
    def process_condom_data(condom_data):
        if sc.isnumber(condom_data) or isinstance(condom_data, dict):
            return condom_data
        if isinstance(condom_data, pd.DataFrame):
            df = condom_data.melt(id_vars=['partnership'])
            dd = dict()
            for pcombo in df.partnership.unique():
                key = tuple(map(int, pcombo[1:-1].split(','))) if pcombo != '(fsw,client)' else ('fsw', 'client')
                thisdf = df.loc[df.partnership == pcombo]
                dd[key] = dict()
                dd[key]['year'] = thisdf.variable.values.astype(int)
                dd[key]['val'] = thisdf.value.values
            return dd
        raise Exception(f"Unknown condom data input type: {type(condom_data)}")

    def init_pre(self, sim):
        super().init_pre(sim)
        if self.pars.recall_prior:
            isprior = [isinstance(nw, PriorPartners) for nw in self.sim.networks.values()]
            if not any(isprior):
                raise ValueError('PriorPartners network is required if recall_prior is True.')
        # Convert any DataFrame condom_data to dict; smooth-interp time-varying values
        cd = self.pars.condom_data
        if cd is not None and not sc.isnumber(cd):
            cd = self.process_condom_data(cd)
            self.pars.condom_data = cd
            for rgtuple, valdict in cd.items():
                valdict['simvals'] = sc.smoothinterp(self.t.yearvec, valdict['year'], valdict['val'],
                                                     smoothness=self.pars.condom_smoothness)
        return

    def init_post(self):
        super().init_post(add_pairs=False)
        self.set_network_states()
        return

    def set_network_states(self, upper_age=None):
        """Default: only set debut. Subclasses extend with risk groups, sex work, etc."""
        self.set_debut(upper_age=upper_age)

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
        return within_age.uids

    def set_debut(self, upper_age=None):
        uids = self._get_uids(upper_age=upper_age, by_sex=False)
        female = self.sim.people.female[uids]
        f_uids = uids[female]
        m_uids = uids[~female]
        self.debut[f_uids] = self.pars.debut_f.rvs(f_uids)
        self.debut[m_uids] = self.pars.debut_m.rvs(m_uids)
        return

    def net_beta(self, disease_beta=None, uids=None, disease=None):
        if uids is None: uids = Ellipsis
        p_condom = self.edges.condoms[uids]
        eff_condom = disease.pars.eff_condom
        p_trans_condom = (1 - disease_beta*(1-eff_condom))**(self.edges.acts[uids]*p_condom)
        p_trans_no_condom = (1 - disease_beta)**(self.edges.acts[uids]*(1-p_condom))
        return (1 - p_trans_condom * p_trans_no_condom) * self.edges.beta[uids]

    def step(self):
        self.end_pairs()
        self.set_network_states(upper_age=self.t.dt_year)
        self.add_pairs()
        self.set_condom_use()
        return

    def add_pairs(self):
        """Subclasses override."""
        pass

    def set_condom_use(self):
        """Apply ``self.pars.condom_data`` to per-edge condom values.

        ``condom_data`` may be:
        - ``None``: no-op.
        - A scalar: applied uniformly to all edges.
        - A dict, with keys interpreted by what the network exposes:
          - ``(rgm, rgf)`` (ints, MF risk-group pair) — applied to edges
            where ``risk_group[p1]==rgm`` and ``risk_group[p2]==rgf``.
            Skipped on networks without a ``risk_group`` state.
          - ``('fsw', 'client')`` — applied to all SW edges. Skipped on
            networks that don't define a ``'sw'`` edge type.
        """
        cd = self.pars.condom_data
        if cd is None:
            return
        if sc.isnumber(cd):
            self.edges.condoms[:] = cd
            return
        if not isinstance(cd, dict):
            raise Exception("Unknown condom data input type")

        for key, valdict in cd.items():
            val = valdict['simvals'][self.ti]
            if key == ('fsw', 'client'):
                if 'sw' in self.edge_types:
                    sw_mask = self.edges.edge_type == self.edge_types['sw']
                    self.edges.condoms[sw_mask] = val
            elif hasattr(self, 'risk_group'):
                rgm, rgf = key
                pair_mask = (self.risk_group[self.p1] == rgm) & (self.risk_group[self.p2] == rgf)
                self.edges.condoms[pair_mask] = val
        return

    def end_pairs(self):
        people = self.sim.people
        self.edges.dur = self.edges.dur - 1
        alive = people.alive[ss.uids(self.edges.p1)] & people.alive[ss.uids(self.edges.p2)]
        active = (self.edges.dur > 0) & alive
        self._on_edge_dissolution(active)
        # Drop expired edges
        for k in self.meta_keys():
            self.edges[k] = self.edges[k][active]
        return

    def _on_edge_dissolution(self, active):
        """Subclass hook — runs after the active mask is computed and before
        expired edges are removed. Default is a no-op."""
        pass


# %% Network classes — concrete subclasses of BaseNetwork

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

        lam = np.full(uids.shape, fill_value=np.nan, dtype=ss_float)
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

    # def bin_agents_by_age(self, agents_df, bin_width):
    #     ppl = self.sim.people
    #
    #     min_age = np.floor(ppl.age[agents_df].min())
    #     max_age = np.ceil(ppl.age[agents_df].max())
    #     max_age_bin_tolerated = max_age + bin_width/2  # the bin covering the extreme high age potential (likely empty)
    #     bin_offset = bin_width / 2
    #
    #     # bin_centers = [age + bin_offset for age]
    #     # for bin_age in range(min_age + bin_offset, )
    #     age_bins = {}
    #     for age in np.arange(min_age + bin_offset, max_age_bin_tolerated, bin_width):
    #         age_bins[age] = agents_df & (ppl.age >= age - bin_offset) & (ppl.age < age + bin_offset)
    #     return age_bins

    def bin_agents_by_age(self, agent_uids, bin_width):
        """
        return an {age: uids} keyed dict, where age is the center of bins bin_width years wide and uids is the
        uids from agent_uids matching the keyed age bin.
        """
        ppl = self.sim.people

        min_age = np.floor(ppl.age[agent_uids].min())
        max_age = np.ceil(ppl.age[agent_uids].max())
        max_age_bin_tolerated = max_age + bin_width/2  # the bin covering the extreme high age potential (likely empty)
        bin_offset = bin_width / 2

        age_bins = {}
        for age in np.arange(min_age + bin_offset, max_age_bin_tolerated, bin_width):
            age_bins[age] = agent_uids & (ppl.age >= age - bin_offset) & (ppl.age < age + bin_offset)
        return age_bins

    def bin_female_uids_by_target_age(self, female_uids, target_age_gaps, male_age_bins, bin_width):
        """female_uids and target_age_gaps are assumed to be co-sorted"""
        ppl = self.sim.people
        target_age_bins = {}
        bin_offset = bin_width / 2
        target_ages = ppl.age[female_uids] + target_age_gaps
        all_matched_indicies = []
        for age in male_age_bins:
            indicies_matching = np.where( (target_ages >= age - bin_offset) & (target_ages < age + bin_offset) )[0]
            uids_matching = female_uids[indicies_matching]
            target_age_bins[age] = uids_matching
            # print(f"target: {age} uids: {uids_matching} ages: {ppl.age[uids_matching]} gap: {target_age_gaps[indicies_matching]} selected_target: {ppl.age[uids_matching] + target_age_gaps[indicies_matching]}")
            all_matched_indicies.extend(indicies_matching)
        return target_age_bins

    # def match_pairs_redesign3(self):
    def match_pairs(self):
        # TODO: do we need to consider multi-matching in this method? e.g., men who are still underpartnered after calling this method once?

        ppl = self.sim.people
        # Find people eligible for a relationship
        active = self.over_debut
        underpartnered = self.partners < self.concurrency
        f_eligible = active & ppl.female & underpartnered
        m_eligible_uids = (active & ppl.male & underpartnered).uids
        f_looking_uids = self.pars.p_pair_form.filter(f_eligible.uids)  # ss.uids of women looking for partners

        paired_females = []
        paired_males = []
        i = 0
        male_bin_widths = [1, 3, 5]
        while len(f_looking_uids) > 0 and i < len(male_bin_widths):
            # on each successive pass, any females left looking for matches will re-select a target age for use
            # with a wider age binning of remaining males (representing females getting less choosy)
            loc, scale = self.get_age_risk_pars(f_looking_uids, self.pars.age_diff_pars)
            self.pars.age_diffs.set(loc=loc, scale=scale)
            target_age_gaps = self.pars.age_diffs.rvs(f_looking_uids)  # Sample the age differences

            male_bin_width = male_bin_widths[i]
            # divvy up men into bins
            male_bins = self.bin_agents_by_age(m_eligible_uids, bin_width=male_bin_width)  # returns a {age: df} dict
            # create female queues to male bins
            female_uid_bins = self.bin_female_uids_by_target_age(female_uids=f_looking_uids,
                                                                 target_age_gaps=target_age_gaps,
                                                                 male_age_bins=male_bins.keys(),
                                                                 bin_width=male_bin_width)
            # assign men to women by male queue. For each age queue, assign pairings until men or women run out.
            for male_age, male_uids in male_bins.items():
                female_uids = female_uid_bins[male_age]
                limit = min(len(male_uids), len(female_uids))
                # arbitrarily assign; no sorting has been done
                selected_male_uids = male_uids[:limit]
                selected_female_uids = female_uids[:limit]
                if len(selected_male_uids) != len(selected_female_uids):
                    raise Exception(
                        f"Uh oh, mismatched m/f uid lengths male age: {male_age} M: {len(selected_male_uids)} F: {len(selected_female_uids)}")
                paired_males.extend(selected_male_uids)
                paired_females.extend(selected_female_uids)
            i += 1
            # now update the males and females still looking for the next round of pairing.
            f_looking_uids = f_looking_uids - selected_female_uids
            m_eligible_uids = m_eligible_uids - selected_male_uids

        m_selected = uids(paired_males)
        f_selected = uids(paired_females)

        return m_selected, f_selected

    def match_pairs_redesign2(self):
    # def match_pairs(self):

        ppl = self.sim.people
        # TODO: do we need to consider multi-matching in this method? e.g., men who are still underpartnered after calling this method once?
        # Find people eligible for a relationship
        active = self.over_debut
        underpartnered = self.partners < self.concurrency
        f_eligible = active & ppl.female & underpartnered
        m_eligible = active & ppl.male & underpartnered
        f_looking = self.pars.p_pair_form.filter(f_eligible.uids)  # ss.uids of women looking for partners
        f_eligible_and_looking = f_eligible & f_looking

        # Get mean age differences and desired ages
        # desired_ages = ppl.age[f_looking] + age_gaps    # Desired ages of the male partners
        # ppl.desired_ages[f_looking] = ppl.age[f_looking] + age_gaps    # Desired ages of the male partners

        male_bin_widths = [1, 3, 5]  # TODO: fix code below to remove paired ppl for following rounds, 3, 5]
        i = 0
        paired_females = []
        paired_males = []
        while len(f_looking) > 0 and i < len(male_bin_widths):
            loc, scale = self.get_age_risk_pars(f_looking, self.pars.age_diff_pars)
            self.pars.age_diffs.set(loc=loc, scale=scale)
            age_gaps = self.pars.age_diffs.rvs(f_eligible_and_looking)  # Sample the age differences

            male_bin_width = male_bin_widths[i]
            # divvy up men into bins
            male_bins = self.bin_agents_by_age(m_eligible, bin_width=male_bin_width)  # returns a {age: df} dict
            # create female queues to male bins
            female_uid_bins = self.bin_female_uids_by_target_age(female_uids=f_eligible_and_looking.uids, target_age_gaps=age_gaps,
                                                             male_age_bins=male_bins.keys(), bin_width=male_bin_width)
            # assign men to women by male queue. For each queue, assign until men or women run out.
            for male_age, men in male_bins.items():
                female_uids = female_uid_bins[male_age]
                limit = min(len(men.uids), len(female_uids))
                # arbitrarily assign; no sorting has been done
                selected_male_uids = men.uids[:limit]
                selected_female_uids = female_uids[:limit]
                if len(selected_male_uids) != len(selected_female_uids):
                    raise Exception(f"Uh oh, mismatched m/f uid lengths male age: {male_age} M: {len(selected_male_uids)} F: {len(selected_female_uids)}")
                paired_males.extend(selected_male_uids)
                paired_females.extend(selected_female_uids)
            i += 1
            # now update the males and females still looking for the next round of pairing.
            # remove newly-ineligible males
            # remove newly-ineligible females
            f_eligible_and_looking[selected_female_uids] = False
            f_looking = f_looking - selected_female_uids
            m_eligible[selected_male_uids] = False


        m_selected = uids(paired_males)
        f_selected = uids(paired_females)  # uids(f_selected)

        p1 = m_selected
        p2 = f_selected
        # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        return p1, p2


    def match_pairs_redesign1(self):
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

        # >>>
        # Attempt 2 sketch
        """
        - divvy up men into 1 year age bins, centered on 0.5 ages
        - queue up the women into queues for the male bins
        - for bin B in bins:
        -   distribute men to women until either men or women run out for bin B
        - if men and women remain, repeat from A with age_bin_width += 1
        - Do this for age_bin_width <= maximum SD from input pars.age_diffs
        - Return with selected matched men and women.
        """
        # <<<

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # Attempt 1: proof of concept works much better, but horrendous algorithm
        from starsim.arrays import uids
        m_selected = []
        f_selected = []
        m_ages_list = list(m_ages[ind_m])
        for i in range(len(desired_ages)):
            done = False
            delta = 0
            while not done:
                index = bisect_left(m_ages_list, desired_ages[i]) + delta
                try:
                    if m_eligible.uids[ind_m][index] not in m_selected:
                        done = True
                except IndexError as e:
                    done = True
                delta += 1
            try:
                selected = m_eligible.uids[ind_m][index]
                m_selected.append(selected)
                f_selected.append(f_looking[i])
                # m_ages_list.pop(selected)
                # print(f"Looking for age: {desired_ages[i]} found age: {ppl.age[selected]} len(m_ages_list): {len(m_ages_list)}")
            except IndexError as e:
                pass  # female could not find a male of target age
        # if any men are competed over, drop duplicates
        m_selected = uids(m_selected)
        f_selected = uids(f_selected)  # uids(f_selected)

        p1 = m_selected
        p2 = f_selected
        # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        if len(ind_m) == 0 or len(ind_f) == 0:
            raise NoPartnersFound()

        return p1, p2

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


class SWNetwork(BaseNetwork):
    """Standalone sex-work contact network (FSW–client edges only).

    Tracks female sex workers (``fsw``) and clients of sex workers (``client``)
    with one-timestep partnerships weighted by ``sw_intensity``. Can be used on
    its own or alongside :class:`MFNetwork`. For backward-compatible MF + SW
    bundling on a single network, see :class:`StructuredSexual`.

    Args:
        pars (dict): Parameter overrides (see :class:`SWPars` for defaults).
        condom_data: Optional condom-use data (DataFrame, dict, or scalar).
        name (str): Network name (default: auto-assigned).
        **kwargs: Additional parameter overrides forwarded to ``update_pars``.
    """

    def __init__(self, pars=None, condom_data=None, name=None, **kwargs):
        super().__init__(name=name)
        self.define_pars(**SWPars())
        self.update_pars(pars, **kwargs)
        self.edge_types = {'sw': 0}
        self.define_states(*_sw_states())
        return

    def set_network_states(self, upper_age=None):
        super().set_network_states(upper_age=upper_age)  # set_debut
        self.set_sex_work(upper_age=upper_age)
        return

    def set_sex_work(self, upper_age=None):
        """Draw lifetime SW fate and per-agent age window for new agents.

        Falls back to ``debut`` for ``age_*_start`` and to ``inf`` for ``dur_*``
        when the corresponding distribution is ``None``, so the ``fsw``/``client``
        properties always see valid window bounds. ``debut`` must be set before
        this method runs.
        """
        f_uids, m_uids = self._get_uids(upper_age=upper_age)
        new_fsw = self.pars.fsw_shares.rvs(f_uids)
        new_client = self.pars.client_shares.rvs(m_uids)
        self.ever_fsw[f_uids] = new_fsw
        self.ever_client[m_uids] = new_client

        fsw_uids = f_uids[new_fsw]
        if len(fsw_uids):
            if self.pars.age_sw_start is not None:
                drawn = self.pars.age_sw_start.rvs(fsw_uids)
                # Clamp entry age to debut — can't be SW before sexual debut
                self.age_sw_start[fsw_uids] = np.maximum(drawn, self.debut[fsw_uids])
            else:
                self.age_sw_start[fsw_uids] = self.debut[fsw_uids]
            if self.pars.dur_sw is not None:
                self.dur_sw[fsw_uids] = self.pars.dur_sw.rvs(fsw_uids)
            else:
                self.dur_sw[fsw_uids] = np.inf

        client_uids = m_uids[new_client]
        if len(client_uids):
            if self.pars.age_client_start is not None:
                drawn = self.pars.age_client_start.rvs(client_uids)
                self.age_client_start[client_uids] = np.maximum(drawn, self.debut[client_uids])
            else:
                self.age_client_start[client_uids] = self.debut[client_uids]
            if self.pars.dur_client is not None:
                self.dur_client[client_uids] = self.pars.dur_client.rvs(client_uids)
            else:
                self.dur_client[client_uids] = np.inf
        return

    @property
    def age_sw_stop(self):
        return self.age_sw_start + self.dur_sw

    @property
    def age_client_stop(self):
        return self.age_client_start + self.dur_client

    @property
    def fsw(self):
        """Currently a female sex worker."""
        age = self.sim.people.age
        return (age >= self.age_sw_start) & (age < self.age_sw_stop) & self.ever_fsw

    @property
    def client(self):
        """Currently a client of sex workers."""
        age = self.sim.people.age
        return (age >= self.age_client_start) & (age < self.age_client_stop) & self.ever_client

    def match_pairs(self):
        """ Match sex workers to clients """
        active = self.over_debut
        active_fsw = active & self.fsw
        active_clients = active & self.client
        self.sw_intensity[active_fsw.uids] = self.pars.sw_intensity.rvs(active_fsw.uids)

        self.pars.sw_seeking_dist.pars.p = self.pars.sw_seeking_rate.to_prob()
        m_looking = self.pars.sw_seeking_dist.filter(active_clients.uids)

        if len(m_looking) == 0 or len(active_fsw.uids) == 0:
            raise NoPartnersFound()

        # Repeat-sample FSW weighted by intensity to assign each client a partner
        if len(m_looking) > len(active_fsw.uids):
            n_repeats = (self.sw_intensity[active_fsw]*10).astype(int)+1
            fsw_repeats = np.repeat(active_fsw.uids, n_repeats)
            if len(fsw_repeats) < len(m_looking):
                fsw_repeats = np.repeat(fsw_repeats, 10)

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

    def add_pairs(self):
        """ Match and add sex worker partnerships for this timestep.

        Each partnership has duration = 1 timestep. Updates ``lifetime_sw_partners``.
        Uses an explicit ``SWNetwork.match_pairs`` reference so subclasses that
        also inherit ``MFNetwork.match_pairs`` (e.g. :class:`StructuredSexual`)
        still pick the SW matcher here.
        """
        ppl = self.sim.people

        try:
            p1, p2 = SWNetwork.match_pairs(self)
        except NoPartnersFound:
            return

        match_count = len(p1)
        beta = np.ones(match_count, dtype=ss_float)
        condoms = np.zeros(match_count, dtype=ss_float)
        acts = (self.pars.acts.rvs(p2)).astype(int)
        dur = np.full(match_count, fill_value=1)
        age_p1 = ppl.age[p1]
        age_p2 = ppl.age[p2]
        edge_types = np.full(match_count, dtype=ss_float, fill_value=self.edge_types['sw'])

        self.append(p1=p1, p2=p2, beta=beta, condoms=condoms, dur=dur, acts=acts, age_p1=age_p1, age_p2=age_p2, edge_type=edge_types)

        p1_edges, p1_counts = np.unique(p1, return_counts=True)
        p2_edges, p2_counts = np.unique(p2, return_counts=True)

        self.lifetime_sw_partners[p1_edges] += p1_counts
        self.lifetime_sw_partners[p2_edges] += p2_counts

        return


class StructuredSexual(MFNetwork):
    """Heterosexual + sex-work network on a single edge list (backward compatible).

    Bundles :class:`MFNetwork` (risk groups, concurrency, age preferences) with
    sex-work edges. Maintains the historical API so existing consumers can read
    ``fsw``/``client``/``sw_intensity`` directly via ``sim.networks.structuredsexual``.

    For modular use — for example, MF without SW, or SW alone — use
    :class:`MFNetwork` and :class:`SWNetwork` directly.

    Args:
        pars (dict): Parameter overrides (see :class:`NetworkPars` for defaults).
        condom_data: Optional condom-use data (DataFrame or dict).
        name (str): Network name (default: auto-assigned).
        **kwargs: Additional parameter overrides forwarded to ``update_pars``.
    """

    def __init__(self, pars=None, condom_data=None, name=None, **kwargs):
        # MF layer (defaults only — no user pars/kwargs yet)
        super().__init__(name=name)
        # SW layer
        self.define_pars(**SWPars())
        self.define_states(*_sw_states())
        self.edge_types['sw'] = max(self.edge_types.values()) + 1
        # Apply user pars/kwargs against the full pars dict
        self.update_pars(pars, **kwargs)
        return

    # SW behavior aliased from SWNetwork (single source of truth).
    # Note: ``match_pairs`` is intentionally not aliased — MFNetwork's version
    # is used for MF matching; SWNetwork.add_pairs calls SWNetwork.match_pairs explicitly.
    set_sex_work = SWNetwork.set_sex_work
    fsw = SWNetwork.fsw
    client = SWNetwork.client
    age_sw_stop = SWNetwork.age_sw_stop
    age_client_stop = SWNetwork.age_client_stop

    def set_network_states(self, upper_age=None):
        super().set_network_states(upper_age=upper_age)
        self.set_sex_work(upper_age=upper_age)
        return

    def add_pairs(self):
        MFNetwork.add_pairs(self)
        SWNetwork.add_pairs(self)
        return

    # end_pairs and set_condom_use inherited from BaseNetwork. _decrement_partners
    # (in MFNetwork) skips SW edges via the ``'sw' in self.edge_types`` check;
    # set_condom_use's per-key dispatch handles ``(rgm, rgf)`` and ``('fsw','client')``
    # entries automatically.


# %% Auxiliary networks — PriorPartners (recall buffer) and MSM variants

class PriorPartners(ss.DynamicNetwork):
    """Lightweight network that tracks prior sexual partners for partner notification.

    Stores edges representing ended relationships and increments their duration
    each timestep. Edges older than ``dur_recall`` are removed. Used by partner
    notification interventions to trace recent contacts.

    Args:
        pars (dict): Parameter overrides; key parameter is ``dur_recall``
            (default 1 year).
        name (str): Network name (default ``'priorpartners'``).
        **kwargs: Additional parameter overrides.
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


class AgeMatchedMSM(MFNetwork):
    """Men-who-have-sex-with-men network using exact age-sorted matching.

    Extends :class:`StructuredSexual` for MSM partnerships. Eligible males
    are sorted by age and paired sequentially so that partners have similar
    ages. The ``msm_share`` parameter controls what fraction of males
    participate.

    Args:
        pars (dict): Parameter overrides; key parameter is ``msm_share``
            (default ``ss.bernoulli(p=0.015)``).
        **kwargs: Additional parameter overrides.
    """

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


class AgeApproxMSM(MFNetwork):
    """Men-who-have-sex-with-men network using approximate age-preference matching.

    Extends :class:`StructuredSexual` for MSM partnerships. Unlike
    :class:`AgeMatchedMSM`, this variant splits eligible males into two
    arbitrary groups and matches them using the standard age-difference
    preference distributions rather than exact age sorting.

    Args:
        **kwargs: Parameter overrides forwarded to :class:`StructuredSexual`.
    """

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