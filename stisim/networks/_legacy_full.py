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

from .base import NoPartnersFound, BasePars, NetworkPars, BaseNetwork, ss_float, ss_int
from .mf import MFPars, MFNetwork, _mf_states

# Specify all externally visible functions this file defines; see also more definitions below
__all__ = ['SWPars',
           'StructuredSexual', 'SWNetwork',
           'PriorPartners', 'AgeMatchedMSM', 'AgeApproxMSM']


# %% State helpers — lists of agent-level states fed to ``define_states``

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


# %% Network classes — concrete subclasses of BaseNetwork

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