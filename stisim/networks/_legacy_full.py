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
from .fsw import SWPars, SWNetwork, _sw_states

# Specify all externally visible functions this file defines; see also more definitions below
__all__ = ['StructuredSexual',
           'PriorPartners', 'AgeMatchedMSM', 'AgeApproxMSM']


# %% Network classes — concrete subclasses of BaseNetwork

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