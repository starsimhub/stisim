"""Layered networks — composite networks that combine multiple network modules.

A "layered network" is a single class that bundles two or more network
behaviours so a sim can model their combined dynamics without juggling
multiple network objects. The canonical example is :class:`StructuredSexual`,
which bundles heterosexual partnerships (:class:`MFNetwork`) and sex-work
partnerships (:class:`SWNetwork`) on a single edge list.

:class:`PriorPartners` is a supplementary network you add *alongside* an MF
network when ``recall_prior=True`` — it stores recently-dissolved edges so
partner-notification interventions can trace them.
"""
import numpy as np
import starsim as ss
from .base import BaseNetwork, NoPartnersFound, ss_float
from .mf import MFNetwork, _mf_states
from .fsw import SWNetwork, SWPars, _sw_states

__all__ = ['StructuredSexual', 'PriorPartners']


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
        # Multiplier on FSW non-sex-work (MF) concurrency: values <1 mean
        # active sex workers have fewer non-sex-work partners. Default 1.0 =
        # no effect. Lives here rather than on MFNetwork because StructuredSexual
        # (MF + SW combined) is the only consumer.
        self.define_pars(fsw_mf_conc_mult=1.0)
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
        # Scale FSW agents' MF (non-sex-work) concurrency, now that self.fsw is
        # populated. This affects only sex workers' non-commercial partnerships;
        # their sex-work edges come from SWNetwork.add_pairs and are untouched.
        # No-op at default 1.0.
        mult = self.pars.fsw_mf_conc_mult
        if mult != 1.0:
            fsw_uids = self.fsw.uids
            if len(fsw_uids) > 0:
                scaled = np.maximum(1, np.round(self.concurrency[fsw_uids] * mult)).astype(int)
                self.concurrency[fsw_uids] = scaled
        return

    def add_pairs(self):
        MFNetwork.add_pairs(self)
        SWNetwork.add_pairs(self)
        return

    # end_pairs and set_condom_use inherited from BaseNetwork. _decrement_partners
    # (in MFNetwork) skips SW edges via the ``'sw' in self.edge_types`` check;
    # set_condom_use's per-key dispatch handles ``(rgm, rgf)`` and ``('fsw','client')``
    # entries automatically.


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
