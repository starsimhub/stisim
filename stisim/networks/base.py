"""Base classes shared by all STIsim sexual networks."""

import starsim as ss
import sciris as sc
import numpy as np
import pandas as pd

ss_float = ss.dtypes.float
ss_int = ss.dtypes.int

__all__ = ['NoPartnersFound', 'BasePars', 'NetworkPars', 'BaseNetwork', 'ss_float', 'ss_int']


class NoPartnersFound(Exception):
    # Raise this exception if the matching algorithm wasn't able to match any partners
    pass


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


class NetworkPars(ss.Pars):
    """Combined defaults — union of base + MF + SW pars.

    Used by :class:`sti.Sim` to filter network kwargs into the right bucket.

    Args:
        **kwargs: Any parameter name-value pairs to override defaults.
    """
    def __init__(self, **kwargs):
        super().__init__()
        from .mf import MFPars
        from .fsw import SWPars
        for src in (BasePars(), MFPars(), SWPars()):
            self.update(dict(src), create=True)
        self.update(kwargs)
        return


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
        self.meta.formation_ti = ss_int  # timestep at which the edge was formed
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
            from .layered_networks import PriorPartners
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
        # TODO: consider rejecting a (p1, p2) pairing if that pair already has an
        # active edge in this or any other known network, to enforce the
        # no-concurrent-duplicate-edge invariant that partner-uniqueness
        # reporting (e.g. PartnershipFormationAnalyzer) assumes.
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
