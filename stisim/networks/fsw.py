"""Sex-work (FSW–client) sexual network and its parameters.

The class is named ``SWNetwork`` (unchanged); the file is named ``fsw.py``
for discoverability.
"""
import numpy as np
import starsim as ss
from .base import BaseNetwork, NoPartnersFound, ss_float

__all__ = ['SWPars', 'SWNetwork']


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
        ss.BoolArr('paused'),                  # Reversible exclusion from FSW status (e.g. during pregnancy); leaves ever_fsw untouched
    ]


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
        """Currently a female sex worker.

        Derived from the lifetime ``ever_fsw`` flag, the per-agent
        ``age_sw_start``/``age_sw_stop`` window, AND a reversible
        ``paused`` flag (default ``False`` for everyone). Interventions
        that want to temporarily exclude an agent from FSW status
        (e.g. ``sti.PregnancyRiskReduction``) toggle ``paused`` instead
        of mutating any of the lifetime fields.
        """
        age = self.sim.people.age
        return (age >= self.age_sw_start) & (age < self.age_sw_stop) & self.ever_fsw & ~self.paused

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
        # TODO: consider rejecting a (p1, p2) pairing if that pair already has an
        # active edge in this or any other known network, to enforce the
        # no-concurrent-duplicate-edge invariant that partner-uniqueness
        # reporting (e.g. PartnershipFormationAnalyzer) assumes.
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

        ti_formed = np.full(match_count, self.ti, dtype=int)
        self.append(p1=p1, p2=p2, beta=beta, condoms=condoms, dur=dur, acts=acts, age_p1=age_p1, age_p2=age_p2, edge_type=edge_types, ti_formed=ti_formed)

        p1_edges, p1_counts = np.unique(p1, return_counts=True)
        p2_edges, p2_counts = np.unique(p2, return_counts=True)

        self.lifetime_sw_partners[p1_edges] += p1_counts
        self.lifetime_sw_partners[p2_edges] += p2_counts

        return
