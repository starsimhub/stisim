"""Men-who-have-sex-with-men (MSM) sexual networks.

Three MSM variants live here:

- :class:`AgeMatchedMSM` and :class:`AgeApproxMSM` extend :class:`MFNetwork`
  and use the standard discrete pair-matching paradigm (one match per step).
- :class:`MSMScaleFreeNetwork` extends :class:`BaseNetwork` directly and uses
  a continuous-time preferential-attachment rate kernel (Whittles 2019 S2).
  It is not branching-stable under starsim CRN — see its docstring.
"""
import numpy as np
import starsim as ss
from .base import NoPartnersFound, BaseNetwork
from .mf import MFNetwork

__all__ = ['AgeMatchedMSM', 'AgeApproxMSM', 'MSMScaleFreeNetwork']


class AgeMatchedMSM(MFNetwork):
    """Men-who-have-sex-with-men network using exact age-sorted matching.

    Extends :class:`MFNetwork` for MSM partnerships. Eligible males
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
        super().set_network_states(upper_age=upper_age)  # debut, risk groups, concurrency
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
        m_eligible = active & ppl.male & underpartnered & self.participant
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


class MSMScaleFreeNetwork(BaseNetwork):
    """Preferential-attachment MSM sexual network with continuous-time
    formation and Markovian deletion (Whittles-2019 S2 kernel).

    Contributed by Stephen Attwood (Postdoctoral Researcher, Department of
    Biology, University of Oxford), June 2026. Originally prototyped as
    ``BespokeNet`` in the ``starsim_x`` repository.

    Eligible agents are post-debut males. Edges form via a rate kernel
    proportional to ``g(λ_i) g(λ_j) × pair_weight``, where the default
    pair weight is a degree-based multiplier (rich-get-richer). Edges are
    deleted uniformly at random subject to a hard duration cap.

    Subclasses can plug in age- or risk-weighted formation by overriding
    ``_mix_node_arrays`` and ``_mix_weights_row``, or change the eligible
    pool (e.g. heterosexual) by overriding ``_get_pool``.

    **CRN caveat.** Unlike most stisim networks, this class is not
    branching-stable under starsim's Common Random Numbers framework.
    Edge formation uses global weighted-categorical draws over a pair
    catalog, which couples every draw to the entire current pool — there
    is no per-agent slot the way ``ss.bernoulli.filter(uids)`` provides.
    Reproducibility is preserved per ``(rand_seed, ti)``, but adding an
    intervention that perturbs the population will reshuffle edges from
    that point on. This matches the precedent set by the
    ``desired_age_bucket`` matcher in :mod:`stisim.networks.matchers`.

    Args:
        target_mean_degree (float): target mean number of concurrent
            partners per agent in the pool (default 2.0).
        target_mean_dur (ss.dur or int): target mean edge duration. Accepts
            an ``ss.dur`` (converted to integer steps at ``init_pre`` via
            ``self.t.dt``) or a raw integer step count. Default ``ss.years(2)``.
        max_edge_dur (ss.dur or int): hard cap on edge persistence. Default
            ``ss.years(10)``.
        phi (float): Whittles-2019 turnover parameter. Default 1.0. Sets
            the q0 initial-network density via ``q0 = 1/(1+phi)``.
        name (str): network name (default: auto-assigned).
        **kwargs: forwarded to ``update_pars``.
    """

    def __init__(self, pars=None, name=None, **kwargs):
        super().__init__(name=name)
        self.define_pars(
            msm_share=ss.bernoulli(p=0.015),  # fraction of males in the MSM pool
            target_mean_degree=2.0,
            target_mean_dur=ss.years(2),
            max_edge_dur=ss.years(10),
            phi=1.0,
        )
        self.update_pars(pars, **kwargs)
        self._target_mean_dur_steps = None
        self._max_edge_dur_steps = None
        self._kernel_pool_uids = None
        self._kernel_pairs_i = None
        self._kernel_pairs_j = None
        self._kernel_A_w = None
        self._kernel_sel_w = None
        self._kernel_sel_cdf = None  # cached cdf for rng.choice fast-path
        self._kernel_hat_Ra = None
        self._kernel_q0 = None
        # BespokeNet defaults — promote to user pars only when needed.
        self._k1 = 1.0  # consumed by _g_of_lambda
        self._rebuild_every = 10  # consumed by step()
        # Subclass extension knobs — read by overridden _mix_weights_row.
        self._age_a = 0.0
        self._partnering_assortivity = 1.0
        return

    def init_pre(self, sim):
        """Convert ``ss.dur`` parameters to integer step counts.

        Uses ``ss.dur / ss.dur`` arithmetic (sciris-native unit conversion).
        """
        super().init_pre(sim)
        self._target_mean_dur_steps = int(round(self.pars.target_mean_dur / self.t.dt))
        self._max_edge_dur_steps = int(round(self.pars.max_edge_dur / self.t.dt))
        if self._max_edge_dur_steps < self._target_mean_dur_steps:
            raise ValueError(
                f'max_edge_dur ({self._max_edge_dur_steps} steps) must be >= '
                f'target_mean_dur ({self._target_mean_dur_steps} steps)'
            )
        return

    def set_network_states(self, upper_age=None):
        super().set_network_states(upper_age=upper_age)  # set debut
        self.set_msm(upper_age=upper_age)
        return

    def set_msm(self, upper_age=None):
        _, m_uids = self._get_uids(upper_age=upper_age)
        self.participant[m_uids] = self.pars.msm_share.rvs(m_uids)
        return

    def _get_rng(self):
        """Step-local numpy Generator. Not CRN — see class docstring."""
        seed = int(self.sim.pars.rand_seed) + 2003 + int(self.ti)
        return np.random.default_rng(seed)

    # ---- Subclass extension points -----------------------------------------
    def _get_pool(self):
        """Eligible agents the kernel operates on. Default: post-debut MSM males."""
        return self.over_debut & self.sim.people.male & self.participant

    def _mix_node_arrays(self):
        """Per-node arrays consumed by ``_mix_weights_row``.

        Base implementation provides ``nodes`` (sorted uids) and
        ``log1p_deg`` (1 + log1p of current degree). Subclasses extend
        the dict with age/risk vectors as needed.
        """
        pool_uids = self._get_pool().uids
        n = pool_uids.size
        if n == 0:
            return {'nodes': pool_uids, 'log1p_deg': np.empty(0, dtype=float)}
        endpoints = np.concatenate([self.edges.p1, self.edges.p2])
        if endpoints.size:
            idx = np.searchsorted(pool_uids, endpoints)
            in_bounds = idx < n
            hit = np.zeros_like(endpoints, dtype=bool)
            hit[in_bounds] = pool_uids[idx[in_bounds]] == endpoints[in_bounds]
            deg_vec = np.bincount(idx[hit], minlength=n).astype(float)
        else:
            deg_vec = np.zeros(n, dtype=float)
        return {'nodes': pool_uids, 'log1p_deg': 1.0 + np.log1p(deg_vec)}

    def _mix_weights_row(self, i, mix_arrays):
        """Vectorised pair weights between node i and all j > i.

        Returns 1D array of length (n - 1 - i) of non-negative weights.
        Default = preferential-attachment degree multiplier.
        """
        log1p_deg = mix_arrays['log1p_deg']
        return log1p_deg[i] * log1p_deg[i+1:]

    # ---- S2 kernel ---------------------------------------------------------
    def _g_of_lambda(self, lam):
        """``g(λ) = λ / (λ + k1)`` — smooth monotone saturation."""
        lam = np.asarray(lam, dtype=float)
        return lam / (lam + self._k1)

    def _build_kernel(self):
        """Build the S2 formation kernel against the current pool.

        Populates ``_kernel_pairs_{i,j}``, ``_kernel_A_w``, ``_kernel_sel_w``,
        ``_kernel_sel_cdf`` (normalised cumsum for fast sampling),
        ``_kernel_hat_Ra``, and ``_kernel_q0``. Returns the pool size ``n``;
        0 on a degenerate pool.

        Deterministic given the pool — no RNG needed. Mirrors
        ``starsim_x/networks.py:3106 (_s2_build_A_and_q0)``.
        """
        mix_arrays = self._mix_node_arrays()
        pool_uids = mix_arrays['nodes']
        n = pool_uids.size
        phi = float(self.pars.phi)
        if n < 2:
            self._kernel_pool_uids = pool_uids
            self._kernel_pairs_i = np.empty(0, dtype=int)
            self._kernel_pairs_j = np.empty(0, dtype=int)
            self._kernel_A_w = np.empty(0, dtype=float)
            self._kernel_sel_w = np.empty(0, dtype=float)
            self._kernel_sel_cdf = np.empty(0, dtype=float)
            self._kernel_hat_Ra = 0.0
            self._kernel_q0 = 1.0 / (1.0 + phi)
            return 0

        lambdas = np.ones(n, dtype=float)
        g_raw = self._g_of_lambda(lambdas)
        g_sum = float(np.sum(g_raw))
        g_den = max(float(np.sqrt(g_sum / max(phi, 1.0))), 1e-12)
        g = g_raw / g_den

        # Row loop (upper triangle). S2 step 5 (Whittles 2019):
        #   A_ij = -ln[(1 + 1/phi) - g(λ_i) g(λ_j) (1 + 1/phi)]
        idx_i_list, idx_j_list, A_w_list, sel_w_list = [], [], [], []
        g_fac = 1.0 + 1.0 / phi
        sqrt_g_fac = np.sqrt(g_fac)
        for i in range(n - 1):
            gi = float(g[i])
            gj_all = g[i+1:]
            gl_i = gi * sqrt_g_fac
            gl_all = gj_all * sqrt_g_fac
            inner = np.clip(g_fac - gl_i * gl_all, 1e-15, None)
            base = np.maximum(-np.log(inner), 0.0)

            # Selection weight = base × normalised mix_weight.
            w = self._mix_weights_row(i, mix_arrays)
            wm = float(w.mean()) if w.size else 1.0
            w_norm = w / max(wm, 1e-12)
            sel_base = base * w_norm

            if np.any(base > 0):
                thresh = 1e-4 * float(np.median(base[base > 0]))
            else:
                thresh = 0.0
            keep = base > thresh
            if np.any(keep):
                jj = np.nonzero(keep)[0] + (i + 1)
                idx_i_list.append(np.full(jj.size, i, dtype=int))
                idx_j_list.append(jj.astype(int))
                A_w_list.append(base[keep].astype(float))
                sel_w_list.append(sel_base[keep].astype(float))

        if idx_i_list:
            pairs_i = np.concatenate(idx_i_list)
            pairs_j = np.concatenate(idx_j_list)
            A_w = np.concatenate(A_w_list)
            sel_w = np.concatenate(sel_w_list)
        else:
            pairs_i = np.empty(0, dtype=int)
            pairs_j = np.empty(0, dtype=int)
            A_w = np.empty(0, dtype=float)
            sel_w = np.empty(0, dtype=float)

        A_w[A_w < 0.0] = 0.0
        sel_w[sel_w < 0.0] = 0.0

        # Fallback for uniform-λ case: when the Whittles formula collapses
        # (no biological heterogeneity), build a dense kernel from mix_weights
        # alone with a constant rate baseline. Lets the preferential-attachment
        # hook drive the network structure without an external λ vector.
        if pairs_i.size == 0:
            ii, jj, ww = [], [], []
            for i in range(n - 1):
                row = self._mix_weights_row(i, mix_arrays)
                if row.size and row.max() > 0:
                    ii.append(np.full(row.size, i, dtype=int))
                    jj.append(np.arange(i+1, n, dtype=int))
                    ww.append(row.astype(float))
            if ii:
                pairs_i = np.concatenate(ii)
                pairs_j = np.concatenate(jj)
                sel_w = np.concatenate(ww)
                A_w = np.ones_like(sel_w)

        if sel_w.size > 0 and sel_w.sum() > 0.0:
            cs = np.cumsum(sel_w)
            sel_cdf = cs / cs[-1]
            sel_cdf[-1] = 1.0
        else:
            sel_cdf = np.empty(0, dtype=float)

        # Global add-event clock from target stock + duration.
        # Note: Whittles' q0 = 1/(1+phi) is the stationary density over the
        # *coupled* (non-zero A_ij) pair set. With heterogeneous λ that's
        # naturally sparse; with the uniform-λ fallback the catalog spans
        # all upper-triangle pairs, so we derive q0 from the target stock
        # directly to avoid over-densification at init.
        E_star = 0.5 * float(n) * float(self.pars.target_mean_degree)
        m_pairs = pairs_i.size
        if m_pairs > 0:
            q0 = min(1.0, E_star / float(m_pairs))
        else:
            q0 = 1.0 / (1.0 + phi)
        D = max(float(self._target_mean_dur_steps), 1e-12)
        mu_turnover = E_star / D
        hat_Ra = float(mu_turnover / max(1.0 - q0, 1e-12))

        self._kernel_pool_uids = pool_uids
        self._kernel_pairs_i = pairs_i
        self._kernel_pairs_j = pairs_j
        self._kernel_A_w = A_w
        self._kernel_sel_w = sel_w
        self._kernel_sel_cdf = sel_cdf
        self._kernel_hat_Ra = hat_Ra
        self._kernel_q0 = q0
        return n

    def _sample_pair_index(self, rng):
        """Draw one pair index from the cached selection CDF.

        Returns -1 if the kernel holds no mass (caller decides what to do).
        """
        if self._kernel_sel_cdf is None or self._kernel_sel_cdf.size == 0:
            return -1
        u = float(rng.random())
        k = int(np.searchsorted(self._kernel_sel_cdf, u, side='left'))
        if k >= self._kernel_sel_cdf.size:
            k = self._kernel_sel_cdf.size - 1
        return k

    # ---- Initial network ---------------------------------------------------
    def _sample_P0(self, rng):
        """Sample initial-network edges via Bernoulli(q0) per S2 Step-3."""
        m = self._kernel_pairs_i.size
        empty = ss.uids()
        if m == 0:
            return empty, empty
        draws = rng.random(m) < self._kernel_q0
        if not draws.any():
            return empty, empty
        pool = self._kernel_pool_uids
        return ss.uids(pool[self._kernel_pairs_i[draws]]), ss.uids(pool[self._kernel_pairs_j[draws]])

    def _append_new_edges(self, p1, p2):
        """Append edges with default attributes + ``max_edge_dur`` lifetime.

        ``BaseNetwork.end_pairs`` decrements ``dur`` each step, so setting
        ``dur = max_edge_dur_steps`` at creation gives the hard cap for free.
        """
        n = len(p1)
        if n == 0:
            return
        ages = self.sim.people.age
        self.append(
            p1=p1, p2=p2,
            beta=np.ones(n, dtype=float),
            dur=np.full(n, float(self._max_edge_dur_steps), dtype=float),
            acts=np.ones(n, dtype=int),
            condoms=np.zeros(n, dtype=float),
            age_p1=ages[p1],
            age_p2=ages[p2],
            edge_type=np.zeros(n, dtype=float),
        )
        return

    def init_post(self):
        """Build kernel once + sample the q0 initial network."""
        super().init_post()
        n = self._build_kernel()
        if n < 2:
            return
        rng = np.random.default_rng(int(self.sim.pars.rand_seed) + 1001)
        p1, p2 = self._sample_P0(rng)
        if len(p1):
            self._append_new_edges(p1, p2)
        return

    # ---- Step dynamics -----------------------------------------------------
    def step(self):
        """One step of S2 dynamics.

        Rebuilds the kernel every ``_rebuild_every`` steps, then defers to
        ``BaseNetwork.step()`` which runs ``end_pairs`` (dur decrement →
        ``max_edge_dur`` hard cap), ``set_network_states``, ``add_pairs``
        (our S2 add + Markov delete), and ``set_condom_use``.
        """
        if self.ti > 0 and (self.ti % self._rebuild_every == 0):
            self._build_kernel()
        super().step()
        return

    def add_pairs(self):
        """S2 add + Markov delete, called from ``BaseNetwork.step()``.

        Adds are Poisson-counted draws from the sel-weight CDF; deletes
        are Poisson-counted uniform-at-random removals from the current
        edge list. Both target ``mu_turnover = E* / D``, balancing the
        stationary edge stock.
        """
        if self._kernel_sel_cdf is None:
            return
        rng = self._get_rng()
        rate = float(self._kernel_hat_Ra) * (1.0 - float(self._kernel_q0))
        if rate <= 0.0:
            return

        # Add events.
        n_adds = int(rng.poisson(rate))
        if n_adds and self._kernel_sel_cdf.size > 0:
            pool = self._kernel_pool_uids
            alive = self.sim.people.alive
            new_p1, new_p2 = [], []
            for _ in range(n_adds):
                k = self._sample_pair_index(rng)
                if k < 0:
                    break
                u1 = int(pool[self._kernel_pairs_i[k]])
                u2 = int(pool[self._kernel_pairs_j[k]])
                # Pool may be stale between rebuilds; skip dead endpoints.
                if alive[u1] and alive[u2]:
                    new_p1.append(u1)
                    new_p2.append(u2)
            if new_p1:
                self._append_new_edges(ss.uids(new_p1), ss.uids(new_p2))

        # Markov deletes: uniform among current edges.
        n_cur = len(self.edges.p1)
        if n_cur > 0:
            n_dels = min(int(rng.poisson(rate)), n_cur)
            if n_dels > 0:
                keep = np.ones(n_cur, dtype=bool)
                idx = rng.choice(n_cur, size=n_dels, replace=False)
                keep[idx] = False
                for k in self.meta_keys():
                    self.edges[k] = self.edges[k][keep]
        return


class AgeApproxMSM(MFNetwork):
    """Men-who-have-sex-with-men network using approximate age-preference matching.

    Extends :class:`MFNetwork` for MSM partnerships. Unlike
    :class:`AgeMatchedMSM`, this variant splits eligible males into two
    groups and matches them using the standard age-difference preference
    distributions rather than exact age sorting.

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
        super().set_network_states(upper_age=upper_age)  # debut, risk groups, concurrency
        self.set_msm(upper_age=upper_age)
        return

    def set_msm(self, upper_age=None):
        _, m_uids = self._get_uids(upper_age=upper_age)
        self.participant[m_uids] = self.pars.msm_share.rvs(m_uids)
        return

    def match_pairs(self):
        """ Match pairs using age preferences """
        ppl = self.sim.people

        # Find males eligible for a relationship
        active = self.over_debut
        underpartnered = self.partners < self.concurrency
        m_eligible = active & ppl.male & underpartnered & self.participant
        m_looking = self.pars.p_pair_form.filter(m_eligible.uids)
        if len(m_looking) < 2:
            raise NoPartnersFound()

        # Split the males looking for partners into two groups, then pair
        # group1 with group2 by age preference (closest-age matching).
        group1 = m_looking[::2]
        group2 = m_looking[1::2]
        loc, scale = self.get_age_risk_pars(group1, self.pars.age_diff_pars)
        self.pars.age_diffs.set(loc=loc, scale=scale)
        age_gaps = self.pars.age_diffs.rvs(group1)
        desired_ages = ppl.age[group1] + age_gaps
        g2_ages = ppl.age[group2]
        p1 = group1[np.argsort(desired_ages)]
        p2 = group2[np.argsort(g2_ages)]
        maxlen = min(len(p1), len(p2))
        return p1[:maxlen], p2[:maxlen]
