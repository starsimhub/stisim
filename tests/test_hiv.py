import numpy as np
import sciris as sc
import matplotlib.pyplot as plt
import starsim as ss
import stisim as sti


class TrackValues(ss.Analyzer):
    # Track outputs for viral load and CD4 counts
    # Assumes no births; for diagnostic/debugging purposes only
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.n = len(sim.people)

        self.hiv_rel_sus = np.empty((sim.t.npts, self.n), dtype=ss.dtypes.float)
        self.hiv_rel_trans = np.empty((sim.t.npts, self.n), dtype=ss.dtypes.float)

        self.syph_rel_sus = np.empty((sim.t.npts, self.n), dtype=ss.dtypes.float)
        self.syph_rel_trans = np.empty((sim.t.npts, self.n), dtype=ss.dtypes.float)

        self.cd4 = np.empty((sim.t.npts, self.n), dtype=ss.dtypes.float)

        self.care_seeking = np.empty((sim.t.npts, self.n), dtype=ss.dtypes.float)

    @property
    def has_hiv(self):
        return 'hiv' in self.sim.diseases

    @property
    def has_syph(self):
        return isinstance(self.sim.diseases.get('syph'), sti.Syphilis)

    def step(self):
        ti = self.ti

        if self.has_hiv:
            hiv = self.sim.diseases.hiv
            self.hiv_rel_sus[ti, :self.n] = hiv.rel_sus.raw
            self.hiv_rel_trans[ti, :self.n] = hiv.rel_trans.raw
            self.cd4[ti, :self.n] = hiv.cd4.raw
            self.care_seeking[ti, :self.n] = hiv.care_seeking.raw

        if self.has_syph:
            syph = self.sim.diseases.syphilis
            self.syph_rel_sus[ti, :self.n] = syph.rel_sus.raw
            self.syph_rel_trans[ti, :self.n] = syph.rel_trans.raw

        return

    def plot(self, agents: dict):
        """
        :param agents: Dictionary of events per agent {'agent_description':[('event_type', ti),...]}
        :return: Matplotlib figure
        """

        def plot_with_events(ax, x, y, agents, title):
            h = ax.plot(x, y)
            x_ev = []
            y_ev = []
            for i, events in enumerate(agents.values()):
                for event in events:
                    x_ev.append(self.sim.timevec[event[1]])
                    y_ev.append(y[event[1], i])
            ax.scatter(x_ev, y_ev, marker='*', color='yellow', edgecolor='red', s=100, linewidths=0.5, zorder=100)
            ax.set_title(title)
            return h

        if self.has_syph:
            fig, ax = plt.subplots(2, 4)
        else:
            fig, ax = plt.subplots(1, 2)

        ax = ax.ravel()

        h = plot_with_events(ax[0], self.sim.timevec, self.cd4, agents, 'CD4')
        # h = plot_with_events(ax[1], self.sim.timevec, self.hiv_rel_sus, agents, 'HIV rel_sus')
        h = plot_with_events(ax[1], self.sim.timevec, self.hiv_rel_trans, agents, 'HIV rel_trans')
        # h = plot_with_events(ax[3], self.sim.timevec, self.care_seeking, agents, 'HIV care seeking')

        if self.has_syph:
            h = plot_with_events(ax[4], self.sim.timevec, self.syph_rel_sus, agents, 'Syphilis rel_sus')
            h = plot_with_events(ax[5], self.sim.timevec, self.syph_rel_trans, agents, 'Syphilis rel_trans')

        # fig.legend(h, agents.keys(), loc='upper right', bbox_to_anchor=(1.1, 1))

        return fig


class PerformTest(ss.Intervention):

    def __init__(self, events=None):
        """
        :param events: List of (uid, 'event', ti) to apply events to an agent
        """
        super().__init__()
        self.define_pars(
            dt='month',
        )
        self.hiv_infections = sc.ddict(list)
        self.syphilis_infections = sc.ddict(list)
        self.art_start = sc.ddict(list)
        self.art_stop = sc.ddict(list)
        self.pregnant = sc.ddict(list)

        if events:
            for uid, event, ti in events:
                if event == 'hiv_infection':
                    self.hiv_infections[ti].append(uid)
                elif event == 'syphilis_infection':
                    self.syphilis_infections[ti].append(uid)
                elif event == 'art_start':
                    self.art_start[ti].append(uid)
                elif event == 'art_stop':
                    self.art_stop[ti].append(uid)
                elif event == 'pregnant':
                    self.pregnant[ti].append(uid)
                else:
                    raise Exception(f'Unknown event "{event}"')

    def initiate_ART(self, uids):
        if len(uids):
            self.sim.diseases.hiv.start_art(ss.uids(uids))

    def end_ART(self, uids):
        if len(uids):
            self.sim.diseases.hiv.stop_art(ss.uids(uids))

    def step(self):
        sim = self.sim
        ti = self.ti
        self.initiate_ART(self.art_start[ti])
        self.end_ART(self.art_stop[ti])
        if 'hiv' in sim.diseases:
            self.sim.diseases.hiv.set_prognoses(ss.uids(self.hiv_infections[ti]))
        if 'syphilis' in sim.diseases:
            self.sim.diseases.syphilis.set_prognoses(ss.uids(self.syphilis_infections[ti]))


def test_hiv():
    # AGENTS
    agents = sc.odict()
    agents['No infection'] = []
    agents['Infection without ART'] = [('hiv_infection', 5)]
    agents['Goes onto ART early (CD4 > 200) and stays on forever'] = [('hiv_infection', 4), ('art_start', 1 * 12)]
    agents['Goes onto ART late (CD4 < 200) and stays on forever'] = [('hiv_infection', 3), ('art_start', 10 * 12)]
    agents['Goes off ART with CD4 > 200'] = [('hiv_infection', 2), ('art_start', 5 * 12), ('art_stop', 12 * 12)]
    agents['Goes off ART with CD4 < 200'] = [('hiv_infection', 2), ('art_start', 12 * 12), ('art_stop', 13 * 12)]
    agents['pregnant'] = [('pregnant', 5), ('hiv_infection', 10)]

    events = []
    for i, x in enumerate(agents.values()):
        for y in x:
            events.append((i,) + y)

    pars = {}
    pars['n_agents'] = len(agents)
    pars['start'] = 2020
    pars['stop'] = 2040
    hiv = sti.HIV(init_prev=0, p_hiv_death=0, include_aids_deaths=False)
    pars['diseases'] = [hiv]
    pars['demographics'] = [ss.Pregnancy(fertility_rate=0), ss.Deaths(death_rate=0)]
    pars['interventions'] = PerformTest(events)
    output = TrackValues()
    pars['analyzers'] = output

    sim = sti.Sim(pars).run()
    sim.analyzers.trackvalues.plot(agents)
    return sim


def test_hiv_syph():
    # AGENTS
    agents = sc.odict()
    agents['No infection'] = []
    agents['HIV only'] = [('hiv_infection', 5)]
    agents['HIV before syphilis'] = [('syphilis_infection', 1), ('hiv_infection', 12)]
    agents['HIV after syphilis'] = [('hiv_infection', 1), ('syphilis_infection', 12)]

    events = []
    for i, x in enumerate(agents.values()):
        for y in x:
            events.append((i,) + y)

    pars = {}
    pars['n_agents'] = len(agents)
    pars['start'] = 2020
    pars['stop'] = 2040
    pars['dt'] = 1 / 12
    hiv = sti.HIV(init_prev=0, p_hiv_death=0, include_aids_deaths=False, beta={'structuredsexual': [0, 0], 'maternal': [0, 0]})
    syphilis = sti.SyphilisPlaceholder(prevalence=None)

    pars['diseases'] = [hiv, syphilis]
    pars['networks'] = [sti.StructuredSexual(), ss.MaternalNet()]
    pars['demographics'] = [ss.Pregnancy(fertility_rate=0), ss.Deaths(death_rate=0)]
    pars['interventions'] = PerformTest(events)
    output = TrackValues()
    pars['analyzers'] = output
    pars['connectors'] = sti.hiv_syph(hiv, syphilis, rel_sus_hiv_syph=100, rel_trans_hiv_syph=100)

    sim = ss.Sim(pars, copy_inputs=False).run()

    sim.analyzers.trackvalues.plot(agents)

    return sim


if __name__ == '__main__':
    s0 = test_hiv()
    s1 = test_hiv_syph()
    plt.show()
