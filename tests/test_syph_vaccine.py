import starsim as ss
import sys
x = 'c:\\users\\alina.muellenmeister\\documents\\github\\gavi-outbreaks'
y = 'c:\\users\\alina.muellenmeister\\documents\\github\\syphilis_analyses'
if x in sys.path:
    sys.path.remove(x)
if y in sys.path:
    sys.path.remove(y)
import stisim as sti
import pandas as pd
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import sciris as sc


class TrackValues(ss.Analyzer):
    # Track outputs for viral load and CD4 counts
    # Assumes no births; for diagnostic/debugging purposes only
    def init_pre(self, sim):
        super().init_pre(sim)
        self.n = len(sim.people)

        self.hiv_rel_sus = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)
        self.hiv_rel_trans = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)

        self.syph_rel_sus = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)
        self.syph_rel_trans = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)

        self.syph_immunity_inf = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)
        self.syph_immunity_trans = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)

    @property
    def has_hiv(self):
        return 'hiv' in self.sim.diseases

    @property
    def has_syph(self):
        return isinstance(self.sim.diseases.get('syphilis'), sti.Syphilis)

    def apply(self, sim):

        if self.has_hiv:
            self.hiv_rel_sus[sim.ti, :self.n] = sim.diseases.hiv.rel_sus.values[:self.n]
            self.hiv_rel_trans[sim.ti, :self.n] = sim.diseases.hiv.rel_trans.values[:self.n]

        if self.has_syph:
            self.syph_rel_sus[sim.ti, :self.n] = sim.diseases.syphilis.rel_sus.values[:self.n]
            self.syph_rel_trans[sim.ti, :self.n] = sim.diseases.syphilis.rel_trans.values[:self.n]

        self.syph_immunity_inf[sim.ti, :self.n] = sim.interventions.syph_vaccine.immunity_inf.values[:self.n]
        self.syph_immunity_trans[sim.ti, :self.n] = sim.interventions.syph_vaccine.immunity_trans.values[:self.n]


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
                    x_ev.append(self.sim.yearvec[event[1]])
                    y_ev.append(y[event[1], i])
            ax.scatter(x_ev, y_ev, marker='*', color='yellow', edgecolor='red', s=100, linewidths=0.5, zorder=100)
            ax.set_title(title)
            return h

        if self.has_syph:
            fig, ax = plt.subplots(2, 2)
        else:
            fig, ax = plt.subplots(1, 2)

        ax = ax.ravel()

        h = plot_with_events(ax[0], self.sim.yearvec, self.syph_immunity_inf, agents, 'Syphilis imm inf')
        h = plot_with_events(ax[1], self.sim.yearvec, self.syph_immunity_trans, agents, 'Syphilis imm trans')

        if self.has_syph:
            h = plot_with_events(ax[2], self.sim.yearvec, self.syph_rel_sus, agents, 'Syphilis rel_sus')
            h = plot_with_events(ax[3], self.sim.yearvec, self.syph_rel_trans, agents, 'Syphilis rel_trans')

        fig.legend(h, agents.keys(), loc='upper right', bbox_to_anchor=(1.1, 1))

        return fig


class PerformTest(ss.Intervention):

    def __init__(self, events=None):
        """
        :param events: List of (uid, 'event', ti) to apply events to an agent
        """
        super().__init__()
        self.hiv_infections = defaultdict(list)
        self.syphilis_infections = defaultdict(list)
        self.syph_vaccine = defaultdict(list)
        self.pregnant = defaultdict(list)

        if events:
            for uid, event, ti in events:
                if event == 'hiv_infection':
                    self.hiv_infections[ti].append(uid)
                elif event == 'syphilis_infection':
                    self.syphilis_infections[ti].append(uid)
                elif event == 'syph_vaccine':
                    self.syph_vaccine[ti].append(uid)
                elif event == 'pregnant':
                    self.pregnant[ti].append(uid)
                else:
                    raise Exception(f'Unknown event "{event}"')

    def administer_vaccine(self, uids):
        if len(uids):
            self.sim.interventions.syph_vaccine.vaccinate(self.sim, ss.uids(uids))


    def set_pregnancy(self, uids):
        self.sim.demographics.pregnancy.pregnant[ss.uids(uids)] = True
        self.sim.demographics.pregnancy.ti_pregnant[ss.uids(uids)] = self.sim.ti

    def apply(self, sim):
        self.administer_vaccine(self.syph_vaccine[sim.ti])
        if 'hiv' in sim.diseases:
            self.sim.diseases.hiv.set_prognoses(ss.uids(self.hiv_infections[sim.ti]))
        if 'syphilis' in sim.diseases:
            self.sim.diseases.syphilis.set_prognoses(ss.uids(self.syphilis_infections[sim.ti]))

        # Set pregnancies:
        self.set_pregnancy(self.pregnant[sim.ti])


def test_syph_vacc():
    # AGENTS
    agents = sc.odict()
    agents['No infection'] = []
    agents['Gets vaccine at start of infection'] = [('syphilis_infection', 1), ('syph_vaccine', 2)]
    agents['Gets vaccine later'] = [('syphilis_infection', 1), ('syph_vaccine', 12)]
    agents['Gets vaccine while naive'] = [('syph_vaccine', 12)]

    events = []
    for i, x in enumerate(agents.values()):
        for y in x:
            events.append((i,) + y)

    pars = {}
    pars['n_agents'] = len(agents)
    pars['start'] = 2020
    pars['end'] = 2040
    pars['dt'] = 1 / 12
    syph = sti.Syphilis(init_prev=0, beta={'structuredsexual': [0, 0], 'maternal': [0, 0]})
    pars['diseases'] = [syph]
    pars['networks'] = [sti.StructuredSexual(), ss.MaternalNet()]
    pars['demographics'] = [ss.Pregnancy(fertility_rate=0), ss.Deaths(death_rate=0)]

    # Add Vaccine Intervention
    def vaccine_eligible(sim):
        eligible = sim.people.age >= 0 # Everyone is eligible
        return eligible
    syph_vaccine = sti.SyphVaccine(
        start_year=1981,
        eligibility=vaccine_eligible,
        target_coverage=0.75,
        efficacy=0.95,  # Peak
        dur_reach_peak=6,  # Reaches efficacy after 6 months
        dur_protection=18,  # Assume 18 months
    )

    pars['interventions'] = [PerformTest(events), syph_vaccine]
    output = TrackValues()
    pars['analyzers'] = output

    sim = ss.Sim(pars, copy_inputs=False).run()
    fig = output.plot(agents)
    fig.show()
    return sim



if __name__ == '__main__':
    s0 = test_syph_vacc()
