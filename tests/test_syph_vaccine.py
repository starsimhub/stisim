"""
Test functionality of syphilis vaccine
- Vaccinate susceptible agents
- Vaccinate infected agents
- Infect vaccinated agents
- Vaccinate pregnant women, infected or susceptible
- Vaccinate multiple times
- Vaccinate before or after re-infection 
"""

#% Imports
import starsim as ss
import stisim as sti
import pandas as pd
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import sciris as sc
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
from stisim.interventions import SyphTx


class TrackValues(ss.Analyzer):
    # Track outputs for rel_trans, rel_sus
    # Assumes no births; for diagnostic/debugging purposes only
    def init_pre(self, sim):
        super().init_pre(sim)
        self.n = len(sim.people)

        self.hiv_rel_sus = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)
        self.hiv_rel_trans = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)

        self.syph_rel_sus = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)
        self.syph_rel_trans = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)
        self.rel_trans_maternal = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)

        self.syph_immunity = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)
        self.syph_rel_sus_immunity = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)
        self.syph_rel_trans_immunity = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)
        self.syph_rel_trans_immunity_maternal = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)

        self.syph_state = np.empty((sim.npts, self.n))
        self.ti_secondary = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)
        self.ti_latent = np.empty((sim.npts, self.n), dtype=ss.dtypes.float)

    @property
    def has_hiv(self):
        return 'hiv' in self.sim.diseases

    @property
    def has_syph(self):
        return isinstance(self.sim.diseases.get('syphilis'), sti.Syphilis)

    def step(self):
        sim = self.sim
        if len(sim.people.alive.uids):
            if self.has_hiv:
                self.hiv_rel_sus[sim.ti, :self.n] = sim.diseases.hiv.rel_sus.values[:self.n]
                self.hiv_rel_trans[sim.ti, :self.n] = sim.diseases.hiv.rel_trans.values[:self.n]

            if self.has_syph:
                self.syph_rel_sus[sim.ti, :self.n] = sim.diseases.syphilis.rel_sus.values[:self.n]
                self.syph_rel_trans[sim.ti, :self.n] = sim.diseases.syphilis.rel_trans.values[:self.n]
                self.rel_trans_maternal[sim.ti, :self.n] = sim.diseases.syphilis.rel_trans_maternal.values[:self.n]

                # Immunity Inf and trans
                self.syph_immunity[sim.ti, :self.n] = sim.interventions.syph_vaccine.immunity.values[:self.n]
                self.syph_rel_sus_immunity[sim.ti, :self.n] = sim.interventions.syph_vaccine.rel_sus_immunity.values[:self.n]
                self.syph_rel_trans_immunity[sim.ti, :self.n] = sim.interventions.syph_vaccine.rel_trans_immunity.values[:self.n]
                self.syph_rel_trans_immunity_maternal[sim.ti, :self.n] = sim.interventions.syph_vaccine.rel_trans_immunity_maternal.values[:self.n]

                # State of each agent
                susceptible_agents = sim.diseases.syphilis.susceptible.uids
                primary_agents = sim.diseases.syphilis.primary.uids
                secondary_agents = sim.diseases.syphilis.secondary.uids
                tertiary_agents = sim.diseases.syphilis.tertiary.uids
                latent_agents = sim.diseases.syphilis.latent.uids
                dead = sim.people.dead.uids
                self.syph_state[sim.ti, sim.people.auids] = -1
                self.syph_state[sim.ti, susceptible_agents] = 0
                self.syph_state[sim.ti, primary_agents] = 1
                self.syph_state[sim.ti, secondary_agents] = 2
                self.syph_state[sim.ti, tertiary_agents] = 3
                self.syph_state[sim.ti, latent_agents] = 4

                # Duration of infection
                self.ti_secondary[sim.ti, :self.n] = sim.diseases.syphilis.ti_secondary.values[:self.n]
                self.ti_latent[sim.ti, :self.n] = sim.diseases.syphilis.ti_latent.values[:self.n]

    def plot(self, agents: dict):
        """
        :param agents: Dictionary of events per agent {'agent_description':[('event_type', ti),...]}
        :return: Matplotlib figure
        """
        colors = {0: 'tab:green', 1: 'tab:blue', 2: 'tab:orange', 3: 'tab:red', 4: 'tab:grey', -1: 'black'}
        def plot_with_events(ax, x, y, agents, title, colors):
            for idx, agent in enumerate(agents):
                syph_state_array = self.syph_state[:, idx].astype(int)
                h = ax.scatter(x, y[:, idx], c= [colors[k] for k in syph_state_array])

            x_ev = []
            y_ev = []
            colors = []
            labels = []
            for i, events in enumerate(agents.values()):
                for event in events:
                    x_ev.append(self.sim.timevec[event[1]])
                    y_ev.append(y[event[1], i])
                    if event[0] == 'syphilis_infection':
                        colors.append('red')
                    if event[0] == 'hiv_infection':
                        colors.append('orange')
                    if event[0] == 'syph_vaccine':
                        colors.append('blue')
                    if event[0] == 'pregnant':
                        colors.append('yellow')
                    if event[0] == 'syphilis_treatment':
                        colors.append('black')
            ax.scatter(x_ev, y_ev, marker='*', color=colors, edgecolor='red', s=150, linewidths=0.5, zorder=100)
            ax.set_title(title)
            # ax.set_ylim([0, 1])
            return h

        if self.has_syph:
            fig, ax = plt.subplots(2, 4, figsize=(16, 8))
        else:
            fig, ax = plt.subplots(1, 2)

        ax = ax.ravel()
        h = plot_with_events(ax[0], self.sim.timevec, self.syph_immunity, agents, 'Immunity', colors)
        h = plot_with_events(ax[1], self.sim.timevec, self.syph_rel_sus_immunity, agents, 'Syphilis rel_sus_immunity', colors)
        h = plot_with_events(ax[2], self.sim.timevec, self.syph_rel_trans_immunity, agents, 'Syphilis rel_trans_immunity', colors)
        h = plot_with_events(ax[3], self.sim.timevec, self.syph_rel_trans_immunity_maternal, agents, 'Syphilis rel_trans_immunity_maternal', colors)

        if self.has_syph:
            h = plot_with_events(ax[5], self.sim.timevec, self.syph_rel_sus, agents, 'Syphilis rel_sus', colors)
            h = plot_with_events(ax[6], self.sim.timevec, self.syph_rel_trans, agents, 'Syphilis rel_trans', colors)
            h = plot_with_events(ax[7], self.sim.timevec, self.rel_trans_maternal, agents, 'Syphilis maternal rel_trans', colors)

        # for axis in ax:
        #   axis.set_xlim([2020, 2024])

        # Add Legend
        ax[4].axis('off')
        infection = Line2D([0], [0], label='Syphilis Infection', linestyle='', marker='*', color='red')
        hiv_infection = Line2D([0], [0], label='HIV Infection', linestyle='', marker='*', color='tab:orange')
        vaccination = Line2D([0], [0], label='Vaccination', linestyle='', marker='*', color='blue')
        pregnant = Line2D([0], [0], label='Pregnancy', linestyle='', marker='*', color='yellow')
        treatment = Line2D([0], [0], label='Treatment', linestyle='', marker='*', color='black')
        patches = []
        for state, label in zip([0, 1, 2, 3, 4], ['susceptible', 'primary', 'secondary', 'tertiary', 'latent']):
            state_patch = mpatches.Patch(color=colors[state], label=label)
            patches.append(state_patch)
        ax[4].legend(frameon=False, handles=[infection, hiv_infection, vaccination, pregnant, treatment] + patches)
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
        self.syphilis_treatment = defaultdict(list)
        self.pregnant = defaultdict(list)

        if events:
            for uid, event, ti in events:
                if event == 'hiv_infection':
                    self.hiv_infections[ti].append(uid)
                elif event == 'syphilis_infection':
                    self.syphilis_infections[ti].append(uid)
                elif event == 'syphilis_treatment':
                    self.syphilis_treatment[ti].append(uid)
                elif event == 'syph_vaccine':
                    self.syph_vaccine[ti].append(uid)
                elif event == 'pregnant':
                    self.pregnant[ti].append(uid)
                else:
                    raise Exception(f'Unknown event "{event}"')

    def administer_vaccine(self, uids):
        if len(uids):
            self.sim.interventions.syph_vaccine.vaccinate(self.sim, ss.uids(uids), update_immunity_by_vaccination=False)
            # if (self.sim.interventions.syph_vaccine.doses[ss.uids(uids)] > 1).any():
            #     uids_boost = ss.uids(uids)[self.sim.interventions.syph_vaccine.doses[ss.uids(uids)] > 1]
            #     self.sim.interventions.syph_vaccine.boost_immunity_by_vaccination(self.sim, ss.uids(uids_boost))

    def administer_treatment(self, uids):
        if len(uids):
            self.sim.interventions.treat.change_states('syphilis', ss.uids(uids))

    def set_pregnancy(self, uids):
        self.sim.demographics.pregnancy.pregnant[ss.uids(uids)] = True
        self.sim.demographics.pregnancy.ti_pregnant[ss.uids(uids)] = self.sim.ti

    def step(self):
        sim = self.sim
        if len(sim.people.alive.uids):
            if 'hiv' in sim.diseases:
                self.sim.diseases.hiv.set_prognoses(ss.uids(self.hiv_infections[sim.ti]))
            if 'syphilis' in sim.diseases:
                self.sim.diseases.syphilis.set_prognoses(ss.uids(self.syphilis_infections[sim.ti]))
                self.administer_vaccine(self.syph_vaccine[sim.ti])
                self.administer_treatment(self.syphilis_treatment[sim.ti])

            # Set pregnancies:
            self.set_pregnancy(self.pregnant[sim.ti])


def test_syph_vacc():
    # AGENTS
    agents = sc.odict()
    # agents['Gets vaccine at start of infection'] = [('syphilis_infection', 1), ('syph_vaccine', 1)]
    # agents['Infection, vaccine after infection'] = [('syphilis_infection', 2), ('syph_vaccine', 2)]
    agents['No infection, vaccine'] = [('syph_vaccine', 20)]
    # agents['Infection, no vaccine'] = [('syphilis_infection', 50)]
    # agents['Infection after vaccine'] = [('syphilis_infection', 30), ('syph_vaccine', 20)]
    # agents['Pregnancy'] = [('syphilis_infection', 15), ('pregnant', 50), ('syph_vaccine', 2)]
    # agents['HIV Infection'] = [('syphilis_infection', 85), ('syph_vaccine', 90), ('hiv_infection', 150)]
    # agents['HIV Infection'] = [ ('syph_vaccine', 90), ('hiv_infection', 100), ('syphilis_infection', 120)]
    # agents['No HIV Infection'] = [('syphilis_infection', 10), ('syph_vaccine', 11)]
    # agents['Reinfection'] = [('syphilis_infection', 10), ('syphilis_treatment', 20)] #, ('syph_vaccine', 25), ('syphilis_infection', 30)]
    # agents['Reinfection'] = [('syph_vaccine', 1), ('syphilis_infection', 10), ('syphilis_treatment', 20), ('syphilis_infection', 100)]
    # agents['Reinfection'] = [('syph_vaccine', 1), ('syphilis_infection', 5), ('syph_vaccine', 7), ('syphilis_treatment', 20), ('syphilis_infection', 30)]

    events = []
    for i, x in enumerate(agents.values()):
        for y in x:
            events.append((i,) + y)

    pars = {}
    pars['n_agents'] = len(agents)
    pars['start'] = 2020
    pars['stop'] = 2050
    pars['dt'] = 1 / 12
    syph = sti.Syphilis(init_prev=0, beta={'structuredsexual': [0, 0], 'maternal': [0, 0]})
    hiv = sti.HIV(init_prev=0, beta={'structuredsexual': [0, 0], 'maternal': [0, 0]})
    pars['diseases'] = [syph, hiv]
    pars['networks'] = [sti.StructuredSexual(), ss.MaternalNet()]
    pars['demographics'] = [ss.Pregnancy(fertility_rate=0), ss.Deaths(death_rate=0)]

    # Add Vaccine Intervention          
    def vaccine_eligible(sim):
        eligible = sim.people.age < 0 # Noone is eligible
        return eligible

    syph_vaccine = sti.interventions.SyphVaccine(
        start_year=1981,
        eligibility=vaccine_eligible,
        first_dose_efficacy=0.3,
        second_dose_efficacy=0.5,
        third_dose_efficacy=0.8,
        dose_interval=7,
        # daily_num_doses=500,
        p_second_dose=1,
        p_third_dose=1,
        dur_reach_peak=0.5,  # Reaches efficacy after 6 months
        dur_protection=5,  # Assume 5 years
    )

    # # Add Syphilis Treatment
    treat = SyphTx(rel_treat_prob=0, name='treat', label='treat')

    pars['interventions'] = [PerformTest(events), treat, syph_vaccine]

    output = TrackValues()
    pars['analyzers'] = output

    sim = ss.Sim(pars, copy_inputs=False).run()
    fig = output.plot(agents)
    fig.show()
    return sim


if __name__ == '__main__':
    s0 = test_syph_vacc()
