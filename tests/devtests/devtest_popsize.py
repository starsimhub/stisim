"""
Test that the population size doesn't depend on the number of agents in the simulation
"""

# Imports
import sciris as sc
import starsim as ss
import stisim as sti
import pandas as pd
import numpy as np
import pylab as pl
 
def make_sim(n_agents=500, dt=1):

    fertility_rates = {'fertility_rate': pd.read_csv('test_data/zimbabwe_asfr.csv')}
    pregnancy = ss.Pregnancy(pars=fertility_rates)
    death_rates = {'death_rate': pd.read_csv('test_data/zimbabwe_deaths.csv'), 'units': 1}
    death = ss.Deaths(death_rates)
    ppl = ss.People(n_agents, age_data=pd.read_csv('test_data/zimbabwe_age.csv'))

    sim = ss.Sim(
        dt=dt,
        start=1990,
        total_pop=9980999,
        dur=35,
        people=ppl,
        demographics=[pregnancy, death],
    )

    return sim

def test_n_agents():
    sc.heading('Test pop sizes with varying n_agents')
    results = dict()
    sims = sc.autolist()
    n_agent_list = [5e3, 10e3, 20e3]
    for n_agents in n_agent_list:
        sim = make_sim(n_agents=n_agents, dt=1)
        sim.run()
        results[n_agents] = sim.results.n_alive
        sims += sim

    fig, ax = pl.subplots(1, 1)
    for n_agents in n_agent_list:
        ax.plot(sim.timevec, results[n_agents], label=int(n_agents))
    ax.legend()

    return sims


if __name__ == '__main__':
    # sims = test_n_agents()


    import matplotlib.pyplot as plt
    import numpy as np

    # Define time in months
    months = np.arange(0, 37, 1)

    # Define the cumulative percentage for each scenario
    cumulative_percent1 = np.piecewise(months,
                                       [months <= 2, (months > 2) & (months <= 8), months > 8],
                                       [lambda months: 25*months,
                                        lambda months: 50 + (25/6)*(months-2),
                                        lambda months: 75 + (25/27)*(months-8)])

    cumulative_percent2 = np.piecewise(months,
                                       [months <= 2, (months > 2) & (months <= 8), months > 8],
                                       [lambda months: 10*months,
                                        lambda months: 20 + (60/6)*(months-2),
                                        lambda months: 80 + (20/27)*(months-8)])

    cumulative_percent3 = np.piecewise(months,
                                       [months <= 2, (months > 2) & (months <= 8), months > 8],
                                       [lambda months: 5*months,
                                        lambda months: 10 + (80/6)*(months-2),
                                        lambda months: 90 + (10/27)*(months-8)])

    # Plot the cumulative percentages
    plt.figure(figsize=(10, 6))
    plt.plot(months, cumulative_percent1, label="Scenario 1")
    plt.plot(months, cumulative_percent2, label="Scenario 2")
    plt.plot(months, cumulative_percent3, label="Scenario 3")

    # Shade areas for different periods
    plt.axvspan(0, 2, color='gray', alpha=0.3, label="Primary (0-2 months)")
    plt.axvspan(2, 6, color='blue', alpha=0.1, label="Secondary (2-6 months)")
    plt.axvspan(6, 18, color='green', alpha=0.1, label="Early Latent (6-18 months)")
    plt.axvspan(18, 36, color='red', alpha=0.1, label="Late Latent (18-36 months)")

    # Add labels and title
    plt.xlabel("Time (months)")
    plt.ylabel("Cumulative Percentage of Transmissions")
    plt.title("Cumulative Percentage of Transmissions Over Time")
    plt.legend()
    plt.grid(True)
    plt.show()