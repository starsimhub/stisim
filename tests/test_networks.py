"""
Network dynamics validation: verify that network parameters (concurrency, pair
formation, relationship duration, debut age, MSM) affect behaviour as expected.
"""
import stisim as sti
import numpy as np
import sciris as sc


def test_age_differences(n_agents=10000):
    """
    Sexual network with default relationship type distribution: run for 1 month and
    analyze partner age differences. Produces a scatterplot of female vs male partner
    ages with a best-fit regression line, and prints the mean and standard deviation
    of the age gap.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    new_upper_age = 8  # Adam, modify this number and try rerunning, noticing the plot and statistical differences

    sdA = 3
    sdB = 2
    sdC = 1
    age_diff_pars_older = {'teens': [(new_upper_age, sdA), (new_upper_age, sdA), (new_upper_age, sdC)],
                           'young': [(new_upper_age, sdA), (new_upper_age, sdA), (new_upper_age, sdB)],
                           'adult': [(new_upper_age, sdA), (new_upper_age, sdA), (new_upper_age, sdB)]}
    network = sti.StructuredSexual(age_diff_pars=age_diff_pars_older)
    sim = sti.Sim(n_agents=n_agents, networks=[network], dur=1/12)
    sim.run()

    nw = sim.networks.structuredsexual
    male_ages = np.array(nw.edges.age_p1, dtype=float)
    female_ages = np.array(nw.edges.age_p2, dtype=float)

    # Restrict both axes to 15–60 years
    mask = (female_ages >= 15) & (female_ages <= 60) & (male_ages >= 15) & (male_ages <= 60)
    female_ages = female_ages[mask]
    male_ages = male_ages[mask]

    assert len(female_ages) > 0, "No partnerships found after 1 month — check network setup"

    age_diffs = male_ages - female_ages
    mean_diff = np.mean(age_diffs)
    std_diff = np.std(age_diffs)
    print(f"Average partner age difference (male − female): {mean_diff:.2f} years")
    print(f"Standard deviation of partner age difference:   {std_diff:.2f} years")

    # Linear regression via numpy
    coeffs = np.polyfit(female_ages, male_ages, 1)
    slope, intercept = coeffs

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(female_ages, male_ages, alpha=0.4, s=8, label='Partnerships')
    x_fit = np.array([15.0, 60.0])
    ax.plot(x_fit, np.polyval(coeffs, x_fit), 'r-', linewidth=2,
            label=f'Linear fit  slope={slope:.2f}')
    ax.set_xlabel('Female partner age (years)')
    ax.set_ylabel('Male partner age (years)')
    ax.set_title('Relationship partner ages')
    ax.set_xlim(15, 60)
    ax.set_ylim(15, 60)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    outfile = sc.thisdir(__file__, 'relationship_age_differences.png')
    plt.savefig(outfile, dpi=100)
    plt.close(fig)
    print(f"Scatterplot saved to {outfile}")

    return sim


if __name__ == '__main__':
    test_age_differences()