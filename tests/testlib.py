import stisim as sti


def build_testing_sim(diseases: list,
                      n_agents: int = 500,
                      demographics: list = None,
                      interventions: list = None,
                      networks: list = None):
    sim = sti.Sim(
        start=1990,
        dur=40,
        n_agents=n_agents,
        diseases=diseases,
        networks=networks,
        demographics=demographics,
        interventions=interventions
    )
    return sim
