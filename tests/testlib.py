import stisim as sti


def build_testing_sim(diseases: list,
                      n_agents: int = 500,
                      demographics: list = None,
                      interventions: list = None,
                      networks: list = None,
                      analyzers: list = None,
                      duration: int = 40):
    sim = sti.Sim(
        start=1990,
        dur=duration,
        n_agents=n_agents,
        diseases=diseases,
        networks=networks,
        demographics=demographics,
        interventions=interventions,
        analyzers=analyzers
    )
    return sim
