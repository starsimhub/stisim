"""
STIsim HIV Shiny Web App - Python Backend
Python wrapper functions for HIV simulation
"""

import sys
import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import sciris as sc
import starsim as ss
import stisim as sti

def run_hiv_simulation(params):
    """
    Run HIV simulation with given parameters
    
    Args:
        params (dict): Dictionary of simulation parameters
        
    Returns:
        dict: Simulation results and plots
    """
    
    try:
        # Create HIV disease module
        hiv = sti.HIV(
            init_prev=params['init_prev'],
            beta_m2f=params['beta_m2f'],
            beta_m2c=params['beta_m2c'],
            beta_m2m=params['beta_m2m'],
            rel_beta_f2m=params['rel_beta_f2m'],
            eff_condom=params['eff_condom'],
            dur_acute=ss.months(params['dur_acute']),
            dur_latent=ss.years(params['dur_latent']),
            dur_falling=ss.years(params['dur_falling']),
            rel_trans_acute=params['rel_trans_acute'],
            rel_trans_falling=params['rel_trans_falling']
        )
        
        # Create demographics
        pregnancy = ss.Pregnancy(fertility_rate=10)
        death = ss.Deaths(death_rate=10)
        
        # Create networks based on type
        networks = []
        if params['network_type'] == 'structured':
            # Add PriorPartners network for recall_prior functionality
            prior_partners = sti.PriorPartners()
            sexual = sti.StructuredSexual(recall_prior=True)
            networks.extend([prior_partners, sexual])
        elif params['network_type'] == 'msm':
            msm = sti.AgeMatchedMSM()
            networks.append(msm)
        else:  # mixed
            # Add PriorPartners network for recall_prior functionality
            prior_partners = sti.PriorPartners()
            sexual = sti.StructuredSexual(recall_prior=True)
            msm = sti.AgeMatchedMSM()
            networks.extend([prior_partners, sexual, msm])
        
        # Add maternal network
        maternal = ss.MaternalNet()
        networks.append(maternal)
        
        # Create interventions
        interventions = []
        if params['include_testing']:
            testing = sti.HIVTest(test_prob_data=0.2, start=params['start'])
            interventions.append(testing)
        
        if params['include_art']:
            # Create ART coverage data
            years = np.arange(params['start'], params['start'] + params['dur'] + 1)
            art_coverage = np.linspace(0, 0.9, len(years))
            art_data = pd.DataFrame(index=years, data={'p_art': art_coverage})
            art = sti.ART(coverage_data=art_data)
            interventions.append(art)
        
        if params['include_vmmc']:
            # Create VMMC coverage data
            years = np.arange(params['start'], params['start'] + params['dur'] + 1)
            vmmc_coverage = np.linspace(0.025, 0.125, len(years))
            vmmc_data = pd.DataFrame(index=years, data={'p_vmmc': vmmc_coverage})
            vmmc = sti.VMMC(coverage_data=vmmc_data)
            interventions.append(vmmc)
        
        if params['include_prep']:
            prep = sti.Prep()
            interventions.append(prep)
        
        # Create simulation
        sim = sti.Sim(
            start=params['start'],
            dur=params['dur'],
            n_agents=params['n_agents'],
            diseases=hiv,
            networks=networks,
            demographics=[pregnancy, death],
            interventions=interventions
        )
        
        # Run simulation
        sim.run(verbose=1/12)
        
        # Extract results
        results = extract_simulation_results(sim, params)
        
        return results
        
    except Exception as e:
        raise Exception(f"Simulation failed: {str(e)}")

def extract_simulation_results(sim, params):
    """
    Extract results from simulation object
    
    Args:
        sim: STIsim simulation object
        params: Original parameters
        
    Returns:
        dict: Extracted results
    """
    
    # Basic results
    timevec = sim.timevec
    years = sim.t.yearvec
    
    # HIV prevalence
    hiv_prev = sim.results.hiv.prevalence.values if hasattr(sim.results, 'hiv') else np.zeros_like(timevec)
    
    # HIV incidence
    hiv_inc = sim.results.hiv.incidence.values if hasattr(sim.results, 'hiv') else np.zeros_like(timevec)
    
    # Population summary
    total_pop = sim.pars.n_agents
    hiv_infected = int(np.sum(sim.people.hiv.infected)) if hasattr(sim.people, 'hiv') else 0
    hiv_prevalence = hiv_infected / total_pop if total_pop > 0 else 0
    
    # ART coverage
    on_art = int(np.sum(sim.people.hiv.on_art)) if hasattr(sim.people, 'hiv') else 0
    art_coverage = on_art / hiv_infected if hiv_infected > 0 else 0
    
    # Age-specific prevalence
    age_groups = ['15-24', '25-34', '35-44', '45-54', '55+']
    age_prev = {}
    if hasattr(sim.people, 'hiv'):
        for i, age_group in enumerate(age_groups):
            age_min = 15 + i * 10
            age_max = 24 + i * 10 if i < len(age_groups) - 1 else 100
            age_mask = (sim.people.age >= age_min) & (sim.people.age <= age_max)
            hiv_mask = sim.people.hiv.infected & age_mask
            age_prev[age_group] = np.sum(hiv_mask) / np.sum(age_mask) if np.sum(age_mask) > 0 else 0
    else:
        age_prev = {age_group: 0 for age_group in age_groups}
    
    # CD4 distribution
    cd4_counts = sim.people.hiv.cd4[sim.people.hiv.infected] if hasattr(sim.people, 'hiv') else np.array([])
    
    # Network analysis
    network_stats = {}
    if hasattr(sim, 'networks') and 'structuredsexual' in sim.networks:
        net = sim.networks.structuredsexual
        if hasattr(net, 'partners'):
            network_stats['mean_partners'] = np.mean(net.partners[net.partners > 0])
            network_stats['max_partners'] = np.max(net.partners)
        else:
            network_stats['mean_partners'] = 0
            network_stats['max_partners'] = 0
    
    return {
        'timevec': timevec,
        'years': years,
        'hiv_prevalence': hiv_prev,
        'hiv_incidence': hiv_inc,
        'population_summary': {
            'total_pop': total_pop,
            'hiv_infected': hiv_infected,
            'hiv_prevalence': hiv_prevalence,
            'on_art': on_art,
            'art_coverage': art_coverage
        },
        'age_prevalence': age_prev,
        'cd4_counts': cd4_counts,
        'network_stats': network_stats,
        'parameters': params
    }

def plot_hiv_prevalence(results):
    """
    Create HIV prevalence plot
    
    Args:
        results: Simulation results dictionary
        
    Returns:
        plotly figure
    """
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=results['years'],
        y=results['hiv_prevalence'] * 100,
        mode='lines',
        name='HIV Prevalence',
        line=dict(color='red', width=3)
    ))
    
    fig.update_layout(
        title='HIV Prevalence Over Time',
        xaxis_title='Year',
        yaxis_title='Prevalence (%)',
        hovermode='x unified',
        template='plotly_white'
    )
    
    return fig

def plot_hiv_incidence(results):
    """
    Create HIV incidence plot
    
    Args:
        results: Simulation results dictionary
        
    Returns:
        plotly figure
    """
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=results['years'],
        y=results['hiv_incidence'],
        mode='lines',
        name='HIV Incidence',
        line=dict(color='orange', width=3)
    ))
    
    fig.update_layout(
        title='HIV Incidence Over Time',
        xaxis_title='Year',
        yaxis_title='New Infections',
        hovermode='x unified',
        template='plotly_white'
    )
    
    return fig

def plot_cd4_distribution(results):
    """
    Create CD4 count distribution plot
    
    Args:
        results: Simulation results dictionary
        
    Returns:
        plotly figure
    """
    
    fig = go.Figure()
    
    if len(results['cd4_counts']) > 0:
        fig.add_trace(go.Histogram(
            x=results['cd4_counts'],
            nbinsx=30,
            name='CD4 Count Distribution',
            marker_color='lightblue'
        ))
    else:
        fig.add_annotation(
            text="No CD4 data available",
            x=0.5, y=0.5,
            showarrow=False
        )
    
    fig.update_layout(
        title='CD4 Count Distribution (HIV Infected)',
        xaxis_title='CD4 Count (cells/μL)',
        yaxis_title='Frequency',
        template='plotly_white'
    )
    
    return fig

def plot_art_coverage(results):
    """
    Create ART coverage plot
    
    Args:
        results: Simulation results dictionary
        
    Returns:
        plotly figure
    """
    
    fig = go.Figure()
    
    # Calculate ART coverage over time (simplified)
    art_coverage = np.linspace(0, results['population_summary']['art_coverage'], len(results['years']))
    
    fig.add_trace(go.Scatter(
        x=results['years'],
        y=art_coverage * 100,
        mode='lines',
        name='ART Coverage',
        line=dict(color='green', width=3)
    ))
    
    fig.update_layout(
        title='ART Coverage Over Time',
        xaxis_title='Year',
        yaxis_title='Coverage (%)',
        hovermode='x unified',
        template='plotly_white'
    )
    
    return fig

def plot_age_prevalence(results):
    """
    Create age-specific prevalence plot
    
    Args:
        results: Simulation results dictionary
        
    Returns:
        plotly figure
    """
    
    fig = go.Figure()
    
    age_groups = list(results['age_prevalence'].keys())
    prevalences = list(results['age_prevalence'].values())
    
    fig.add_trace(go.Bar(
        x=age_groups,
        y=[p * 100 for p in prevalences],
        name='Age-specific Prevalence',
        marker_color='purple'
    ))
    
    fig.update_layout(
        title='HIV Prevalence by Age Group',
        xaxis_title='Age Group',
        yaxis_title='Prevalence (%)',
        template='plotly_white'
    )
    
    return fig

def plot_network_analysis(results):
    """
    Create network analysis plot
    
    Args:
        results: Simulation results dictionary
        
    Returns:
        plotly figure
    """
    
    fig = go.Figure()
    
    # Create a simple network visualization
    network_stats = results['network_stats']
    
    fig.add_trace(go.Bar(
        x=['Mean Partners', 'Max Partners'],
        y=[network_stats['mean_partners'], network_stats['max_partners']],
        name='Network Statistics',
        marker_color='teal'
    ))
    
    fig.update_layout(
        title='Sexual Network Statistics',
        xaxis_title='Metric',
        yaxis_title='Number of Partners',
        template='plotly_white'
    )
    
    return fig

def plot_parameter_sensitivity(results):
    """
    Create parameter sensitivity plot
    
    Args:
        results: Simulation results dictionary
        
    Returns:
        plotly figure
    """
    
    fig = go.Figure()
    
    # Create a simple sensitivity analysis visualization
    params = results['parameters']
    param_names = list(params.keys())
    param_values = list(params.values())
    
    # Normalize values for display
    normalized_values = []
    for val in param_values:
        if isinstance(val, (int, float)):
            if val > 1:
                normalized_values.append(val / 10)  # Scale down large values
            else:
                normalized_values.append(val * 10)  # Scale up small values
        else:
            normalized_values.append(1)
    
    fig.add_trace(go.Bar(
        x=param_names,
        y=normalized_values,
        name='Parameter Values (Normalized)',
        marker_color='lightcoral'
    ))
    
    fig.update_layout(
        title='Parameter Values',
        xaxis_title='Parameters',
        yaxis_title='Normalized Values',
        template='plotly_white',
        xaxis_tickangle=-45
    )
    
    return fig
