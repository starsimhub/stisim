"""
Development tests for the chlamydia model
"""
 
# Imports
import sciris as sc
import numpy as np
import pylab as pl
 
def plot_bacterial_load(n_years=2, dt=1/52, init_load=1, peak_load=10e7):

    time_to_peak = 8/52
    half_life = 2.5/52
    ti_peak = int(time_to_peak/dt)
    growth_rate = np.log(peak_load/init_load)/(ti_peak-1)
    decay_rate = np.log(2) / (half_life/dt)
    nt = int(n_years/dt)

    def f1(t1, init_load, peak_load, ti_peak):
        return init_load + (peak_load - init_load) / (ti_peak - 1) * t1

    def f11(t1, init_val, growth_rate):
        return init_val * np.exp(growth_rate * t1)

    def f2(t2, init_val, decay_rate):
        return init_val * np.exp(-decay_rate * t2)

    t1 = np.arange(ti_peak)
    t2 = np.arange(1, nt-ti_peak)
    # y1 = f1(t1, init_load, peak_load, ti_peak)
    y1 = f11(t1, init_val=init_load, growth_rate=growth_rate)
    y2 = f2(t2, init_val=peak_load, decay_rate=decay_rate)
    y = np.maximum(np.concatenate([y1, y2]), 1)
    t = np.concatenate([t1, t2+ti_peak])

    sc.options(fontsize=25)

    fig, axes = pl.subplots(1, 3, figsize=(25, 8))
    ax = axes[0]
    ax.semilogy(t, y)
    ax.axhline(y=50, color='k', ls='--')
    inf_start = sc.findfirst(y > 50)
    inf_end = sc.findfirst(y[8:] < 50) + 8
    ax.axvline(x=inf_start, color='k', ls='--')
    ax.axvline(x=inf_end, color='k', ls='--')
    ax.set_title('Chlamydial load')
    ax.set_xlabel('Weeks since infection')
    ax.set_ylabel('IFUs/mL of mucosal cervicovaginal secretion')

    max_beta = 0.25
    k = 0.25
    ytheor = np.linspace(0, 7, 100)
    yvals = 10**ytheor
    beta = 2*max_beta / (1 + np.exp(-k*np.log10(yvals))) - max_beta
    ax = axes[1]
    ax.semilogx(yvals, beta)
    ax.set_title('Transmissibility and bacterial load')
    ax.set_xlabel('IFUs/mL of mucosal cervicovaginal secretion')
    ax.set_ylabel('Probability of transmission ')

    ax = axes[2]
    betat = 2*max_beta / (1 + np.exp(-k*np.log10(y))) - max_beta
    ax.plot(t, betat)
    ax.set_title('Transmissibility over time')
    ax.set_xlabel('Weeks since infection')
    ax.set_ylabel('Probability of transmission ')


    sc.figlayout()
    pl.savefig(f"chlamydia.png", dpi=100)

    return t, y




if __name__ == '__main__':

    t, y = plot_bacterial_load()
