'''
Load data
'''

#%% Housekeeping
import numpy as np
import pandas as pd
import sciris as sc


__all__ = ['get_age_distribution', 'get_rates']


thisdir = sc.thispath(__file__)
filesdir = thisdir / 'files'
files = sc.objdict()
files.age_dist_sex = 'populations_by_sex.obj'
files.birth = 'cbr.csv'
files.asfr = 'asfr.csv'
files.death = 'deaths.csv'
files.migration = 'migration.csv'


def get_age_distribution(location=None, year=None, datafolder=None):
    """ Load age distribution for a given location & year"""
    if datafolder is None: datafolder = filesdir
    filepath = sc.makefilepath(filename= f'{location}_age_{year}.csv', folder=datafolder)
    try:
        raw_df = pd.read_csv(filepath)
    except Exception as E:
        errormsg = f'Could not locate datafile for age distribution.'
        raise ValueError(errormsg) from E

    # Remap column names
    raw_df = raw_df.rename(columns={'Time': 'year', 'AgeStart': 'age', 'Value': 'value'})
    if 'year' not in raw_df.columns:
        raw_df['year'] = year
    raw_df = raw_df[['year', 'age', 'value']]

    return raw_df


def get_rates(location, which, datafolder=None):
    """ Load birth/death/fertility/migration rates for a given location, with a row for every year (see interp_years()) """
    if datafolder is None: datafolder = filesdir
    filepath = sc.makefilepath(filename=f'{location}_{files[which]}', folder=datafolder)
    try:
        df = pd.read_csv(filepath)
    except Exception as E:
        errormsg = f'Could not locate datafile for {which}.'
        raise ValueError(errormsg) from E
    return interp_years(df)


def interp_years(df, year='Time', keys=('Sex', 'AgeStart')):
    """
    Interpolate rates given for some years (e.g. every decade, as downloaded) linearly to every year in between

    The rates are estimates for single years, and Starsim holds a year's rate until the next year in the data,
    so without this the rate for 1990 would apply until 2000.
    """
    if year not in df.columns:
        return df
    keys = [k for k in keys if k in df.columns]
    years = np.arange(df[year].min(), df[year].max() + 1)
    if len(years) == df[year].nunique(): # Already yearly
        return df
    out = []
    for vals, group in (df.groupby(keys) if keys else [((), df)]):
        group = group.drop(columns=keys).set_index(year).sort_index().reindex(years).interpolate(method='index')
        group[keys] = vals
        out.append(group.reset_index(names=year))
    return pd.concat(out, ignore_index=True)[df.columns].astype({col:df[col].dtype for col in [year] + keys})

