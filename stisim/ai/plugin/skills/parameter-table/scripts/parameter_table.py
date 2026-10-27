"""
Assemble a supplementary parameter table ("Table S1") from the STIsim user-guide pages.

Each user-guide page documents its module's parameters in tables with the header
``| Parameter | Default | Description | Source |``. This script collects those tables
into one, grouped by module and section, and optionally adds one value column per
simulation, marking values that differ from the STIsim default.

Usage (command line):

    python parameter_table.py                                 # markdown to stdout
    python parameter_table.py --format csv -o parameter_table.csv
    python parameter_table.py --docs /path/to/stisim/docs     # use a local docs checkout

Usage (Python):

    import parameter_table
    rows = parameter_table.parse()
    rows = parameter_table.add_sim(rows, sim, label='Zimbabwe')
    print(parameter_table.to_markdown(rows))

Only the standard library is required unless ``add_sim`` is used.
"""
import re
import csv
import sys
import argparse
import inspect
import urllib.request
from pathlib import Path

__all__ = ['PAGES', 'parse', 'add_sim', 'to_markdown', 'to_csv', 'check_coverage']

# Pages in table order: (path relative to docs/, module kind, module name in the sim)
# A module name of None means "search every module of that kind for the parameter".
PAGES = [
    ('user_guide/diseases/hiv.md',              'diseases',   'hiv'),
    ('user_guide/diseases/syphilis.md',         'diseases',   'syph'),
    ('user_guide/diseases/chlamydia.md',        'diseases',   'ct'),
    ('user_guide/diseases/gonorrhea.md',        'diseases',   'ng'),
    ('user_guide/diseases/trichomoniasis.md',   'diseases',   'tv'),
    ('user_guide/diseases/bv.md',               'diseases',   'bv'),
    ('user_guide/diseases/gud.md',              'diseases',   'gud'),
    ('user_guide/networks/structured_sexual.md', 'networks',  'structuredsexual'),
    ('user_guide/networks/msm.md',              'networks',   None),
    ('user_guide/connectors.md',                'connectors', None),
]

HEADER = ['Parameter', 'Default', 'Description', 'Source']
NOTES = ('Distributions: `normal(mean, sd)`, `lognorm(mean, sd)`, `uniform(low, high)`; `[F, M]` gives female and male '
         'values. Sources: a citation is the published source recorded for the default; "Derived" is computed from the '
         'cited source; "Calibrated" is taken from a previous calibration; "Assumption" is a modeling assumption; '
         '"Unpublished evidence review" is an internal review not yet published; "Not recorded" means no source has '
         'been documented; "Not applicable" marks switches and scaling factors that are neutral by default.')
RAW_URL = 'https://raw.githubusercontent.com/starsimhub/stisim/{ref}/docs/{path}'
LINK = re.compile(r'\[([^\]]+)\]\(([^)]+)\)')


def _read(path, docs, ref):
    """ Read a page from a local docs folder, or from GitHub at the given ref """
    if docs is not None:
        return (Path(docs) / path).read_text(encoding='utf-8')
    with urllib.request.urlopen(RAW_URL.format(ref=ref, path=path)) as f:
        return f.read().decode('utf-8')


def _cells(line):
    """ Split a markdown table row into stripped cells, respecting escaped pipes """
    return [c.strip().replace('\\|', '|') for c in re.split(r'(?<!\\)\|', line.strip())[1:-1]]


def parse(docs=None, ref='main'):
    """
    Parse the parameter tables from the user-guide pages.

    Args:
        docs (str/Path): local ``docs/`` folder; if None, pages are fetched from GitHub
        ref (str): git ref (branch or tag, e.g. ``'v1.7.0'``) to fetch from when ``docs`` is None

    Returns:
        list of dicts with keys module, section, parameter, name, default, description, source, kind, modname
    """
    rows = []
    for path, kind, modname in PAGES:
        module = section = None
        in_table = False
        for line in _read(path, docs, ref).splitlines():
            if line.startswith('# '):
                module = line[2:].strip()
            elif line.startswith('## ') or line.startswith('### '):
                section = line.lstrip('#').strip()
            elif line.startswith('|'):
                cells = _cells(line)
                if cells == HEADER:
                    in_table = True
                elif in_table and not set(''.join(cells)) <= set('-: '):
                    par, default, desc, source = cells
                    name = re.search(r'`([^`]+)`', par)
                    rows.append(dict(module=module, section=section, parameter=par, name=name.group(1) if name else par,
                                     default=default, description=desc, source=source, kind=kind, modname=modname))
            else:
                in_table = False
    return rows


def _fmt(value):
    """ Compact string for a parameter value, e.g. ``ss.normal(loc=800, scale=50)`` -> ``normal(loc=800, scale=50)``; for a declared parameter (``ss.par``), its value """
    import starsim as ss
    if isinstance(value, ss.par):
        value = value.raw
    return re.sub(r'\bss\.', '', repr(value)).replace('\n', ' ')


def _candidates(sim, kind, modname):
    mods = getattr(sim, kind)
    if modname is not None:
        return [mods[modname]] if modname in mods else []
    return list(mods.values())


def _default_pars(mod):
    """ Default pars for a module's class; None if the class can't be built without arguments """
    for args in [(), (None, None)]:
        try:
            return type(mod)(*args).pars
        except Exception:
            continue
    return None


def add_sim(rows, sim, label=None):
    """
    Add a value column for an initialized sim. Values that differ from the class default are marked with ``*``.

    Args:
        rows (list): output of ``parse()`` (or of a previous ``add_sim()``)
        sim (sti.Sim): initialized sim
        label (str): column label (default: ``sim.label`` or ``'Value'``)

    Returns:
        rows, with a ``values`` dict keyed by column label
    """
    if not sim.initialized:
        sim.init()
    label = label or getattr(sim, 'label', None) or 'Value'
    defaults = {}
    for row in rows:
        vals = []
        for mod in _candidates(sim, row['kind'], row['modname']):
            if row['name'] not in mod.pars:
                continue
            if id(mod) not in defaults:
                defaults[id(mod)] = _default_pars(mod)
            value = _fmt(mod.pars[row['name']])
            dflt = defaults[id(mod)]
            changed = dflt is not None and row['name'] in dflt and _fmt(dflt[row['name']]) != value
            prefix = f'{mod.name}: ' if row['modname'] is None else ''
            vals.append(f'{prefix}{value}{"*" if changed else ""}')
        row.setdefault('values', {})[label] = '; '.join(vals) if vals else 'not in sim'
    return rows


def _labels(rows):
    return list(rows[0].get('values', {})) if rows else []


def to_markdown(rows, caption=None, notes=True):
    """
    Render rows as a markdown table grouped by module and section, with a caption and footnotes.

    Args:
        rows (list): output of ``parse()`` / ``add_sim()``
        caption (str): table caption; a generic one is used if None
        notes (bool): whether to append the notation and source footnotes
    """
    labels = _labels(rows)
    cols = ['Parameter', 'Description', 'Default'] + labels + ['Source']
    caption = caption or 'Table S1. STIsim model parameters, default values, and sources.'
    out = [f'**{caption}**', '', '| ' + ' | '.join(cols) + ' |', '|' + '---|' * len(cols)]
    group = None
    for row in rows:
        if (row['module'], row['section']) != group:
            group = (row['module'], row['section'])
            title = row['module'] if row['section'] in (None, 'Parameters') else f'{row["module"]}: {row["section"]}'
            out.append(f'| **{title}** |' + ' |' * (len(cols) - 1))
        vals = [row['values'][k] for k in labels]
        cells = [row['parameter'], row['description'], row['default']] + vals + [row['source']]
        out.append('| ' + ' | '.join(c.replace('|', '\\|') for c in cells) + ' |')
    if notes:
        out += ['', NOTES]
        if labels:
            out += ['', '\\* Value differs from the STIsim default.']
    return '\n'.join(out)


def to_csv(rows, file):
    """ Write rows to a CSV file (or file-like object); links in sources are written as ``text <url>`` """
    labels = _labels(rows)
    fields = ['module', 'section', 'name', 'description', 'default'] + labels + ['source']
    f = open(file, 'w', newline='', encoding='utf-8') if isinstance(file, (str, Path)) else file
    writer = csv.writer(f)
    writer.writerow(fields)
    for row in rows:
        source = LINK.sub(r'\1 <\2>', row['source'])
        writer.writerow([row[k] for k in fields[:5]] + [row['values'][k] for k in labels] + [source])
    if f is not file:
        f.close()


def check_coverage(rows, exclude=()):
    """
    Return parameters defined in the code but missing from the docs tables, and vice versa.

    Args:
        rows (list): output of ``parse()``
        exclude (iterable): parameter names that are intentionally undocumented (internal placeholders)

    Returns:
        dict of {module class name: (undocumented, unknown)}
    """
    import starsim as ss
    import stisim as sti
    classes = dict(hiv=sti.HIV, syph=sti.Syphilis, ct=sti.Chlamydia, ng=sti.Gonorrhea, tv=sti.Trichomoniasis,
                   bv=sti.BV, gud=sti.GUD, structuredsexual=sti.StructuredSexual)
    out = {}
    for modname, cls in classes.items():
        documented = {r['name'] for r in rows if r['modname'] == modname}
        defined = {k for k, v in cls().pars.items() if not getattr(v, 'slot_named', False)} # Not those a declared model names after a slot (e.g. "dur" for a dur=0)
        if issubclass(cls, ss.Disease): # And the keywords of a declared disease, e.g. init_prev or birth_outcomes
            not_pars = ['pars', 'name', 'label', 'kwargs', 'age_bins', 'age_range', 'sex_keys', 'include_care', 'init_prev_data'] # What a disease is made with that is not one of its parameters
            defined |= {'init_prev'} | {k for k in inspect.signature(cls).parameters if k not in not_pars}
        out[cls.__name__] = (sorted(defined - documented - set(exclude)), sorted(documented - defined))
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--docs', help='local docs/ folder (default: fetch from GitHub)')
    parser.add_argument('--ref', default='main', help='git ref to fetch from GitHub (default: main)')
    parser.add_argument('--format', choices=['md', 'csv'], default='md')
    parser.add_argument('-o', '--output', help='output file (default: stdout)')
    args = parser.parse_args()
    rows = parse(docs=args.docs, ref=args.ref)
    if args.format == 'csv':
        to_csv(rows, args.output or sys.stdout)
    else:
        text = to_markdown(rows)
        Path(args.output).write_text(text, encoding='utf-8') if args.output else print(text)
