"""
Reproduce the docs CI build locally to catch render-time errors before they
hit main.

This test runs `quarto render` over the full docs site, exactly like the
publish_docs workflow does. It is slow (minutes) and depends on quarto being
installed, so it's marked manual: pytest skips it unless the
``RUN_DOCS_BUILD`` environment variable is set.

It also clears ``docs/_inv/`` first, so the render matches CI's fresh checkout.
That directory is a gitignored cache of the external inventories quartodoc
fetches for interlinks (python, numpy, pandas, matplotlib, sciris, starsim). A
stale entry means the local build never touches the network for that source,
so an unreachable one passes locally and fails on CI — exactly what happened at
the v1.6.1 release, when docs.scipy.org went down between the local check and
the publish (#594).

Run before merging any rcX.X.X branch into main:

    RUN_DOCS_BUILD=1 pytest tests/devtests/devtest_docs.py -s

Skip otherwise — CI will run it on push/merge.
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest


DOCS_DIR = Path(__file__).resolve().parent.parent.parent / 'docs'
INV_DIR = DOCS_DIR / '_inv'  # quartodoc's cache of external interlinks inventories


def _quarto_available():
    return shutil.which('quarto') is not None


@pytest.mark.skipif(
    not os.environ.get('RUN_DOCS_BUILD'),
    reason='Set RUN_DOCS_BUILD=1 to run the docs render test (slow)',
)
@pytest.mark.skipif(not _quarto_available(), reason='quarto not installed')
def test_docs_render():
    """Render the docs site end-to-end; any tutorial/example error will fail."""
    # Force interlinks to re-fetch, so a dead source fails here rather than on CI
    shutil.rmtree(INV_DIR, ignore_errors=True)

    result = subprocess.run(
        ['quarto', 'render'],
        cwd=DOCS_DIR,
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print('---- quarto stdout ----')
        print(result.stdout)
        print('---- quarto stderr ----')
        print(result.stderr)
    assert result.returncode == 0, (
        f'quarto render failed with exit {result.returncode}. '
        f'Inspect output above for the failing notebook/cell.'
    )


if __name__ == '__main__':
    os.environ['RUN_DOCS_BUILD'] = '1'
    test_docs_render()
    print('Docs build OK.')
