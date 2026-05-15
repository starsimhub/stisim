"""
Reproduce the docs CI build locally to catch render-time errors before they
hit main.

This test runs `quarto render` over the full docs site, exactly like the
publish_docs workflow does. It is slow (minutes) and depends on quarto being
installed, so it's marked manual: pytest skips it unless the
``RUN_DOCS_BUILD`` environment variable is set.

Run before merging any rcX.X.X branch into main:

    RUN_DOCS_BUILD=1 pytest tests/test_docs.py -s

Skip otherwise — CI will run it on push/merge.
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest


DOCS_DIR = Path(__file__).resolve().parent.parent / 'docs'


def _quarto_available():
    return shutil.which('quarto') is not None


@pytest.mark.skipif(
    not os.environ.get('RUN_DOCS_BUILD'),
    reason='Set RUN_DOCS_BUILD=1 to run the docs render test (slow)',
)
@pytest.mark.skipif(not _quarto_available(), reason='quarto not installed')
def test_docs_render():
    """Render the docs site end-to-end; any tutorial/example error will fail."""
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
