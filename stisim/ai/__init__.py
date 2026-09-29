"""Bootstrap for the STIsim Claude Code plugin.

The plugin lives under ``stisim/ai/plugin/`` and is activated inside
Claude Code via the marketplace-add-local mechanism. Run
``python -m stisim.ai install`` for the activation commands.
"""

from pathlib import Path

PLUGIN_ROOT = Path(__file__).parent / "plugin"
