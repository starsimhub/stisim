"""CLI for the STIsim Claude Code plugin.

Usage:

    python -m stisim.ai install     # print the two slash commands to run in Claude Code
    python -m stisim.ai uninstall   # print the removal commands
    python -m stisim.ai status      # check whether the plugin is currently registered
"""

import argparse
import json
import sys
from pathlib import Path

from . import PLUGIN_ROOT

MARKETPLACE_NAME = "stisim-local"
PLUGIN_SPEC = f"stisim@{MARKETPLACE_NAME}"
INSTALLED_STATE = Path.home() / ".claude" / "plugins" / "installed_plugins.json"


def _verify():
    manifest = PLUGIN_ROOT / ".claude-plugin" / "plugin.json"
    marketplace = PLUGIN_ROOT / ".claude-plugin" / "marketplace.json"
    for path in (manifest, marketplace):
        if not path.exists():
            sys.exit(f"Plugin file missing: {path}. Reinstall stisim.")


def install(_args):
    _verify()
    print(
        "To activate the STIsim Claude Code plugin, run these two commands\n"
        "inside a Claude Code session (once per user, not per project):\n"
    )
    print(f"  /plugin marketplace add {PLUGIN_ROOT}")
    print(f"  /plugin install {PLUGIN_SPEC}")
    print(
        "\nAfter activation, the /stisim:* skills will be available in every\n"
        "Claude Code session. Upgrading stisim via pip will update the plugin\n"
        "automatically on Claude Code's next marketplace refresh."
    )


def uninstall(_args):
    print(
        "To remove the STIsim Claude Code plugin, run these commands inside\n"
        "a Claude Code session:\n"
    )
    print(f"  /plugin uninstall {PLUGIN_SPEC}")
    print(f"  /plugin marketplace remove {MARKETPLACE_NAME}")


def status(_args):
    if not INSTALLED_STATE.exists():
        print(f"No Claude Code plugin state found at {INSTALLED_STATE}.")
        return
    data = json.loads(INSTALLED_STATE.read_text())
    matches = [
        (key, record)
        for key, records in data.get("plugins", {}).items()
        if key.startswith("stisim@")
        for record in records
    ]
    if not matches:
        print("STIsim plugin not currently registered.")
        print("Run 'python -m stisim.ai install' for activation instructions.")
        return
    for key, record in matches:
        version = record.get("version", "unknown")
        path = record.get("installPath", "unknown")
        print(f"Installed: {key}  version {version}  at {path}")


def main():
    parser = argparse.ArgumentParser(prog="python -m stisim.ai")
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("install", help="print activation commands").set_defaults(func=install)
    sub.add_parser("uninstall", help="print removal commands").set_defaults(func=uninstall)
    sub.add_parser("status", help="check current registration").set_defaults(func=status)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
