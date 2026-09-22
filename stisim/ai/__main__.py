"""CLI for the STIsim Claude Code plugin.

Usage:

    python -m stisim.ai install     # register the plugin with Claude Code
    python -m stisim.ai uninstall   # remove the registration
    python -m stisim.ai status      # check whether the plugin is currently registered
"""

import argparse
import json
import os
import sys
import tempfile
from pathlib import Path

from . import PLUGIN_ROOT

MARKETPLACE_NAME = "stisim-local"
PLUGIN_KEY = f"stisim@{MARKETPLACE_NAME}"
SETTINGS_PATH = Path.home() / ".claude" / "settings.json"


def _verify_plugin():
    for name in ("plugin.json", "marketplace.json"):
        path = PLUGIN_ROOT / ".claude-plugin" / name
        if not path.exists():
            sys.exit(f"Plugin file missing: {path}. Reinstall stisim.")


def _read_settings() -> dict:
    if not SETTINGS_PATH.exists():
        return {}
    try:
        return json.loads(SETTINGS_PATH.read_text())
    except json.JSONDecodeError as exc:
        sys.exit(f"Cannot parse {SETTINGS_PATH}: {exc}. Fix the file and retry.")


def _write_settings(settings: dict) -> None:
    SETTINGS_PATH.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(
        prefix="settings.json.", dir=str(SETTINGS_PATH.parent)
    )
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(settings, f, indent=2)
            f.write("\n")
        os.replace(tmp_path, SETTINGS_PATH)
    except Exception:
        Path(tmp_path).unlink(missing_ok=True)
        raise


def install(_args):
    _verify_plugin()
    settings = _read_settings()
    enabled = settings.setdefault("enabledPlugins", {})
    marketplaces = settings.setdefault("extraKnownMarketplaces", {})

    desired_marketplace = {
        "source": {"source": "directory", "path": str(PLUGIN_ROOT)}
    }
    changed = (
        enabled.get(PLUGIN_KEY) is not True
        or marketplaces.get(MARKETPLACE_NAME) != desired_marketplace
    )

    enabled[PLUGIN_KEY] = True
    marketplaces[MARKETPLACE_NAME] = desired_marketplace

    if not changed:
        print(f"STIsim plugin already registered in {SETTINGS_PATH}.")
    else:
        _write_settings(settings)
        print(f"STIsim plugin registered in {SETTINGS_PATH}.")
        print(f"Marketplace: {MARKETPLACE_NAME} -> {PLUGIN_ROOT}")
        print(f"Plugin key:  {PLUGIN_KEY}")

    print(
        "\nReload your Claude Code session to pick up the change:\n"
        "  - VS Code / Positron:  Command Palette -> Developer: Reload Window\n"
        "  - CLI:                 exit and relaunch `claude`\n"
        "\nAfter reload, /stisim:* skills are available in every Claude Code session."
    )


def uninstall(_args):
    settings = _read_settings()
    enabled = settings.get("enabledPlugins", {})
    marketplaces = settings.get("extraKnownMarketplaces", {})

    removed_plugin = enabled.pop(PLUGIN_KEY, None) is not None
    removed_marketplace = marketplaces.pop(MARKETPLACE_NAME, None) is not None

    if not (removed_plugin or removed_marketplace):
        print("STIsim plugin is not currently registered; nothing to remove.")
        return

    _write_settings(settings)
    print(f"STIsim plugin removed from {SETTINGS_PATH}.")
    print("Reload your Claude Code session for the change to take effect.")


def status(_args):
    if not SETTINGS_PATH.exists():
        print(f"No Claude Code settings file at {SETTINGS_PATH}.")
        return
    settings = _read_settings()
    enabled = settings.get("enabledPlugins", {}).get(PLUGIN_KEY) is True
    marketplace = settings.get("extraKnownMarketplaces", {}).get(MARKETPLACE_NAME)

    if enabled and marketplace:
        path = marketplace.get("source", {}).get("path", "unknown")
        current = str(PLUGIN_ROOT)
        note = "" if path == current else f"  (stale; current install is {current})"
        print(f"Registered: {PLUGIN_KEY}  path {path}{note}")
    elif enabled or marketplace:
        print(
            "STIsim plugin registration is inconsistent "
            f"(enabled={enabled}, marketplace={'present' if marketplace else 'missing'}). "
            "Run 'python -m stisim.ai install' to reconcile."
        )
    else:
        print("STIsim plugin is not registered.")
        print("Run 'python -m stisim.ai install' to register it.")


def main():
    parser = argparse.ArgumentParser(prog="python -m stisim.ai")
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("install", help="register the plugin with Claude Code").set_defaults(func=install)
    sub.add_parser("uninstall", help="remove the registration").set_defaults(func=uninstall)
    sub.add_parser("status", help="check current registration").set_defaults(func=status)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
