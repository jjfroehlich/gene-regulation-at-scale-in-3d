#!/usr/bin/env python3
"""Install and enable Blender MCP add-on in Blender user preferences."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import bpy


ROOT = Path(os.environ.get("GENE_SCENE_ROOT", Path(__file__).resolve().parents[1]))
SOURCE = ROOT / "assets" / "blender_mcp" / "addon.py"
STAGED = ROOT / "assets" / "blender_mcp" / "blender_mcp.py"


def main() -> None:
    if not SOURCE.exists():
        raise FileNotFoundError(SOURCE)
    shutil.copyfile(SOURCE, STAGED)
    bpy.ops.preferences.addon_install(filepath=str(STAGED), overwrite=True, target="DEFAULT")
    bpy.ops.preferences.addon_enable(module="blender_mcp")
    addon = bpy.context.preferences.addons.get("blender_mcp")
    if addon and hasattr(addon.preferences, "telemetry_consent"):
        addon.preferences.telemetry_consent = False
    scene = bpy.context.scene
    if hasattr(scene, "blendermcp_port"):
        scene.blendermcp_port = 9877
    for name in (
        "blendermcp_use_polyhaven",
        "blendermcp_use_hyper3d",
        "blendermcp_use_sketchfab",
        "blendermcp_use_hunyuan3d",
    ):
        if hasattr(scene, name):
            setattr(scene, name, False)
    bpy.ops.wm.save_userpref()
    print("Installed and enabled Blender MCP add-on as module 'blender_mcp'.")
    print("Telemetry consent in the Blender add-on preferences is set to false.")


if __name__ == "__main__":
    main()
