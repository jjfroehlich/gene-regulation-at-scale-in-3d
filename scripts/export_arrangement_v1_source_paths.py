#!/usr/bin/env python3
"""Export edited arrangement V1 source curves for a subsequent rebuild."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import bpy


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import arrangement_v1_config as arrangement  # noqa: E402


BLEND_PATH = arrangement.OUTPUT_DIR / "gene_expression_arrangement_v1.blend"
OUTPUT_PATH = arrangement.EDITED_PATHS_PATH
CURVE_NAMES = ("DNA_source_path", "MRNA_source_path")


def spline_points_mm(obj: bpy.types.Object) -> list[list[float]]:
    points: list[list[float]] = []
    for spline in obj.data.splines:
        if spline.type == "POLY":
            for point in spline.points:
                world = obj.matrix_world @ point.co.to_3d()
                points.append([world.x, world.y, world.z])
        elif spline.type == "BEZIER":
            for point in spline.bezier_points:
                world = obj.matrix_world @ point.co
                points.append([world.x, world.y, world.z])
    return points


def main() -> None:
    if not BLEND_PATH.exists():
        raise FileNotFoundError(f"Missing arrangement blend: {BLEND_PATH}")
    bpy.ops.wm.open_mainfile(filepath=str(BLEND_PATH))
    data = {
        "source_blend": str(BLEND_PATH),
        "note": "Generated from editable arrangement source curves. The next arrangement workflow run will use these points if present.",
        "curves": {},
    }
    for name in CURVE_NAMES:
        obj = bpy.data.objects.get(name)
        if obj is None or obj.type != "CURVE":
            raise RuntimeError(f"Missing source curve {name}")
        points = spline_points_mm(obj)
        if len(points) < 2:
            raise RuntimeError(f"Source curve {name} has fewer than two points")
        data[name] = {"points_mm": points, "point_count": len(points)}
        data["curves"][name] = {"point_count": len(points)}
    nucleosome = bpy.data.objects.get("ATTACH_DNA_Nucleosome")
    if nucleosome is not None:
        data["nucleosome_loop"] = {
            "center_mm": [nucleosome.location.x, nucleosome.location.y, nucleosome.location.z],
            "roll_deg": math.degrees(nucleosome.rotation_euler.z),
        }
    OUTPUT_PATH.write_text(json.dumps(data, indent=2), encoding="utf-8")
    print(f"Wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
