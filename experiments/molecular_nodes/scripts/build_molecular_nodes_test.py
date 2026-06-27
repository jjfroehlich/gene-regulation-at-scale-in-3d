#!/usr/bin/env python3
"""Build a small Molecular Nodes surface comparison scene."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import addon_utils
import bpy
from mathutils import Vector


ROOT = Path(os.environ.get("GENE_SCENE_ROOT", Path(__file__).resolve().parents[3])).resolve()
sys.path.insert(0, str(ROOT / "scripts"))
import build_gene_expression_scene as base  # noqa: E402

RCSB_DIR = ROOT / "assets" / "rcsb"
OUTPUT_DIR = ROOT / "experiments" / "molecular_nodes" / "outputs"
BLEND_PATH = OUTPUT_DIR / "molecular_nodes_test.blend"
PREVIEW_PATH = OUTPUT_DIR / "preview_molecular_nodes_test.png"
REPORT_PATH = OUTPUT_DIR / "molecular_nodes_test_report.json"
ANGSTROM_TO_MM = 0.04
BACKGROUND_COLOR = (0.985, 0.985, 0.965, 1.0)


def clean_scene() -> None:
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete()
    for collection in list(bpy.data.collections):
        bpy.data.collections.remove(collection)


def configure_scene() -> None:
    scene = bpy.context.scene
    scene.unit_settings.system = "METRIC"
    scene.unit_settings.length_unit = "MILLIMETERS"
    scene.unit_settings.scale_length = 0.001
    scene.render.engine = "BLENDER_WORKBENCH"
    scene.render.resolution_x = 1800
    scene.render.resolution_y = 950
    scene.display.shading.light = "STUDIO"
    scene.display.shading.color_type = "MATERIAL"
    scene.display.shading.background_type = "VIEWPORT"
    scene.display.shading.background_color = BACKGROUND_COLOR[:3]
    scene.view_settings.view_transform = "Standard"
    world = scene.world or bpy.data.worlds.new("World")
    scene.world = world
    world.color = BACKGROUND_COLOR[:3]
    world.use_nodes = True
    background = world.node_tree.nodes.get("Background")
    if background:
        background.inputs["Color"].default_value = BACKGROUND_COLOR
        background.inputs["Strength"].default_value = 1.0


def create_material(name: str, color: tuple[float, float, float, float]) -> bpy.types.Material:
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = color
    return mat


def create_text(name: str, body: str, location: tuple[float, float, float], size: float, material: bpy.types.Material) -> None:
    curve = bpy.data.curves.new(name, "FONT")
    curve.body = body
    curve.align_x = "CENTER"
    curve.align_y = "CENTER"
    curve.size = size
    curve.materials.append(material)
    obj = bpy.data.objects.new(name, curve)
    obj.location = location
    bpy.context.scene.collection.objects.link(obj)


def bounds(objects: list[bpy.types.Object]) -> tuple[float, float, float]:
    coords = []
    for obj in objects:
        if obj.type != "MESH":
            continue
        for corner in obj.bound_box:
            coords.append(obj.matrix_world @ Vector(corner))
    if not coords:
        return (0.0, 0.0, 0.0)
    return (
        max(c.x for c in coords) - min(c.x for c in coords),
        max(c.y for c in coords) - min(c.y for c in coords),
        max(c.z for c in coords) - min(c.z for c in coords),
    )


def expected_bbox_mm(pdb_id: str) -> tuple[float, float, float]:
    atoms = base.parse_atom_site(RCSB_DIR / f"{pdb_id}.cif")
    points = base.residue_points(atoms)
    vectors = [Vector(point["pos_A"]) for point in points]
    return (
        (max(v.x for v in vectors) - min(v.x for v in vectors)) * ANGSTROM_TO_MM,
        (max(v.y for v in vectors) - min(v.y for v in vectors)) * ANGSTROM_TO_MM,
        (max(v.z for v in vectors) - min(v.z for v in vectors)) * ANGSTROM_TO_MM,
    )


def import_with_molecular_nodes(pdb_id: str, target_x: float) -> dict:
    before = {obj.name for obj in bpy.data.objects}
    filepath = str((RCSB_DIR / f"{pdb_id}.cif").resolve())
    result = bpy.ops.mn.import_local(
        filepath=filepath,
        style="surface",
        node_setup=True,
        centre=True,
        centre_type="mass",
        remove_solvent=True,
        assembly=False,
    )
    created = [obj for obj in bpy.data.objects if obj.name not in before]
    raw_bbox = bounds(created)
    expected_bbox = expected_bbox_mm(pdb_id)
    raw_max = max(raw_bbox)
    expected_max = max(expected_bbox)
    scale_factor = expected_max / raw_max if raw_max > 0 else ANGSTROM_TO_MM
    for obj in created:
        obj.scale = (scale_factor, scale_factor, scale_factor)
        obj.location.x += target_x
        obj["angstrom_to_mm"] = ANGSTROM_TO_MM
        obj["molecular_nodes_empirical_scale_factor"] = scale_factor
        obj["source_pdb_id"] = pdb_id
        obj["style_source"] = "Molecular Nodes surface"
    bpy.context.view_layer.update()
    scaled_bbox = bounds(created)
    return {
        "pdb_id": pdb_id,
        "operator_result": list(result),
        "objects": [obj.name for obj in created],
        "object_count": len(created),
        "raw_molecular_nodes_bbox": raw_bbox,
        "expected_bbox_mm_from_cif": expected_bbox,
        "applied_scale_factor": scale_factor,
        "bbox_after_scale_mm": scaled_bbox,
    }


def add_camera() -> None:
    light_data = bpy.data.lights.new("softbox", "AREA")
    light_data.energy = 500.0
    light_data.size = 55.0
    light = bpy.data.objects.new("softbox", light_data)
    light.location = (0.0, -10.0, 60.0)
    bpy.context.scene.collection.objects.link(light)

    camera_data = bpy.data.cameras.new("Camera")
    camera = bpy.data.objects.new("Camera", camera_data)
    camera.location = (0.0, 0.0, 65.0)
    camera.rotation_euler = (0.0, 0.0, 0.0)
    camera_data.type = "ORTHO"
    camera_data.ortho_scale = 14.0
    bpy.context.scene.collection.objects.link(camera)
    bpy.context.scene.camera = camera


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    clean_scene()
    configure_scene()
    enabled = addon_utils.check("bl_ext.blender_org.molecularnodes")[1]
    report = {
        "title": "Molecular Nodes surface trial",
        "addon_module": "bl_ext.blender_org.molecularnodes",
        "addon_enabled": enabled,
        "angstrom_to_mm": ANGSTROM_TO_MM,
        "entries": [],
        "errors": [],
    }
    text_mat = create_material("text_grey", (0.22, 0.22, 0.22, 1.0))
    try:
        if not enabled:
            addon_utils.enable("bl_ext.blender_org.molecularnodes", default_set=False, persistent=False)
            report["addon_enabled_after_enable"] = addon_utils.check("bl_ext.blender_org.molecularnodes")[1]
        report["entries"].append(import_with_molecular_nodes("1NKP", -4.0))
        report["entries"].append(import_with_molecular_nodes("1EHZ", 4.0))
        create_text("label_1NKP", "Molecular Nodes surface: 1NKP", (-4.0, 5.0, 0.0), 0.38, text_mat)
        create_text("label_1EHZ", "Molecular Nodes surface: 1EHZ tRNA", (4.0, 5.0, 0.0), 0.38, text_mat)
        create_text("label_scale", "Empirically matched to 1 A = 0.04 mm bbox scale", (0.0, -5.7, 0.0), 0.32, text_mat)
    except Exception as exc:  # noqa: BLE001 - comparison trial should report blocker and still produce files.
        report["errors"].append(repr(exc))
        create_text("label_error", f"Molecular Nodes trial failed: {exc!r}", (0.0, 0.0, 0.0), 1.1, text_mat)
    add_camera()
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_PATH))
    bpy.context.scene.render.filepath = str(PREVIEW_PATH)
    bpy.ops.render.render(write_still=True)
    print(f"Wrote {BLEND_PATH}")
    print(f"Wrote {PREVIEW_PATH}")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
