#!/usr/bin/env python3
"""Build a Molecular Nodes DNA/RNA comparison scene."""

from __future__ import annotations

import importlib
import json
import math
import os
import sys
from pathlib import Path

import addon_utils
import bpy
from mathutils import Vector


ROOT = Path(os.environ.get("GENE_SCENE_ROOT", Path(__file__).resolve().parents[3])).resolve()
sys.path.insert(0, str(ROOT / "scripts"))
import build_gene_expression_scene as base  # noqa: E402


EXPERIMENT_DIR = ROOT / "experiments" / "molecular_nodes"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ASSET_DIR = EXPERIMENT_DIR / "assets"
OXDNA_DIR = ASSET_DIR / "oxdna"
BLEND_PATH = OUTPUT_DIR / "molecular_nodes_dna_rna_comparison.blend"
PREVIEW_PATH = OUTPUT_DIR / "preview_molecular_nodes_dna_rna_comparison.png"
REPORT_PATH = OUTPUT_DIR / "molecular_nodes_dna_rna_comparison_report.json"
RCSB_DIR = ROOT / "assets" / "rcsb"
REDUCED_SURFACE_DIR = ROOT / "assets" / "pymol_exports" / "surface_assets_reduced"
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
    scene.render.resolution_x = 2200
    scene.render.resolution_y = 1250
    scene.display.shading.light = "STUDIO"
    scene.display.shading.color_type = "MATERIAL"
    scene.display.shading.background_type = "VIEWPORT"
    scene.display.shading.background_color = BACKGROUND_COLOR[:3]
    scene.view_settings.view_transform = "Standard"
    world = scene.world or bpy.data.worlds.new("World")
    scene.world = world
    world.color = BACKGROUND_COLOR[:3]


def material(name: str, color: tuple[float, float, float, float]) -> bpy.types.Material:
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = color
    return mat


def add_text(name: str, body: str, location: tuple[float, float, float], size: float, mat: bpy.types.Material) -> None:
    curve = bpy.data.curves.new(name, "FONT")
    curve.body = body
    curve.align_x = "CENTER"
    curve.align_y = "CENTER"
    curve.size = size
    curve.extrude = 0.03
    curve.materials.append(mat)
    obj = bpy.data.objects.new(name, curve)
    obj.location = (location[0], location[1], location[2] + 1.5)
    bpy.context.scene.collection.objects.link(obj)
    bpy.ops.object.select_all(action="DESELECT")
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.convert(target="MESH")
    bpy.ops.object.mode_set(mode="EDIT")
    bpy.ops.mesh.select_all(action="SELECT")
    bpy.ops.mesh.flip_normals()
    bpy.ops.object.mode_set(mode="OBJECT")


def parse_obj(path: Path) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]]]:
    vertices = []
    faces = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            if raw_line.startswith("v "):
                _, x, y, z = raw_line.split()[:4]
                vertices.append((float(x), float(y), float(z)))
            elif raw_line.startswith("f "):
                face = []
                for token in raw_line.split()[1:]:
                    index = token.split("/")[0]
                    if index:
                        face.append(int(index) - 1)
                if len(face) >= 3:
                    faces.append(tuple(face))
    return vertices, faces


def import_reduced_proxy_component(path: Path, name: str, mat: bpy.types.Material, offset: Vector, display_scale: float) -> dict:
    vertices_A, faces = parse_obj(path)
    vertices = []
    for x, y, z in vertices_A:
        world = Vector((x, y, z)) * ANGSTROM_TO_MM * display_scale + offset
        vertices.append((world.x, world.y, world.z))
    mesh = bpy.data.meshes.new(f"{name}_mesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    for polygon in mesh.polygons:
        polygon.use_smooth = True
    mesh.materials.append(mat)
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.scene.collection.objects.link(obj)
    obj["style_source"] = "Canonical procedural PyMOL surface proxy"
    return {"object": obj.name, "source_obj": str(path.relative_to(ROOT)).replace("\\", "/"), "faces": len(faces)}


def add_our_proxy_panel(text_mat: bpy.types.Material) -> dict:
    mats = {
        "dna_a": material("proxy_dna_a", (0.87, 0.48, 0.24, 1.0)),
        "dna_b": material("proxy_dna_b", (0.55, 0.29, 0.15, 1.0)),
        "mrna": material("proxy_mrna", (0.62, 0.58, 0.28, 1.0)),
    }
    display_scale = 0.045
    reports = []
    offset = Vector((-13.5, -2.0, 0.0))
    for component, mat_key in [("strand_A", "dna_a"), ("strand_B", "dna_b"), ("base_pairs", "dna_a")]:
        path = REDUCED_SURFACE_DIR / "DNA_PROXY" / f"DNA_PROXY_surface_{component}.obj"
        reports.append(import_reduced_proxy_component(path, f"Our DNA proxy {component}", mats[mat_key], offset, display_scale))
    for component in ["utr5", "coding", "utr3"]:
        path = REDUCED_SURFACE_DIR / "MRNA_PROXY" / f"MRNA_PROXY_surface_{component}.obj"
        reports.append(import_reduced_proxy_component(path, f"Our mRNA proxy {component}", mats["mrna"], offset + Vector((0.0, 4.6, 0.0)), display_scale))
    add_text("label_our_proxy", "Our repaired PyMOL surface proxy", (-13.5, 8.6, 0.0), 0.42, text_mat)
    add_text("label_our_proxy_scale", "path-controlled geometry, displayed at 4.5%", (-13.5, -8.6, 0.0), 0.28, text_mat)
    return {"display_scale": display_scale, "components": reports}


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


def import_mn_surface(pdb_id: str, target_x: float) -> dict:
    before = {obj.name for obj in bpy.data.objects}
    result = bpy.ops.mn.import_local(
        filepath=str((RCSB_DIR / f"{pdb_id}.cif").resolve()),
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
        obj.location.y += 0.5
        obj["style_source"] = "Molecular Nodes Style Surface"
        obj["molecular_nodes_empirical_scale_factor"] = scale_factor
    bpy.context.view_layer.update()
    return {
        "pdb_id": pdb_id,
        "operator_result": list(result),
        "object_count": len(created),
        "raw_bbox": raw_bbox,
        "expected_bbox_mm": expected_bbox,
        "applied_empirical_scale_factor": scale_factor,
        "bbox_after_scale_mm": bounds(created),
        "canonical_suitability": "not_canonical_scale_requires_empirical_bbox_factor",
    }


def write_oxdna_files() -> tuple[Path, Path, dict]:
    OXDNA_DIR.mkdir(parents=True, exist_ok=True)
    top_path = OXDNA_DIR / "generated_42bp_dna.top"
    dat_path = OXDNA_DIR / "generated_42bp_dna.dat"
    bases_a = "ACGT" * 11
    bases_a = bases_a[:42]
    bases_b = bases_a.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    top_path.write_text(f"84 2 5->3\n{bases_a} type=DNA\n{bases_b} type=DNA\n", encoding="utf-8")
    lines = ["t = 0\n", "b = 300 300 300\n", "E = 0 0 0\n"]
    rise_A = 3.4
    radius_A = 8.5
    for i in range(42):
        theta = 2.0 * math.pi * i / 10.5
        for phase in (0.0, math.pi):
            x = i * rise_A
            y = radius_A * math.cos(theta + phase)
            z = radius_A * math.sin(theta + phase)
            tangent = (1.0, 0.0, 0.0)
            normal = (0.0, math.cos(theta + phase), math.sin(theta + phase))
            lines.append(
                f"{x:.3f} {y:.3f} {z:.3f} "
                f"{tangent[0]:.3f} {tangent[1]:.3f} {tangent[2]:.3f} "
                f"{normal[0]:.3f} {normal[1]:.3f} {normal[2]:.3f} "
                "0 0 0 0 0 0\n"
            )
    dat_path.write_text("".join(lines), encoding="utf-8")
    return top_path, dat_path, {"bp": 42, "expected_axis_length_mm": 42 * rise_A * ANGSTROM_TO_MM}


def import_oxdna_trial(target_x: float) -> dict:
    top_path, dat_path, geometry = write_oxdna_files()
    mn = importlib.import_module("bl_ext.blender_org.molecularnodes")
    trajectory = mn.entities.trajectory
    before = {obj.name for obj in bpy.data.objects}
    result = trajectory.load_oxdna(top=top_path, traj=dat_path, name="MN generated oxDNA 42bp", style="oxdna", world_scale=0.004)
    created = [obj for obj in bpy.data.objects if obj.name not in before]
    for obj in created:
        obj.location.x += target_x
        obj.location.y -= 1.0
        obj["style_source"] = "Molecular Nodes oxDNA"
        obj["world_scale_without_empirical_bbox_fit"] = 0.004
    bpy.context.view_layer.update()
    return {
        "topology": str(top_path.relative_to(ROOT)).replace("\\", "/"),
        "trajectory": str(dat_path.relative_to(ROOT)).replace("\\", "/"),
        "object_count": len(created),
        "created_objects": [obj.name for obj in created],
        "bbox_mm": bounds(created),
        "geometry": geometry,
        "result_type": type(result).__name__,
        "canonical_suitability": "not_surface_style_oxdna_ribbon_coarse_grained",
    }


def add_camera() -> None:
    light_data = bpy.data.lights.new("softbox", "AREA")
    light_data.energy = 800.0
    light_data.size = 60.0
    light = bpy.data.objects.new("softbox", light_data)
    light.location = (0.0, -5.0, 65.0)
    bpy.context.scene.collection.objects.link(light)

    camera_data = bpy.data.cameras.new("Camera")
    camera = bpy.data.objects.new("Camera", camera_data)
    camera.location = (0.0, 0.0, 70.0)
    camera.rotation_euler = (0.0, 0.0, 0.0)
    camera_data.type = "ORTHO"
    camera_data.ortho_scale = 24.0
    bpy.context.scene.collection.objects.link(camera)
    bpy.context.scene.camera = camera


def disable_molecular_nodes_save_handlers(report: dict) -> list:
    """Keep the comparison blend save clean when oxDNA registers a DAT trajectory."""
    removed = []
    for handler in list(bpy.app.handlers.save_post):
        module = getattr(handler, "__module__", "")
        name = getattr(handler, "__name__", repr(handler))
        if "molecularnodes" in module.lower():
            bpy.app.handlers.save_post.remove(handler)
            removed.append(handler)
    if removed:
        report["save_notes"] = {
            "molecular_nodes_save_handlers_temporarily_removed": [
                f"{getattr(handler, '__module__', '')}.{getattr(handler, '__name__', repr(handler))}"
                for handler in removed
            ],
            "reason": (
                "Molecular Nodes oxDNA uses a .dat trajectory that the add-on save hook "
                "tries to reopen through MDAnalysis; the comparison stores the evaluated "
                "Blender objects and keeps Molecular Nodes out of the canonical workflow."
            ),
        }
    return removed


def restore_molecular_nodes_save_handlers(handlers: list) -> None:
    for handler in handlers:
        if handler not in bpy.app.handlers.save_post:
            bpy.app.handlers.save_post.append(handler)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    clean_scene()
    configure_scene()
    text_mat = material("text_grey", (0.22, 0.22, 0.22, 1.0))
    report = {
        "title": "Molecular Nodes DNA/RNA full comparison",
        "addon_module": "bl_ext.blender_org.molecularnodes",
        "entries": {},
        "errors": [],
    }
    try:
        enabled = addon_utils.check("bl_ext.blender_org.molecularnodes")[1]
        if not enabled:
            addon_utils.enable("bl_ext.blender_org.molecularnodes", default_set=False, persistent=False)
        report["addon_enabled"] = addon_utils.check("bl_ext.blender_org.molecularnodes")[1]
        report["entries"]["our_pymol_proxy"] = add_our_proxy_panel(text_mat)
        report["entries"]["mn_style_surface_1NKP"] = import_mn_surface("1NKP", -1.5)
        report["entries"]["mn_style_surface_1EHZ"] = import_mn_surface("1EHZ", 6.0)
        report["entries"]["mn_oxdna_generated"] = import_oxdna_trial(13.5)
        add_text("label_mn_surface", "Molecular Nodes Style Surface", (2.5, 8.0, 0.0), 0.42, text_mat)
        add_text("label_mn_surface_note", "real structures; empirical scale fit", (2.5, -8.0, 0.0), 0.28, text_mat)
        add_text("label_oxdna", "Molecular Nodes oxDNA", (13.5, 8.0, 0.0), 0.42, text_mat)
        add_text("label_oxdna_note", "path-controllable, coarse ribbon, not surface", (13.5, -8.0, 0.0), 0.28, text_mat)
    except Exception as exc:  # noqa: BLE001
        report["errors"].append(repr(exc))
        add_text("label_error", f"Molecular Nodes comparison failed: {exc!r}", (0.0, 0.0, 0.0), 0.48, text_mat)
    add_camera()
    removed_save_handlers = disable_molecular_nodes_save_handlers(report)
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    try:
        bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_PATH))
    finally:
        restore_molecular_nodes_save_handlers(removed_save_handlers)
    bpy.context.scene.render.filepath = str(PREVIEW_PATH)
    bpy.ops.render.render(write_still=True)
    print(f"Wrote {BLEND_PATH}")
    print(f"Wrote {PREVIEW_PATH}")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
