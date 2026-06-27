#!/usr/bin/env python3
"""Compare direct procedural nucleic acids, PyMOL calibrators, and Molecular Nodes trials."""

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
import blender_nucleic_meshes as nucleic_meshes  # noqa: E402


EXPERIMENT_DIR = ROOT / "experiments" / "procedural_nucleic_acids"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
ASSET_DIR = EXPERIMENT_DIR / "assets"
PYMOL_CALIBRATOR_DIR = ASSET_DIR / "pymol_calibrator_surfaces"
OXDNA_DIR = ASSET_DIR / "molecular_nodes_oxdna"
RCSB_DIR = ROOT / "assets" / "rcsb"
MANIFEST_PATH = ROOT / "config" / "scene_manifest.json"
CALIBRATOR_REPORT = OUTPUT_DIR / "nucleic_acid_calibrator_analysis.json"
BLEND_PATH = OUTPUT_DIR / "nucleic_acid_generation_comparison.blend"
PREVIEW_PATH = OUTPUT_DIR / "preview_nucleic_acid_generation_comparison.png"
REPORT_PATH = OUTPUT_DIR / "nucleic_acid_generation_comparison_report.json"
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
    scene.render.resolution_x = 2400
    scene.render.resolution_y = 1500
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
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf:
        bsdf.inputs["Base Color"].default_value = color
        bsdf.inputs["Roughness"].default_value = 0.9
    return mat


def add_text(name: str, body: str, location: tuple[float, float, float], size: float, mat: bpy.types.Material) -> None:
    curve = bpy.data.curves.new(name, "FONT")
    curve.body = body
    curve.align_x = "CENTER"
    curve.align_y = "CENTER"
    curve.size = size
    curve.extrude = 0.02
    curve.materials.append(mat)
    obj = bpy.data.objects.new(name, curve)
    obj.location = location
    bpy.context.scene.collection.objects.link(obj)


def object_bounds(objects: list[bpy.types.Object]) -> tuple[float, float, float]:
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


def translate_objects(objects: list[bpy.types.Object], target_center: Vector) -> None:
    coords = []
    for obj in objects:
        if obj.type != "MESH":
            continue
        for corner in obj.bound_box:
            coords.append(obj.matrix_world @ Vector(corner))
    if not coords:
        return
    center = Vector(
        (
            (min(c.x for c in coords) + max(c.x for c in coords)) * 0.5,
            (min(c.y for c in coords) + max(c.y for c in coords)) * 0.5,
            (min(c.z for c in coords) + max(c.z for c in coords)) * 0.5,
        )
    )
    delta = target_center - center
    for obj in objects:
        obj.location += delta


def build_custom_dna(materials: dict, origin: Vector, bp: int = 80) -> dict:
    bp_to_mm = 0.136
    pitch_bp = 10.5
    strand_center_radius = 0.304
    strand_radius = 0.092
    base_radius = 0.034
    strand_a = []
    strand_b = []
    for i in range(bp + 1):
        x = i * bp_to_mm
        theta = 2.0 * math.pi * i / pitch_bp
        center = origin + Vector((x, 0.0, 0.0))
        radial = Vector((0.0, math.sin(theta), math.cos(theta)))
        strand_a.append(center + radial * strand_center_radius)
        strand_b.append(center - radial * strand_center_radius)
    objects = []
    for name, points, mat in [
        ("custom_DNA_strand_A", strand_a, materials["dna_a"]),
        ("custom_DNA_strand_B", strand_b, materials["dna_b"]),
    ]:
        vertices, faces = nucleic_meshes.tube_mesh(points, strand_radius, 9, 0.08)
        obj, report = nucleic_meshes.create_mesh_object(name, vertices, faces, mat, bpy.context.scene.collection, "procedural_blender_surface_proxy")
        objects.append(obj)
    bridges = [(strand_a[i], strand_b[i]) for i in range(0, len(strand_a), 2)]
    vertices, faces = nucleic_meshes.cylinder_segments_mesh(bridges, base_radius, 6, 0.03)
    obj, bridge_report = nucleic_meshes.create_mesh_object(
        "custom_DNA_base_pairs",
        vertices,
        faces,
        materials["dna_a"],
        bpy.context.scene.collection,
        "procedural_blender_surface_proxy",
    )
    objects.append(obj)
    return {
        "type": "custom_direct_blender_dna",
        "bp": bp,
        "axis_length_mm": bp * bp_to_mm,
        "bp_spacing_mm": bp_to_mm,
        "pitch_bp": pitch_bp,
        "target_envelope_diameter_mm": 0.8,
        "estimated_envelope_diameter_mm": 2.0 * (strand_center_radius + strand_radius),
        "objects": [obj.name for obj in objects],
        "bbox_mm": object_bounds(objects),
    }


def build_custom_rna(materials: dict, origin: Vector, nt: int = 160) -> dict:
    nt_to_mm = 0.12
    tube_radius = 0.22
    target_contour = nt * nt_to_mm

    def make_path(span: float) -> list[Vector]:
        points = []
        for i in range(nt + 1):
            t = i / nt
            x = span * t
            y = 0.75 * math.sin(2.0 * math.pi * i / 34.0)
            points.append(origin + Vector((x, y, 0.0)))
        return points

    def length(points: list[Vector]) -> float:
        return sum((points[i] - points[i - 1]).length for i in range(1, len(points)))

    low = 0.0
    high = target_contour
    while length(make_path(high)) < target_contour:
        high *= 1.5
    for _ in range(50):
        mid = (low + high) * 0.5
        if length(make_path(mid)) < target_contour:
            low = mid
        else:
            high = mid
    path = make_path(high)
    objects = []
    vertices, faces = nucleic_meshes.tube_mesh(path, tube_radius, 9, 0.12)
    obj, tube_report = nucleic_meshes.create_mesh_object(
        "custom_RNA_direct_surface",
        vertices,
        faces,
        materials["rna"],
        bpy.context.scene.collection,
        "procedural_blender_surface_proxy",
    )
    objects.append(obj)
    lobes = []
    for i in range(0, nt, 8):
        center = path[i]
        tangent = (path[min(nt, i + 1)] - path[max(0, i - 1)]).normalized()
        normal, _ = nucleic_meshes.path_frame(tangent)
        lobes.append((center + normal * (0.16 if (i // 8) % 2 == 0 else -0.16), center + normal * (0.162 if (i // 8) % 2 == 0 else -0.162)))
    vertices, faces = nucleic_meshes.cylinder_segments_mesh(lobes, 0.11, 7, 0.05)
    obj, lobe_report = nucleic_meshes.create_mesh_object(
        "custom_RNA_base_lobes",
        vertices,
        faces,
        materials["rna"],
        bpy.context.scene.collection,
        "procedural_blender_surface_proxy",
    )
    objects.append(obj)
    return {
        "type": "custom_direct_blender_rna",
        "nt": nt,
        "target_contour_mm": target_contour,
        "measured_polyline_mm": length(path),
        "nt_spacing_mm": nt_to_mm,
        "tube_radius_mm": tube_radius,
        "objects": [obj.name for obj in objects],
        "bbox_mm": object_bounds(objects),
    }


def parse_obj(path: Path) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]]]:
    vertices = []
    faces = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            if raw_line.startswith("v "):
                _, x, y, z = raw_line.split()[:4]
                vertices.append((float(x) * ANGSTROM_TO_MM, float(y) * ANGSTROM_TO_MM, float(z) * ANGSTROM_TO_MM))
            elif raw_line.startswith("f "):
                face = []
                for token in raw_line.split()[1:]:
                    index = token.split("/")[0]
                    if index:
                        face.append(int(index) - 1)
                if len(face) >= 3:
                    faces.append(tuple(face))
    return vertices, faces


def import_pymol_calibrator(pdb_id: str, target_center: Vector, mat: bpy.types.Material) -> dict:
    path = PYMOL_CALIBRATOR_DIR / pdb_id / f"{pdb_id}_nucleic_surface.obj"
    if not path.exists() or path.stat().st_size == 0:
        return {"pdb_id": pdb_id, "status": "missing_pymol_calibrator_surface", "path": str(path)}
    vertices, faces = parse_obj(path)
    obj, report = nucleic_meshes.create_mesh_object(
        f"PyMOL calibrator {pdb_id} nucleic surface",
        vertices,
        faces,
        mat,
        bpy.context.scene.collection,
        "pymol_calibrator_nucleic_surface",
    )
    translate_objects([obj], target_center)
    report.update({"pdb_id": pdb_id, "status": "imported", "path": str(path.relative_to(ROOT)).replace("\\", "/")})
    return report


def best_rna_calibrator() -> str | None:
    if not CALIBRATOR_REPORT.exists():
        return None
    report = json.loads(CALIBRATOR_REPORT.read_text(encoding="utf-8"))
    candidates = [
        entry
        for entry in report.get("entries", [])
        if entry.get("expected_type") == "rna" and entry.get("status") == "ok" and entry.get("nucleic_residue_count", 0) > 0
    ]
    if not candidates:
        return None
    candidates.sort(key=lambda item: item.get("nucleic_residue_count", 10**9))
    return candidates[0]["pdb_id"]


def import_mn_surface(pdb_id: str, target_center: Vector) -> dict:
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
    translate_objects(created, target_center)
    for obj in created:
        obj["style_source"] = "Molecular Nodes Style Surface"
    return {
        "pdb_id": pdb_id,
        "operator_result": list(result),
        "object_count": len(created),
        "bbox_raw_blender_units": object_bounds(created),
        "canonical_suitability": "comparison_only_real_structure_import_not_path_generation",
    }


def write_oxdna_files(name: str, polymer_type: str, count: int, duplex: bool, rise_A: float) -> tuple[Path, Path, dict]:
    OXDNA_DIR.mkdir(parents=True, exist_ok=True)
    top_path = OXDNA_DIR / f"{name}.top"
    dat_path = OXDNA_DIR / f"{name}.dat"
    bases = ("ACGU" if polymer_type == "RNA" else "ACGT") * (count // 4 + 1)
    bases_a = bases[:count]
    if duplex:
        bases_b = bases_a.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        top_text = f"{count * 2} 2 5->3\n{bases_a} type={polymer_type}\n{bases_b} type={polymer_type}\n"
    else:
        top_text = f"{count} 1 5->3\n{bases_a} type={polymer_type}\n"
    top_path.write_text(top_text, encoding="utf-8")
    lines = ["t = 0\n", "b = 300 300 300\n", "E = 0 0 0\n"]
    phases = (0.0, math.pi) if duplex else (0.0,)
    radius_A = 8.5 if duplex else 0.0
    for i in range(count):
        theta = 2.0 * math.pi * i / 10.5
        for phase in phases:
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
    return top_path, dat_path, {
        "count": count,
        "polymer_type": polymer_type,
        "duplex": duplex,
        "expected_axis_length_mm": (count - 1) * rise_A * ANGSTROM_TO_MM,
    }


def import_mn_oxdna(name: str, polymer_type: str, count: int, duplex: bool, rise_A: float, target_center: Vector) -> dict:
    top_path, dat_path, geometry = write_oxdna_files(name, polymer_type, count, duplex, rise_A)
    trajectory = importlib.import_module("bl_ext.blender_org.molecularnodes.entities.trajectory")
    before = {obj.name for obj in bpy.data.objects}
    result = trajectory.load_oxdna(top=top_path, traj=dat_path, name=name, style="oxdna", world_scale=0.004)
    created = [obj for obj in bpy.data.objects if obj.name not in before]
    translate_objects(created, target_center)
    bbox = object_bounds(created)
    return {
        "topology": str(top_path.relative_to(ROOT)).replace("\\", "/"),
        "trajectory": str(dat_path.relative_to(ROOT)).replace("\\", "/"),
        "object_count": len(created),
        "created_objects": [obj.name for obj in created],
        "bbox_mm": bbox,
        "geometry": geometry,
        "axis_length_error_mm": bbox[0] - geometry["expected_axis_length_mm"] if bbox else None,
        "result_type": type(result).__name__,
        "canonical_suitability": "not_surface_style_coarse_grained_oxdna",
    }


def disable_mn_save_handlers(report: dict) -> list:
    removed = []
    for handler in list(bpy.app.handlers.save_post):
        module = getattr(handler, "__module__", "")
        if "molecularnodes" in module.lower():
            bpy.app.handlers.save_post.remove(handler)
            removed.append(handler)
    if removed:
        report["save_notes"] = {"molecular_nodes_save_handlers_temporarily_removed": len(removed)}
    return removed


def restore_mn_save_handlers(handlers: list) -> None:
    for handler in handlers:
        if handler not in bpy.app.handlers.save_post:
            bpy.app.handlers.save_post.append(handler)


def add_camera() -> None:
    light_data = bpy.data.lights.new("softbox", "AREA")
    light_data.energy = 850.0
    light_data.size = 50.0
    light = bpy.data.objects.new("softbox", light_data)
    light.location = (0.0, -10.0, 55.0)
    bpy.context.scene.collection.objects.link(light)
    camera_data = bpy.data.cameras.new("Camera")
    camera = bpy.data.objects.new("Camera", camera_data)
    camera.location = (0.0, 0.0, 62.0)
    camera.rotation_euler = (0.0, 0.0, 0.0)
    camera_data.type = "ORTHO"
    camera_data.ortho_scale = 33.0
    bpy.context.scene.collection.objects.link(camera)
    bpy.context.scene.camera = camera


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    clean_scene()
    configure_scene()
    materials = {
        "dna_a": material("dna_a", (0.88, 0.48, 0.24, 1.0)),
        "dna_b": material("dna_b", (0.55, 0.28, 0.13, 1.0)),
        "rna": material("rna", (0.62, 0.58, 0.28, 1.0)),
        "pymol_dna": material("pymol_dna", (0.80, 0.58, 0.42, 1.0)),
        "pymol_rna": material("pymol_rna", (0.66, 0.62, 0.42, 1.0)),
        "text": material("text_grey", (0.20, 0.20, 0.20, 1.0)),
    }
    report = {
        "title": "Procedural nucleic-acid generation comparison",
        "entries": {},
        "errors": [],
    }

    report["entries"]["custom_direct_dna"] = build_custom_dna(materials, Vector((-17.0, 7.0, 0.0)))
    report["entries"]["custom_direct_rna"] = build_custom_rna(materials, Vector((-17.0, -5.0, 0.0)))
    add_text("label_custom", "custom direct Blender", (-12.0, 13.3, 1.2), 0.45, materials["text"])

    rna_calibrator = best_rna_calibrator()
    report["entries"]["pymol_calibrator_dna_1BNA"] = import_pymol_calibrator("1BNA", Vector((0.0, 7.0, 0.0)), materials["pymol_dna"])
    if rna_calibrator:
        report["entries"][f"pymol_calibrator_rna_{rna_calibrator}"] = import_pymol_calibrator(rna_calibrator, Vector((0.0, -5.0, 0.0)), materials["pymol_rna"])
    add_text("label_pymol", "PyMOL real-structure calibrators", (0.0, 13.3, 1.2), 0.45, materials["text"])

    try:
        enabled = addon_utils.check("bl_ext.blender_org.molecularnodes")[1]
        if not enabled:
            addon_utils.enable("bl_ext.blender_org.molecularnodes", default_set=False, persistent=False)
        report["molecular_nodes_enabled"] = addon_utils.check("bl_ext.blender_org.molecularnodes")[1]
        report["entries"]["mn_style_surface_dna_1BNA"] = import_mn_surface("1BNA", Vector((10.0, 7.0, 0.0)))
        if rna_calibrator:
            report["entries"][f"mn_style_surface_rna_{rna_calibrator}"] = import_mn_surface(rna_calibrator, Vector((10.0, -5.0, 0.0)))
        report["entries"]["mn_oxdna_generated_dna"] = import_mn_oxdna("MN generated oxDNA 42bp", "DNA", 42, True, 3.4, Vector((20.0, 7.0, 0.0)))
        report["entries"]["mn_oxrna_generated_rna"] = import_mn_oxdna("MN generated oxRNA 80nt", "RNA", 80, False, 3.0, Vector((20.0, -5.0, 0.0)))
    except Exception as exc:  # noqa: BLE001
        report["errors"].append(repr(exc))
    add_text("label_mn", "Molecular Nodes", (15.0, 13.3, 1.2), 0.45, materials["text"])

    if CALIBRATOR_REPORT.exists():
        report["calibrator_analysis"] = str(CALIBRATOR_REPORT.relative_to(ROOT)).replace("\\", "/")
    add_camera()
    removed_handlers = disable_mn_save_handlers(report)
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    try:
        bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_PATH))
    finally:
        restore_mn_save_handlers(removed_handlers)
    bpy.context.scene.render.filepath = str(PREVIEW_PATH)
    bpy.ops.render.render(write_still=True)
    print(f"Wrote {BLEND_PATH}")
    print(f"Wrote {PREVIEW_PATH}")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
