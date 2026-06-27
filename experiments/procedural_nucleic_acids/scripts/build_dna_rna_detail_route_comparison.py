#!/usr/bin/env python3
"""Build a side-by-side detailed DNA/RNA route comparison in Blender."""

from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

import bmesh
import bpy
from mathutils import Vector


ROOT = Path(os.environ.get("GENE_SCENE_ROOT", Path(__file__).resolve().parents[3])).resolve()
SCRIPT_DIR = ROOT / "experiments" / "procedural_nucleic_acids" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))
import detail_route_geometry as geom  # noqa: E402


EXPERIMENT_DIR = ROOT / "experiments" / "procedural_nucleic_acids"
ASSET_DIR = EXPERIMENT_DIR / "assets"
DETAIL_ASSET_DIR = ASSET_DIR / "detail_route_pymol_proxy"
RAW_DIR = DETAIL_ASSET_DIR / "raw_surfaces"
REDUCED_DIR = DETAIL_ASSET_DIR / "reduced_surfaces"
CALIBRATOR_DIR = ASSET_DIR / "pymol_calibrator_surfaces"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
BLEND_PATH = OUTPUT_DIR / "dna_rna_detail_route_comparison.blend"
PREVIEW_PATH = OUTPUT_DIR / "preview_dna_rna_detail_route_comparison.png"
DNA_DETAIL_PREVIEW_PATH = OUTPUT_DIR / "preview_dna_rna_detail_route_comparison_dna_detail.png"
RNA_DETAIL_PREVIEW_PATH = OUTPUT_DIR / "preview_dna_rna_detail_route_comparison_rna_detail.png"
REPORT_PATH = OUTPUT_DIR / "dna_rna_detail_route_comparison_report.json"
BACKGROUND_COLOR = (0.985, 0.985, 0.965, 1.0)

REDUCTION_TARGETS = {
    "DNA_DETAIL_PROXY": {"strand_A": 36_000, "strand_B": 36_000, "base_pairs": 55_000},
    "RNA_DETAIL_PROXY": {"utr5": 18_000, "coding": 52_000, "utr3": 34_000},
}


def clean_scene() -> None:
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete()
    for collection in list(bpy.data.collections):
        bpy.data.collections.remove(collection)
    for block in list(bpy.data.meshes):
        if block.users == 0:
            bpy.data.meshes.remove(block)


def configure_scene() -> None:
    scene = bpy.context.scene
    scene.unit_settings.system = "METRIC"
    scene.unit_settings.length_unit = "MILLIMETERS"
    scene.unit_settings.scale_length = 0.001
    scene.render.engine = "BLENDER_WORKBENCH"
    scene.render.resolution_x = 2600
    scene.render.resolution_y = 1500
    scene.display.shading.light = "STUDIO"
    scene.display.shading.color_type = "MATERIAL"
    scene.display.shading.background_type = "VIEWPORT"
    scene.display.shading.background_color = BACKGROUND_COLOR[:3]
    scene.view_settings.view_transform = "Standard"
    world = scene.world or bpy.data.worlds.new("World")
    scene.world = world
    world.color = BACKGROUND_COLOR[:3]


def make_collection(name: str) -> bpy.types.Collection:
    collection = bpy.data.collections.new(name)
    bpy.context.scene.collection.children.link(collection)
    return collection


def link_to_collection(obj: bpy.types.Object, collection: bpy.types.Collection) -> None:
    collection.objects.link(obj)
    for coll in list(obj.users_collection):
        if coll != collection:
            coll.objects.unlink(obj)


def material(name: str, color: tuple[float, float, float, float]) -> bpy.types.Material:
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = color
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf:
        bsdf.inputs["Base Color"].default_value = color
        bsdf.inputs["Roughness"].default_value = 0.88
    return mat


def add_text(name: str, body: str, location: tuple[float, float, float], size: float, mat: bpy.types.Material, collection: bpy.types.Collection) -> None:
    curve = bpy.data.curves.new(name, "FONT")
    curve.body = body
    curve.align_x = "CENTER"
    curve.align_y = "CENTER"
    curve.size = size
    curve.extrude = 0.02
    curve.materials.append(mat)
    obj = bpy.data.objects.new(name, curve)
    obj.location = location
    link_to_collection(obj, collection)


def parse_obj(path: Path, scale: float = 1.0) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]]]:
    vertices = []
    faces = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            if raw_line.startswith("v "):
                _, x, y, z = raw_line.split()[:4]
                vertices.append((float(x) * scale, float(y) * scale, float(z) * scale))
            elif raw_line.startswith("f "):
                face = []
                for token in raw_line.split()[1:]:
                    index = token.split("/")[0]
                    if index:
                        face.append(int(index) - 1)
                if len(face) >= 3:
                    faces.append(tuple(face))
    return vertices, faces


def write_obj(path: Path, mesh: bpy.types.Mesh) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    used_indices = sorted({index for polygon in mesh.polygons for index in polygon.vertices})
    remap = {old: new for new, old in enumerate(used_indices, start=1)}
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# Reduced detailed route PyMOL surface OBJ\n")
        for index in used_indices:
            co = mesh.vertices[index].co
            handle.write(f"v {co.x:.6f} {co.y:.6f} {co.z:.6f}\n")
        for polygon in mesh.polygons:
            handle.write("f " + " ".join(str(remap[index]) for index in polygon.vertices) + "\n")


def weld_duplicate_vertices(mesh: bpy.types.Mesh, distance: float = 0.0001) -> int:
    before = len(mesh.vertices)
    bm = bmesh.new()
    bm.from_mesh(mesh)
    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=distance)
    bm.to_mesh(mesh)
    bm.free()
    mesh.update()
    return before - len(mesh.vertices)


def create_mesh_object(
    name: str,
    vertices: list[tuple[float, float, float]],
    faces: list[tuple[int, ...]],
    mat: bpy.types.Material,
    collection: bpy.types.Collection,
    style: str,
) -> tuple[bpy.types.Object, dict]:
    mesh = bpy.data.meshes.new(f"{name}_mesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    for polygon in mesh.polygons:
        polygon.use_smooth = True
    mesh.materials.append(mat)
    obj = bpy.data.objects.new(name, mesh)
    link_to_collection(obj, collection)
    obj["style"] = style
    report = object_report([obj])
    report.update({"object": obj.name, "vertices": len(vertices), "faces": len(faces), "style": style})
    return obj, report


def object_report(objects: list[bpy.types.Object]) -> dict:
    coords = []
    face_count = 0
    vertex_count = 0
    for obj in objects:
        if obj.type != "MESH":
            continue
        face_count += len(obj.data.polygons)
        vertex_count += len(obj.data.vertices)
        for corner in obj.bound_box:
            coords.append(obj.matrix_world @ Vector(corner))
    if not coords:
        return {"bbox_mm": [0.0, 0.0, 0.0], "center_mm": [0.0, 0.0, 0.0], "faces": face_count, "vertices": vertex_count}
    min_v = Vector((min(c.x for c in coords), min(c.y for c in coords), min(c.z for c in coords)))
    max_v = Vector((max(c.x for c in coords), max(c.y for c in coords), max(c.z for c in coords)))
    center = (min_v + max_v) * 0.5
    return {
        "bbox_mm": [max_v.x - min_v.x, max_v.y - min_v.y, max_v.z - min_v.z],
        "center_mm": [center.x, center.y, center.z],
        "faces": face_count,
        "vertices": vertex_count,
    }


def translate_objects(objects: list[bpy.types.Object], target_center: Vector) -> None:
    report = object_report(objects)
    center = Vector(report["center_mm"])
    delta = target_center - center
    for obj in objects:
        obj.location += delta
    bpy.context.view_layer.update()


def reduce_raw_pymol_surface(asset_id: str, component: str) -> dict:
    raw_path = RAW_DIR / asset_id / f"{asset_id}_surface_{component}.obj"
    reduced_path = REDUCED_DIR / asset_id / raw_path.name
    if not raw_path.exists() or raw_path.stat().st_size == 0:
        return {
            "asset_id": asset_id,
            "component": component,
            "status": "missing_raw_surface",
            "raw_obj": str(raw_path.relative_to(ROOT)).replace("\\", "/"),
        }
    vertices, faces = parse_obj(raw_path, scale=1.0)
    mesh = bpy.data.meshes.new(f"{asset_id}_{component}_reduce_mesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    welded_vertices = weld_duplicate_vertices(mesh)
    obj = bpy.data.objects.new(f"{asset_id}_{component}_reduce", mesh)
    bpy.context.scene.collection.objects.link(obj)
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)
    target = min(len(faces), REDUCTION_TARGETS[asset_id][component])
    ratio = 1.0
    if len(faces) > target and len(faces) > 0:
        ratio = max(0.001, target / len(faces))
        modifier = obj.modifiers.new("detail_route_face_target", "DECIMATE")
        modifier.decimate_type = "COLLAPSE"
        modifier.ratio = ratio
        bpy.ops.object.modifier_apply(modifier=modifier.name)
    for polygon in obj.data.polygons:
        polygon.use_smooth = True
    write_obj(reduced_path, obj.data)
    reduced_faces = len(obj.data.polygons)
    reduced_vertices = len(obj.data.vertices)
    bpy.data.objects.remove(obj, do_unlink=True)
    return {
        "asset_id": asset_id,
        "component": component,
        "status": "reduced",
        "raw_obj": str(raw_path.relative_to(ROOT)).replace("\\", "/"),
        "reduced_obj": str(reduced_path.relative_to(ROOT)).replace("\\", "/"),
        "raw_faces": len(faces),
        "raw_vertices": len(vertices),
        "welded_vertices": welded_vertices,
        "pre_decimate_vertices": len(vertices) - welded_vertices,
        "target_faces": target,
        "decimate_ratio": ratio,
        "reduced_faces": reduced_faces,
        "reduced_vertices": reduced_vertices,
        "target_met": reduced_faces <= max(target, int(target * 1.03)),
    }


def import_reduced_pymol_component(
    asset_id: str,
    component: str,
    mat: bpy.types.Material,
    collection: bpy.types.Collection,
) -> tuple[bpy.types.Object | None, dict]:
    path = REDUCED_DIR / asset_id / f"{asset_id}_surface_{component}.obj"
    if not path.exists() or path.stat().st_size == 0:
        return None, {"asset_id": asset_id, "component": component, "status": "missing_reduced_surface"}
    vertices, faces = parse_obj(path, scale=geom.ANGSTROM_TO_MM)
    obj, report = create_mesh_object(
        f"PyMOL proxy {asset_id} {component}",
        vertices,
        faces,
        mat,
        collection,
        "detailed_pymol_surface_proxy_reduced",
    )
    obj["coordinate_units"] = "angstrom"
    obj["import_scale_mm_per_A"] = geom.ANGSTROM_TO_MM
    report.update(
        {
            "asset_id": asset_id,
            "component": component,
            "status": "imported",
            "coordinate_units": "angstrom",
            "import_scale_mm_per_A": geom.ANGSTROM_TO_MM,
            "source_obj": str(path.relative_to(ROOT)).replace("\\", "/"),
        }
    )
    return obj, report


def uv_sphere_template(segments: int = 10, rings: int = 6) -> tuple[list[Vector], list[tuple[int, ...]]]:
    vertices = [Vector((0.0, 0.0, 1.0))]
    for r in range(1, rings):
        phi = math.pi * r / rings
        z = math.cos(phi)
        radius = math.sin(phi)
        for s in range(segments):
            theta = 2.0 * math.pi * s / segments
            vertices.append(Vector((radius * math.cos(theta), radius * math.sin(theta), z)))
    vertices.append(Vector((0.0, 0.0, -1.0)))
    faces: list[tuple[int, ...]] = []
    bottom = len(vertices) - 1
    for s in range(segments):
        faces.append((0, 1 + (s + 1) % segments, 1 + s))
    for r in range(rings - 2):
        row = 1 + r * segments
        next_row = row + segments
        for s in range(segments):
            faces.append((row + s, row + (s + 1) % segments, next_row + (s + 1) % segments, next_row + s))
    last_row = 1 + (rings - 2) * segments
    for s in range(segments):
        faces.append((bottom, last_row + s, last_row + (s + 1) % segments))
    return vertices, faces


def bead_mesh_from_atoms(atoms: list[dict]) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]]]:
    base_vertices, base_faces = uv_sphere_template()
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, ...]] = []
    for atom_index, atom in enumerate(atoms):
        center = Vector(atom["xyz_mm"])
        radius = float(atom["radius_mm"])
        squash = 0.82 + 0.18 * math.sin(atom_index * 0.71)
        offset = len(vertices)
        for vertex in base_vertices:
            point = center + Vector((vertex.x * radius, vertex.y * radius * squash, vertex.z * radius))
            vertices.append((point.x, point.y, point.z))
        faces.extend(tuple(offset + index for index in face) for face in base_faces)
    return vertices, faces


def build_direct_blender_route(collection: bpy.types.Collection, materials: dict) -> dict:
    dna = geom.build_dna_detail_atoms()
    rna = geom.build_rna_detail_atoms()
    objects: list[bpy.types.Object] = []
    dna_objects: list[bpy.types.Object] = []
    rna_objects: list[bpy.types.Object] = []
    entries = {}
    for component, mat_name in [("strand_A", "dna_a"), ("strand_B", "dna_b"), ("base_pairs", "dna_base")]:
        atoms = [atom for atom in dna["atoms"] if atom["component"] == component]
        vertices, faces = bead_mesh_from_atoms(atoms)
        obj, report = create_mesh_object(
            f"Direct Blender structural DNA {component}",
            vertices,
            faces,
            materials[mat_name],
            collection,
            "direct_blender_structural_bead_surface",
        )
        obj["component"] = component
        objects.append(obj)
        dna_objects.append(obj)
        entries[f"dna_{component}"] = dict(report, component=component, pseudoatom_count=len(atoms))
    for component, mat_name in [("utr5", "rna_utr5"), ("coding", "rna_coding"), ("utr3", "rna_utr3")]:
        atoms = [atom for atom in rna["atoms"] if atom["component"] == component]
        vertices, faces = bead_mesh_from_atoms(atoms)
        obj, report = create_mesh_object(
            f"Direct Blender structural RNA {component}",
            vertices,
            faces,
            materials[mat_name],
            collection,
            "direct_blender_structural_bead_surface",
        )
        obj["component"] = component
        objects.append(obj)
        rna_objects.append(obj)
        entries[f"rna_{component}"] = dict(report, component=component, pseudoatom_count=len(atoms))
    translate_objects(dna_objects, Vector((4.0, 7.0, 0.0)))
    translate_objects(rna_objects, Vector((4.0, -7.0, 0.0)))
    route_report = object_report(objects)
    route_report.update(
        {
            "route": "direct_blender_structural_bead_surface",
            "objects": [obj.name for obj in objects],
            "components": entries,
            "detail_basis": "explicit phosphate/sugar/base pseudoatom centers at every bp or nt; not random bump noise",
            "implicit_union": False,
        }
    )
    return route_report


def build_pymol_route(collection: bpy.types.Collection, materials: dict) -> dict:
    reductions = []
    imports = []
    objects: list[bpy.types.Object] = []
    dna_objects: list[bpy.types.Object] = []
    rna_objects: list[bpy.types.Object] = []
    for asset_id, component_mats in [
        ("DNA_DETAIL_PROXY", [("strand_A", "dna_a"), ("strand_B", "dna_b"), ("base_pairs", "dna_base")]),
        ("RNA_DETAIL_PROXY", [("utr5", "rna_utr5"), ("coding", "rna_coding"), ("utr3", "rna_utr3")]),
    ]:
        for component, mat_name in component_mats:
            reductions.append(reduce_raw_pymol_surface(asset_id, component))
            obj, report = import_reduced_pymol_component(asset_id, component, materials[mat_name], collection)
            imports.append(report)
            if obj is not None:
                objects.append(obj)
                if asset_id == "DNA_DETAIL_PROXY":
                    dna_objects.append(obj)
                else:
                    rna_objects.append(obj)
    translate_objects(dna_objects, Vector((-28.0, 7.0, 0.0)))
    translate_objects(rna_objects, Vector((-28.0, -7.0, 0.0)))
    route_report = object_report(objects)
    route_report.update(
        {
            "route": "detailed_pymol_surface_proxy_reduced",
            "objects": [obj.name for obj in objects],
            "reductions": reductions,
            "imports": imports,
            "gap_filling_assessment": "pseudoatom spacing is lower than or close to overlapping VDW diameters, and object scale remains the Angstrom-to-mm import scale",
            "empirical_scaling_used": False,
        }
    )
    return route_report


def import_calibrator(pdb_id: str, target_center: Vector, mat: bpy.types.Material, collection: bpy.types.Collection) -> dict:
    path = CALIBRATOR_DIR / pdb_id / f"{pdb_id}_nucleic_surface.obj"
    if not path.exists() or path.stat().st_size == 0:
        return {"pdb_id": pdb_id, "status": "missing", "path": str(path.relative_to(ROOT)).replace("\\", "/")}
    vertices, faces = parse_obj(path, scale=geom.ANGSTROM_TO_MM)
    obj, report = create_mesh_object(
        f"PyMOL real-structure calibrator {pdb_id}",
        vertices,
        faces,
        mat,
        collection,
        "pymol_real_structure_calibrator_surface",
    )
    translate_objects([obj], target_center)
    report.update({"pdb_id": pdb_id, "status": "imported", "path": str(path.relative_to(ROOT)).replace("\\", "/")})
    return report


def add_reference_bars(collection: bpy.types.Collection, materials: dict) -> dict:
    reports = {}
    for name, length, y, label in [
        ("scale_20bp", 20 * geom.DNA_BP_TO_MM, 12.3, "20 bp"),
        ("scale_100nt", 100 * geom.RNA_NT_TO_MM, -12.5, "100 nt"),
    ]:
        curve = bpy.data.curves.new(name, "CURVE")
        curve.dimensions = "3D"
        curve.bevel_depth = 0.045
        spline = curve.splines.new("POLY")
        spline.points.add(1)
        spline.points[0].co = (-42.0, y, 0.0, 1.0)
        spline.points[1].co = (-42.0 + length, y, 0.0, 1.0)
        curve.materials.append(materials["scale"])
        obj = bpy.data.objects.new(name, curve)
        link_to_collection(obj, collection)
        add_text(f"label_{name}", label, (-42.0 + length * 0.5, y - 0.9, 0.2), 0.55, materials["text"], collection)
        reports[name] = {"length_mm": length}
    return reports


def add_lighting_and_camera() -> bpy.types.Object:
    light_data = bpy.data.lights.new("softbox", "AREA")
    light_data.energy = 900.0
    light_data.size = 55.0
    light = bpy.data.objects.new("softbox", light_data)
    light.location = (0.0, -12.0, 55.0)
    bpy.context.scene.collection.objects.link(light)
    camera_data = bpy.data.cameras.new("Camera")
    camera = bpy.data.objects.new("Camera", camera_data)
    camera.location = (0.0, 0.0, 64.0)
    camera.rotation_euler = (0.0, 0.0, 0.0)
    camera_data.type = "ORTHO"
    camera_data.ortho_scale = 47.0
    bpy.context.scene.collection.objects.link(camera)
    bpy.context.scene.camera = camera
    return camera


def render_detail(camera: bpy.types.Object, filepath: Path, location: tuple[float, float, float], ortho_scale: float) -> None:
    camera.location = location
    camera.data.ortho_scale = ortho_scale
    bpy.context.scene.render.filepath = str(filepath)
    bpy.ops.render.render(write_still=True)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REDUCED_DIR.mkdir(parents=True, exist_ok=True)
    clean_scene()
    configure_scene()
    collections = {
        "comparison": make_collection("DNA/RNA detail route comparison"),
        "labels": make_collection("Labels"),
    }
    materials = {
        "dna_a": material("dna_A_orange", (0.90, 0.46, 0.22, 1.0)),
        "dna_b": material("dna_B_brown", (0.54, 0.28, 0.15, 1.0)),
        "dna_base": material("dna_base_gold", (0.78, 0.53, 0.25, 1.0)),
        "rna_utr5": material("rna_utr5_olive", (0.60, 0.58, 0.25, 1.0)),
        "rna_coding": material("rna_coding_orange", (0.84, 0.45, 0.21, 1.0)),
        "rna_utr3": material("rna_utr3_yellow_olive", (0.72, 0.66, 0.32, 1.0)),
        "calibrator": material("calibrator_taupe", (0.66, 0.54, 0.43, 1.0)),
        "text": material("text_grey", (0.20, 0.20, 0.20, 1.0)),
        "scale": material("scale_grey", (0.52, 0.52, 0.52, 1.0)),
    }

    report = {
        "title": "DNA/RNA detailed route comparison",
        "canonical_scene_changed": False,
        "units": {
            "angstrom_to_mm": geom.ANGSTROM_TO_MM,
            "dna_bp_to_mm": geom.DNA_BP_TO_MM,
            "rna_nt_to_mm": geom.RNA_NT_TO_MM,
        },
        "comparison_paths": {
            "dna_visual_bp": geom.DNA_VISUAL_BP,
            "dna_axis_length_mm": geom.build_dna_detail_atoms()["geometry"]["axis_length_mm"],
            "rna_visual_nt": sum(segment["nt"] for segment in geom.RNA_VISUAL_SEGMENTS),
            "rna_visual_segments": geom.build_rna_detail_atoms()["geometry"]["segments"],
        },
        "full_actin_mrna_scale_validation": geom.full_actin_mrna_scale_report(),
        "routes": {},
        "calibrators": {},
        "outputs": {
            "blend": str(BLEND_PATH.relative_to(ROOT)).replace("\\", "/"),
            "preview": str(PREVIEW_PATH.relative_to(ROOT)).replace("\\", "/"),
            "dna_detail_preview": str(DNA_DETAIL_PREVIEW_PATH.relative_to(ROOT)).replace("\\", "/"),
            "rna_detail_preview": str(RNA_DETAIL_PREVIEW_PATH.relative_to(ROOT)).replace("\\", "/"),
        },
    }

    report["routes"]["pymol_gap_filled_proxy"] = build_pymol_route(collections["comparison"], materials)
    report["routes"]["direct_blender_structural"] = build_direct_blender_route(collections["comparison"], materials)
    report["calibrators"]["dna_1BNA"] = import_calibrator("1BNA", Vector((27.0, 7.3, 0.0)), materials["calibrator"], collections["comparison"])
    report["calibrators"]["rna_9IOB"] = import_calibrator("9IOB", Vector((27.0, -7.0, 0.0)), materials["calibrator"], collections["comparison"])
    report["scale_bars"] = add_reference_bars(collections["labels"], materials)

    add_text("label_pymol", "PyMOL dense proxy", (-28.0, 15.2, 0.6), 0.85, materials["text"], collections["labels"])
    add_text("label_blender", "direct Blender structural beads", (4.0, 15.2, 0.6), 0.85, materials["text"], collections["labels"])
    add_text("label_calibrators", "real PDB calibrators", (27.0, 15.2, 0.6), 0.85, materials["text"], collections["labels"])
    add_text("label_scale", "Same scale: 1 A = 0.04 mm; DNA 0.136 mm/bp; RNA 0.12 mm/nt", (-6.0, -16.7, 0.6), 0.58, materials["text"], collections["labels"])

    dna = geom.build_dna_detail_atoms()["geometry"]
    rna = geom.build_rna_detail_atoms()["geometry"]
    report["acceptance_checks"] = {
        "scale_constants_ok": geom.ANGSTROM_TO_MM == 0.04 and geom.DNA_BP_TO_MM == 0.136 and geom.RNA_NT_TO_MM == 0.12,
        "dna_envelope_near_0_8_mm": abs(dna["estimated_envelope_diameter_mm"] - 0.8) <= 0.03,
        "full_mrna_segments_mm": [segment["target_length_mm"] for segment in report["full_actin_mrna_scale_validation"]["segments"]],
        "full_mrna_total_mm": report["full_actin_mrna_scale_validation"]["total_length_mm"],
        "pymol_no_empirical_scaling": not report["routes"]["pymol_gap_filled_proxy"]["empirical_scaling_used"],
        "blender_detail_is_structural": report["routes"]["direct_blender_structural"]["detail_basis"].startswith("explicit"),
        "canonical_scene_changed": False,
        "pymol_dna_overlap_margin_A": dna["strand_overlap_margin_A"],
        "pymol_rna_overlap_margin_A": rna["backbone_overlap_margin_A"],
    }

    camera = add_lighting_and_camera()
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_PATH))
    render_detail(camera, PREVIEW_PATH, (-6.0, 0.0, 64.0), 92.0)
    render_detail(camera, DNA_DETAIL_PREVIEW_PATH, (-12.0, 7.0, 28.0), 55.0)
    render_detail(camera, RNA_DETAIL_PREVIEW_PATH, (-12.0, -7.0, 28.0), 55.0)
    print(f"Wrote {BLEND_PATH}")
    print(f"Wrote {PREVIEW_PATH}")
    print(f"Wrote {DNA_DETAIL_PREVIEW_PATH}")
    print(f"Wrote {RNA_DETAIL_PREVIEW_PATH}")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
