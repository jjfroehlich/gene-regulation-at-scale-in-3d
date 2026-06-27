#!/usr/bin/env python3
"""Build the DNA/RNA optimization V3 comparison scene."""

from __future__ import annotations

import json
import math
import os
import sys
from collections import Counter
from pathlib import Path

import bmesh
import bpy
from mathutils import Vector
from mathutils.kdtree import KDTree


ROOT = Path(os.environ.get("GENE_SCENE_ROOT", Path(__file__).resolve().parents[3])).resolve()
SCRIPT_DIR = ROOT / "experiments" / "procedural_nucleic_acids" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))
import optimization_v3_geometry as opt_geom  # noqa: E402


EXPERIMENT_DIR = ROOT / "experiments" / "procedural_nucleic_acids"
ASSET_DIR = EXPERIMENT_DIR / "assets"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
OPT_DIR = ASSET_DIR / "optimization_v3"
OPT_RAW_DIR = OPT_DIR / "raw_surfaces"
OPT_REDUCED_DIR = OPT_DIR / "reduced_surfaces"
CANONICAL_REDUCED_DIR = ROOT / "assets" / "pymol_exports" / "surface_assets_reduced"
CALIBRATOR_DIR = ASSET_DIR / "pymol_calibrator_surfaces"
ASSET_REPORT_PATH = OPT_DIR / "optimization_v3_assets_report.json"

BLEND_PATH = OUTPUT_DIR / "dna_rna_optimization_v3.blend"
PREVIEW_PATH = OUTPUT_DIR / "preview_dna_rna_optimization_v3.png"
DNA_DETAIL_PATH = OUTPUT_DIR / "preview_dna_rna_optimization_v3_dna_detail.png"
ELONGATED_RNA_DETAIL_PATH = OUTPUT_DIR / "preview_dna_rna_optimization_v3_elongated_rna_detail.png"
COMPACT_RNA_DETAIL_PATH = OUTPUT_DIR / "preview_dna_rna_optimization_v3_compact_rna_detail.png"
CALIBRATOR_DETAIL_PATH = OUTPUT_DIR / "preview_dna_rna_optimization_v3_calibrators_detail.png"
REPORT_PATH = OUTPUT_DIR / "dna_rna_optimization_v3_report.json"

ANGSTROM_TO_MM = 0.04
PANEL_DISPLAY_SCALE = 1.0
IMPORT_SCALE = ANGSTROM_TO_MM * PANEL_DISPLAY_SCALE
BACKGROUND_COLOR = (0.985, 0.985, 0.965, 1.0)
BASELINE_VOXEL_MM = 0.106
OPTIMIZED_VOXEL_MM = 0.075
DNA_BASE_STRICT_MARGIN_MM = 0.040
DNA_BASE_SOFT_MARGIN_MM = 0.060
DNA_STRAND_SYMMETRY_MM = 0.22
DNA_MATERIAL_ISLAND_MIN_FACES = 60

DNA_BASELINE_COMPONENTS = [("strand_A", "dna_a"), ("strand_B", "dna_b"), ("base_pairs", "dna_base")]
RNA_COMPONENTS = [("utr5", "rna_utr5"), ("coding", "rna_coding"), ("utr3", "rna_utr3")]
OPT_TARGETS = {
    opt_geom.DNA_ID: {"strand_A": 80_000, "strand_B": 80_000, "base_pairs": 110_000},
    opt_geom.MRNA_ELONGATED_ID: {"utr5": 20_000, "coding": 105_000, "utr3": 70_000},
    opt_geom.MRNA_COMPACT_ID: {"utr5": 20_000, "coding": 105_000, "utr3": 70_000},
}


def clean_scene() -> None:
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete()
    for collection in list(bpy.data.collections):
        bpy.data.collections.remove(collection)
    for mesh in list(bpy.data.meshes):
        if mesh.users == 0:
            bpy.data.meshes.remove(mesh)


def configure_scene() -> None:
    scene = bpy.context.scene
    scene.unit_settings.system = "METRIC"
    scene.unit_settings.length_unit = "MILLIMETERS"
    scene.unit_settings.scale_length = 0.001
    scene.render.engine = "BLENDER_WORKBENCH"
    scene.render.resolution_x = 3200
    scene.render.resolution_y = 1850
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
        bsdf.inputs["Metallic"].default_value = 0.0
    return mat


def add_text(name: str, body: str, location: tuple[float, float, float], size: float, mat: bpy.types.Material, collection: bpy.types.Collection) -> None:
    curve = bpy.data.curves.new(name, "FONT")
    curve.body = body
    curve.align_x = "CENTER"
    curve.align_y = "CENTER"
    curve.size = size
    curve.extrude = 0.015
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
        handle.write("# DNA/RNA optimization V3 reduced surface OBJ\n")
        for index in used_indices:
            co = mesh.vertices[index].co
            handle.write(f"v {co.x:.6f} {co.y:.6f} {co.z:.6f}\n")
        for polygon in mesh.polygons:
            handle.write("f " + " ".join(str(remap[index]) for index in polygon.vertices) + "\n")


def create_mesh_object(
    name: str,
    vertices: list[tuple[float, float, float]],
    faces: list[tuple[int, ...]],
    mat: bpy.types.Material,
    collection: bpy.types.Collection,
    style: str,
) -> bpy.types.Object:
    mesh = bpy.data.meshes.new(f"{name}_mesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    for polygon in mesh.polygons:
        polygon.use_smooth = True
    mesh.materials.append(mat)
    obj = bpy.data.objects.new(name, mesh)
    obj["style"] = style
    link_to_collection(obj, collection)
    return obj


def weld_duplicate_vertices(mesh: bpy.types.Mesh, distance: float = 0.0001) -> int:
    before = len(mesh.vertices)
    bm = bmesh.new()
    bm.from_mesh(mesh)
    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=distance)
    bm.to_mesh(mesh)
    bm.free()
    mesh.update()
    return before - len(mesh.vertices)


def mesh_data(obj: bpy.types.Object) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]]]:
    return [(v.co.x, v.co.y, v.co.z) for v in obj.data.vertices], [tuple(poly.vertices) for poly in obj.data.polygons]


def topology(vertices: list[tuple[float, float, float]], faces: list[tuple[int, ...]]) -> dict:
    edge_counts: dict[tuple[int, int], int] = {}
    for face in faces:
        for i, a in enumerate(face):
            b = face[(i + 1) % len(face)]
            edge = (a, b) if a < b else (b, a)
            edge_counts[edge] = edge_counts.get(edge, 0) + 1
    return {
        "vertices": len(vertices),
        "faces": len(faces),
        "boundary_edges": sum(1 for count in edge_counts.values() if count == 1),
        "non_manifold_edges": sum(1 for count in edge_counts.values() if count > 2),
    }


def mesh_bounds(vertices: list[tuple[float, float, float]]) -> dict:
    min_v = Vector((min(v[0] for v in vertices), min(v[1] for v in vertices), min(v[2] for v in vertices)))
    max_v = Vector((max(v[0] for v in vertices), max(v[1] for v in vertices), max(v[2] for v in vertices)))
    center = (min_v + max_v) * 0.5
    return {
        "min_mm": [min_v.x, min_v.y, min_v.z],
        "max_mm": [max_v.x, max_v.y, max_v.z],
        "bbox_mm": [max_v.x - min_v.x, max_v.y - min_v.y, max_v.z - min_v.z],
        "center_mm": [center.x, center.y, center.z],
    }


def translate_object_to(obj: bpy.types.Object, target_center: Vector) -> None:
    vertices, _faces = mesh_data(obj)
    center = Vector(mesh_bounds(vertices)["center_mm"])
    delta = target_center - center
    for vertex in obj.data.vertices:
        vertex.co += delta
    obj.data.update()


def apply_weighted_normals(obj: bpy.types.Object) -> None:
    bpy.ops.object.select_all(action="DESELECT")
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.shade_smooth()
    modifier = obj.modifiers.new("weighted_normals", "WEIGHTED_NORMAL")
    modifier.keep_sharp = True
    modifier.weight = 50
    bpy.ops.object.modifier_apply(modifier=modifier.name)


def reduce_raw_component(asset_id: str, component: str) -> dict:
    raw_path = OPT_RAW_DIR / asset_id / f"{asset_id}_surface_{component}.obj"
    if not raw_path.exists() or raw_path.stat().st_size == 0:
        raise FileNotFoundError(raw_path)
    reduced_path = OPT_REDUCED_DIR / asset_id / raw_path.name
    vertices, faces = parse_obj(raw_path)
    raw_stats = topology(vertices, faces)
    obj = create_mesh_object(f"reduce_{asset_id}_{component}", vertices, faces, material("temporary_reduce_mat", (0.8, 0.8, 0.8, 1.0)), bpy.context.scene.collection, "temporary_reduce")
    welded = weld_duplicate_vertices(obj.data)
    target = OPT_TARGETS.get(asset_id, {}).get(component, len(obj.data.polygons))
    ratio = 1.0
    if len(obj.data.polygons) > target:
        ratio = max(0.001, target / len(obj.data.polygons))
        decimate = obj.modifiers.new("component_face_target", "DECIMATE")
        decimate.decimate_type = "COLLAPSE"
        decimate.ratio = ratio
        bpy.ops.object.select_all(action="DESELECT")
        obj.select_set(True)
        bpy.context.view_layer.objects.active = obj
        bpy.ops.object.modifier_apply(modifier=decimate.name)
    write_obj(reduced_path, obj.data)
    vertices_after, faces_after = mesh_data(obj)
    reduced_stats = topology(vertices_after, faces_after)
    bpy.data.objects.remove(obj, do_unlink=True)
    return {
        "asset_id": asset_id,
        "component": component,
        "raw_obj": str(raw_path.relative_to(ROOT)).replace("\\", "/"),
        "reduced_obj": str(reduced_path.relative_to(ROOT)).replace("\\", "/"),
        "raw": raw_stats,
        "welded_vertices": welded,
        "target_faces": target,
        "decimate_ratio": ratio,
        "reduced": reduced_stats,
    }


def reduced_path_for(asset_id: str, component: str, source: str) -> Path:
    if source == "canonical":
        path = CANONICAL_REDUCED_DIR / asset_id / f"{asset_id}_surface_{component}.obj"
    else:
        path = OPT_REDUCED_DIR / asset_id / f"{asset_id}_surface_{component}.obj"
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    return path


def merge_components(asset_id: str, components: list[tuple[str, str]], source: str) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]], dict[str, list[tuple[float, float, float]]], list[dict]]:
    merged_vertices: list[tuple[float, float, float]] = []
    merged_faces: list[tuple[int, ...]] = []
    component_vertices: dict[str, list[tuple[float, float, float]]] = {}
    reports = []
    for component, mat_key in components:
        path = reduced_path_for(asset_id, component, source)
        vertices, faces = parse_obj(path, IMPORT_SCALE)
        offset = len(merged_vertices)
        merged_vertices.extend(vertices)
        merged_faces.extend(tuple(index + offset for index in face) for face in faces)
        component_vertices[mat_key] = vertices
        reports.append({"component": component, "material": mat_key, "source_obj": str(path.relative_to(ROOT)).replace("\\", "/"), **topology(vertices, faces)})
    return merged_vertices, merged_faces, component_vertices, reports


def apply_voxel_polish(obj: bpy.types.Object, voxel_size: float, smooth_iterations: int, smooth_factor: float) -> dict:
    before_vertices, before_faces = mesh_data(obj)
    before = topology(before_vertices, before_faces)
    bpy.ops.object.select_all(action="DESELECT")
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj
    remesh = obj.modifiers.new("voxel_close_small_holes", "REMESH")
    remesh.mode = "VOXEL"
    remesh.voxel_size = voxel_size
    remesh.adaptivity = 0.0
    remesh.use_smooth_shade = True
    bpy.ops.object.modifier_apply(modifier=remesh.name)
    if smooth_iterations:
        smooth = obj.modifiers.new("light_surface_relax", "SMOOTH")
        smooth.factor = smooth_factor
        smooth.iterations = smooth_iterations
        bpy.ops.object.modifier_apply(modifier=smooth.name)
    apply_weighted_normals(obj)
    after_vertices, after_faces = mesh_data(obj)
    return {
        "voxel_size_display_mm": voxel_size,
        "voxel_size_scene_mm": voxel_size / PANEL_DISPLAY_SCALE,
        "smooth_iterations": smooth_iterations,
        "smooth_factor": smooth_factor,
        "before": before,
        "after": topology(after_vertices, after_faces),
    }


def kdtrees_for(component_vertices: dict[str, list[tuple[float, float, float]]]) -> dict[str, KDTree]:
    kdtrees = {}
    for mat_key, vertices in component_vertices.items():
        kd = KDTree(len(vertices))
        for index, vertex in enumerate(vertices):
            kd.insert(Vector(vertex), index)
        kd.balance()
        kdtrees[mat_key] = kd
    return kdtrees


def nearest_distance(kd: KDTree, point: Vector) -> float:
    _co, _index, distance = kd.find(point)
    return float(distance)


def polygon_neighbors(mesh: bpy.types.Mesh) -> list[set[int]]:
    edge_users: dict[tuple[int, int], list[int]] = {}
    for polygon in mesh.polygons:
        vertices = list(polygon.vertices)
        for i, a in enumerate(vertices):
            b = vertices[(i + 1) % len(vertices)]
            edge = (a, b) if a < b else (b, a)
            edge_users.setdefault(edge, []).append(polygon.index)
    neighbors = [set() for _ in mesh.polygons]
    for users in edge_users.values():
        if len(users) > 1:
            for index in users:
                neighbors[index].update(other for other in users if other != index)
    return neighbors


def smooth_assignments(mesh: bpy.types.Mesh, assignments: list[int], iterations: int) -> list[int]:
    neighbors = polygon_neighbors(mesh)
    current = assignments[:]
    for _iteration in range(iterations):
        updated = current[:]
        for index, face_neighbors in enumerate(neighbors):
            if not face_neighbors:
                continue
            counts = Counter(current[neighbor] for neighbor in face_neighbors)
            majority, majority_count = counts.most_common(1)[0]
            if majority != current[index] and majority_count >= max(3, int(len(face_neighbors) * 0.62)):
                updated[index] = majority
        current = updated
    return current


def cleanup_small_material_islands(
    mesh: bpy.types.Mesh,
    assignments: list[int],
    min_faces: int,
    iterations: int,
) -> tuple[list[int], int]:
    current = assignments[:]
    changed_total = 0
    for _iteration in range(iterations):
        neighbors = polygon_neighbors(mesh)
        seen = set()
        updated = current[:]
        changed_this_round = 0
        for start_index, material_index in enumerate(current):
            if start_index in seen:
                continue
            island = []
            stack = [start_index]
            seen.add(start_index)
            while stack:
                index = stack.pop()
                island.append(index)
                for neighbor in neighbors[index]:
                    if neighbor not in seen and current[neighbor] == material_index:
                        seen.add(neighbor)
                        stack.append(neighbor)
            if len(island) >= min_faces:
                continue
            adjacent = Counter(
                current[neighbor]
                for index in island
                for neighbor in neighbors[index]
                if current[neighbor] != material_index
            )
            if not adjacent:
                continue
            replacement = adjacent.most_common(1)[0][0]
            for index in island:
                updated[index] = replacement
            changed_this_round += len(island)
        current = updated
        changed_total += changed_this_round
        if changed_this_round == 0:
            break
    return current, changed_total


def material_islands(mesh: bpy.types.Mesh, assignments: list[int], material_keys: list[str]) -> dict[str, int]:
    neighbors = polygon_neighbors(mesh)
    seen = set()
    counts = {key: 0 for key in material_keys}
    for start_index, material_index in enumerate(assignments):
        if start_index in seen:
            continue
        counts[material_keys[material_index]] += 1
        stack = [start_index]
        seen.add(start_index)
        while stack:
            index = stack.pop()
            for neighbor in neighbors[index]:
                if neighbor not in seen and assignments[neighbor] == material_index:
                    seen.add(neighbor)
                    stack.append(neighbor)
    return counts


def assign_materials(
    obj: bpy.types.Object,
    component_vertices: dict[str, list[tuple[float, float, float]]],
    material_keys: list[str],
    mats: dict[str, bpy.types.Material],
    kind: str,
) -> dict:
    kd = kdtrees_for(component_vertices)
    obj.data.materials.clear()
    for mat_key in material_keys:
        obj.data.materials.append(mats[mat_key])
    raw_assignments: list[int] = []
    raw_counts = {key: 0 for key in material_keys}
    index_by_key = {key: index for index, key in enumerate(material_keys)}
    for polygon in obj.data.polygons:
        center = polygon.center
        if kind == "dna":
            da = nearest_distance(kd["dna_a"], center)
            db = nearest_distance(kd["dna_b"], center)
            dbase = nearest_distance(kd["dna_base"], center)
            closest_strand = "dna_a" if da <= db else "dna_b"
            strand_min = min(da, db)
            strand_symmetry = abs(da - db)
            # Base-pair color is restricted to faces that sit between the two
            # strands; this keeps strand surfaces from picking up central color.
            base_between_strands = strand_symmetry <= DNA_STRAND_SYMMETRY_MM
            base_strict = dbase + DNA_BASE_STRICT_MARGIN_MM < strand_min
            base_soft = dbase < strand_min + DNA_BASE_SOFT_MARGIN_MM
            if base_between_strands and (base_strict or base_soft):
                key = "dna_base"
            else:
                key = closest_strand
        else:
            key = min(((nearest_distance(kd[mat_key], center), mat_key) for mat_key in material_keys), key=lambda item: item[0])[1]
        raw_assignments.append(index_by_key[key])
        raw_counts[key] += 1
    assignments = smooth_assignments(obj.data, raw_assignments, 2 if kind == "dna" else 1)
    cleaned_faces = 0
    if kind == "dna":
        assignments, cleaned_faces = cleanup_small_material_islands(obj.data, assignments, DNA_MATERIAL_ISLAND_MIN_FACES, 3)
    face_counts = {key: 0 for key in material_keys}
    for polygon, material_index in zip(obj.data.polygons, assignments):
        polygon.material_index = material_index
        face_counts[material_keys[material_index]] += 1
    obj.data.update()
    return {
        "method": "strand_symmetry_gated_component_transfer" if kind == "dna" else "nearest_component_transfer",
        "dna_base_between_strands_only": True if kind == "dna" else None,
        "dna_base_strict_margin_mm": DNA_BASE_STRICT_MARGIN_MM if kind == "dna" else None,
        "dna_base_soft_margin_mm": DNA_BASE_SOFT_MARGIN_MM if kind == "dna" else None,
        "dna_strand_symmetry_mm": DNA_STRAND_SYMMETRY_MM if kind == "dna" else None,
        "small_island_min_faces": DNA_MATERIAL_ISLAND_MIN_FACES if kind == "dna" else None,
        "small_island_cleaned_faces": cleaned_faces,
        "raw_face_counts_by_material": raw_counts,
        "face_counts_by_material": face_counts,
        "material_islands": material_islands(obj.data, assignments, material_keys),
    }


def build_proxy_panel(
    label: str,
    asset_id: str,
    components: list[tuple[str, str]],
    source: str,
    mats: dict[str, bpy.types.Material],
    collection: bpy.types.Collection,
    target_center: Vector,
    kind: str,
    voxel_size: float,
    smooth_iterations: int,
    smooth_factor: float,
) -> dict:
    vertices, faces, component_vertices, component_reports = merge_components(asset_id, components, source)
    obj = create_mesh_object(label, vertices, faces, mats[components[0][1]], collection, "optimization_v3_proxy_panel")
    polish = apply_voxel_polish(obj, voxel_size, smooth_iterations, smooth_factor)
    material_report = assign_materials(obj, component_vertices, [mat_key for _component, mat_key in components], mats, kind)
    translate_object_to(obj, target_center)
    final_vertices, final_faces = mesh_data(obj)
    obj["asset_id"] = asset_id
    obj["panel_label"] = label
    obj["coordinate_units"] = "angstrom_scaled_to_display"
    obj["angstrom_to_mm"] = ANGSTROM_TO_MM
    return {
        "object": obj.name,
        "asset_id": asset_id,
        "source": source,
        "components": component_reports,
        "polish": polish,
        "material_transfer": material_report,
        "final_mesh": {**mesh_bounds(final_vertices), **topology(final_vertices, final_faces)},
    }


def import_calibrator(pdb_id: str, target_center: Vector, mats: dict[str, bpy.types.Material], collection: bpy.types.Collection) -> dict:
    path = CALIBRATOR_DIR / pdb_id / f"{pdb_id}_nucleic_surface.obj"
    vertices, faces = parse_obj(path, IMPORT_SCALE)
    obj = create_mesh_object(f"calibrator {pdb_id}", vertices, faces, mats["calibrator"], collection, "real_pdb_pymol_calibrator_surface")
    if len(obj.data.polygons) > 70_000:
        ratio = 70_000 / len(obj.data.polygons)
        decimate = obj.modifiers.new("calibrator_face_target", "DECIMATE")
        decimate.decimate_type = "COLLAPSE"
        decimate.ratio = ratio
        bpy.ops.object.select_all(action="DESELECT")
        obj.select_set(True)
        bpy.context.view_layer.objects.active = obj
        bpy.ops.object.modifier_apply(modifier=decimate.name)
    apply_weighted_normals(obj)
    translate_object_to(obj, target_center)
    final_vertices, final_faces = mesh_data(obj)
    return {
        "object": obj.name,
        "pdb_id": pdb_id,
        "source_obj": str(path.relative_to(ROOT)).replace("\\", "/"),
        "scale": {"angstrom_to_mm": ANGSTROM_TO_MM, "panel_display_scale": PANEL_DISPLAY_SCALE},
        "final_mesh": {**mesh_bounds(final_vertices), **topology(final_vertices, final_faces)},
    }


def add_camera() -> bpy.types.Object:
    camera_data = bpy.data.cameras.new("Camera")
    camera = bpy.data.objects.new("Camera", camera_data)
    camera.location = (0.0, 0.0, 75.0)
    camera.rotation_euler = (0.0, 0.0, 0.0)
    camera_data.type = "ORTHO"
    camera_data.ortho_scale = 42.0
    bpy.context.scene.collection.objects.link(camera)
    bpy.context.scene.camera = camera
    return camera


def render(camera: bpy.types.Object, path: Path, center: tuple[float, float, float], ortho_scale: float) -> None:
    camera.location = center
    camera.data.ortho_scale = ortho_scale
    bpy.context.scene.camera = camera
    bpy.context.scene.render.filepath = str(path)
    bpy.ops.render.render(write_still=True)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    clean_scene()
    configure_scene()
    comparison = make_collection("DNA RNA optimization V3")
    labels = make_collection("Labels")
    mats = {
        "dna_a": material("dna_A_orange", (0.88, 0.48, 0.24, 1.0)),
        "dna_b": material("dna_B_brown", (0.52, 0.28, 0.15, 1.0)),
        "dna_base": material("dna_base_gold", (0.76, 0.58, 0.30, 1.0)),
        "rna_utr5": material("rna_utr5_olive", (0.60, 0.58, 0.25, 1.0)),
        "rna_coding": material("rna_coding_orange", (0.84, 0.45, 0.21, 1.0)),
        "rna_utr3": material("rna_utr3_yellow_olive", (0.72, 0.66, 0.32, 1.0)),
        "calibrator": material("calibrator_taupe", (0.66, 0.56, 0.44, 1.0)),
        "text": material("text_grey", (0.20, 0.20, 0.20, 1.0)),
    }
    report = {
        "title": "DNA/RNA optimization experiment V3",
        "canonical_scene_changed": False,
        "units": {"angstrom_to_mm": ANGSTROM_TO_MM, "panel_display_scale": PANEL_DISPLAY_SCALE, "display_mm_per_A": IMPORT_SCALE},
        "assets": json.loads(ASSET_REPORT_PATH.read_text(encoding="utf-8")) if ASSET_REPORT_PATH.exists() else {},
        "reduction": [],
        "panels": {},
        "calibrators": {},
        "outputs": {
            "blend": str(BLEND_PATH.relative_to(ROOT)).replace("\\", "/"),
            "preview": str(PREVIEW_PATH.relative_to(ROOT)).replace("\\", "/"),
            "dna_detail": str(DNA_DETAIL_PATH.relative_to(ROOT)).replace("\\", "/"),
            "elongated_rna_detail": str(ELONGATED_RNA_DETAIL_PATH.relative_to(ROOT)).replace("\\", "/"),
            "compact_rna_detail": str(COMPACT_RNA_DETAIL_PATH.relative_to(ROOT)).replace("\\", "/"),
            "calibrator_detail": str(CALIBRATOR_DETAIL_PATH.relative_to(ROOT)).replace("\\", "/"),
        },
    }

    for asset_id, targets in OPT_TARGETS.items():
        for component in targets:
            report["reduction"].append(reduce_raw_component(asset_id, component))

    report["panels"]["baseline_dna"] = build_proxy_panel("baseline canonical DNA", "DNA_PROXY", DNA_BASELINE_COMPONENTS, "canonical", mats, comparison, Vector((-260.0, 90.0, 0.0)), "dna", BASELINE_VOXEL_MM, 3, 0.18)
    report["panels"]["optimized_dna"] = build_proxy_panel("optimized right-handed DNA", opt_geom.DNA_ID, DNA_BASELINE_COMPONENTS, "optimization_v3", mats, comparison, Vector((0.0, 90.0, 0.0)), "dna", OPTIMIZED_VOXEL_MM, 1, 0.10)
    report["panels"]["baseline_mrna"] = build_proxy_panel("baseline canonical mRNA", "MRNA_PROXY", RNA_COMPONENTS, "canonical", mats, comparison, Vector((-260.0, -50.0, 0.0)), "rna", BASELINE_VOXEL_MM, 3, 0.18)
    report["panels"]["irregular_elongated_mrna"] = build_proxy_panel("irregular elongated mRNA", opt_geom.MRNA_ELONGATED_ID, RNA_COMPONENTS, "optimization_v3", mats, comparison, Vector((0.0, -50.0, 0.0)), "rna", OPTIMIZED_VOXEL_MM, 1, 0.10)
    report["panels"]["compact_schematic_mrna"] = build_proxy_panel("compact schematic mRNA", opt_geom.MRNA_COMPACT_ID, RNA_COMPONENTS, "optimization_v3", mats, comparison, Vector((0.0, -180.0, 0.0)), "rna", OPTIMIZED_VOXEL_MM, 1, 0.10)
    report["calibrators"]["dna_1BNA"] = import_calibrator("1BNA", Vector((260.0, 90.0, 0.0)), mats, comparison)
    report["calibrators"]["rna_9IOB"] = import_calibrator("9IOB", Vector((260.0, -50.0, 0.0)), mats, comparison)

    add_text("label_baseline", "current canonical", (-260.0, 155.0, 4.0), 8.0, mats["text"], labels)
    add_text("label_optimized", "optimized V3: right-handed B-DNA", (0.0, 155.0, 4.0), 8.0, mats["text"], labels)
    add_text("label_calibrators", "real PDB calibrators", (260.0, 155.0, 4.0), 8.0, mats["text"], labels)
    add_text("label_compact", "compact schematic mRNA", (0.0, -230.0, 4.0), 8.0, mats["text"], labels)
    add_text("label_scale", "All objects use the project scale: 1 A = 0.04 mm", (0.0, -250.0, 4.0), 6.5, mats["text"], labels)

    dna_geom = report["assets"]["assets"][0]["geometry"]
    if not dna_geom["chirality_matches_1BNA"]:
        raise RuntimeError("Generated DNA chirality does not match 1BNA")
    if not dna_geom["screw_handedness_matches_1BNA"] or not dna_geom["generated_screw_handedness"]["right_handed"]:
        raise RuntimeError("Generated DNA screw handedness is not right-handed against the 1BNA reference")
    for asset in report["assets"]["assets"][1:]:
        for segment in asset["geometry"]["segments"]:
            if abs(segment["error_mm"]) > 0.02:
                raise RuntimeError(f"RNA segment length error too high: {asset['asset_id']} {segment}")
    for key, panel in report["panels"].items():
        if panel["polish"]["after"]["boundary_edges"] != 0 or panel["polish"]["after"]["non_manifold_edges"] != 0:
            raise RuntimeError(f"Polished proxy is not closed: {key}")

    camera = add_camera()
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_PATH))
    render(camera, PREVIEW_PATH, (0.0, -40.0, 700.0), 430.0)
    render(camera, DNA_DETAIL_PATH, (0.0, 90.0, 250.0), 12.0)
    render(camera, ELONGATED_RNA_DETAIL_PATH, (0.0, -50.0, 250.0), 30.0)
    render(camera, COMPACT_RNA_DETAIL_PATH, (0.0, -180.0, 250.0), 34.0)
    render(camera, CALIBRATOR_DETAIL_PATH, (260.0, 20.0, 250.0), 180.0)
    print(f"Wrote {BLEND_PATH}")
    print(f"Wrote {PREVIEW_PATH}")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
