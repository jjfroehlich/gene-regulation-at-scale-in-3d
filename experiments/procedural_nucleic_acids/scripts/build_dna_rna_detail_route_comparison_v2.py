#!/usr/bin/env python3
"""Build V2 DNA/RNA generation route comparison."""

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
sys.path.insert(0, str(ROOT / "scripts"))
SCRIPT_DIR = ROOT / "experiments" / "procedural_nucleic_acids" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))
import legacy_gap_proxy_geometry as gap_geom  # noqa: E402
import procedural_nucleic_geometry as geom  # noqa: E402


EXPERIMENT_DIR = ROOT / "experiments" / "procedural_nucleic_acids"
ASSET_DIR = EXPERIMENT_DIR / "assets"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
BLEND_PATH = OUTPUT_DIR / "dna_rna_detail_route_comparison_v2.blend"
PREVIEW_PATH = OUTPUT_DIR / "preview_dna_rna_detail_route_comparison_v2.png"
DNA_DETAIL_PREVIEW_PATH = OUTPUT_DIR / "preview_dna_rna_detail_route_comparison_v2_dna_detail.png"
RNA_DETAIL_PREVIEW_PATH = OUTPUT_DIR / "preview_dna_rna_detail_route_comparison_v2_rna_detail.png"
REPORT_PATH = OUTPUT_DIR / "dna_rna_detail_route_comparison_v2_report.json"
MANIFEST_PATH = ROOT / "config" / "scene_manifest.json"
LEGACY_REDUCED_DIR = ROOT / "assets" / "pymol_exports" / "surface_assets_reduced"
GAP_RAW_DIR = ASSET_DIR / "legacy_gap_pymol_proxy" / "raw_surfaces"
GAP_REDUCED_DIR = ASSET_DIR / "legacy_gap_pymol_proxy" / "reduced_surfaces"
DETAIL_REDUCED_DIR = ASSET_DIR / "detail_route_pymol_proxy" / "reduced_surfaces"
CALIBRATOR_DIR = ASSET_DIR / "pymol_calibrator_surfaces"

ANGSTROM_TO_MM = 0.04
PANEL_DISPLAY_SCALE = 0.045
DNA_ROW_Y = 3.4
RNA_ROW_Y = -3.4
LABEL_Y = 8.6
SCALE_LABEL_Y = -8.4
CALIBRATOR_X = 24.5
PREVIEW_CENTER_X = 2.25
PREVIEW_ORTHO_SCALE = 32.0
BACKGROUND_COLOR = (0.985, 0.985, 0.965, 1.0)
GAP_TARGETS = {
    "DNA_PROXY_GAP": {"strand_A": 55_000, "strand_B": 55_000, "base_pairs": 90_000},
    "MRNA_PROXY_GAP": {"utr5": 15_000, "coding": 90_000, "utr3": 60_000},
}


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
    scene.render.resolution_x = 2800
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
        handle.write("# V2 reduced legacy gap surface OBJ\n")
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


def mesh_island_count(vertices: list[tuple[float, float, float]], faces: list[tuple[int, ...]]) -> int:
    parent: dict[int, int] = {}

    def find(x: int) -> int:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra = find(a)
        rb = find(b)
        if ra != rb:
            parent[rb] = ra

    for face in faces:
        if not face:
            continue
        first = face[0]
        find(first)
        for index in face[1:]:
            union(first, index)
    coordinate_owner: dict[tuple[float, float, float], int] = {}
    for index, vertex in enumerate(vertices):
        key = (round(vertex[0], 5), round(vertex[1], 5), round(vertex[2], 5))
        if key in coordinate_owner:
            union(coordinate_owner[key], index)
        else:
            coordinate_owner[key] = index
    return len({find(index) for index in parent})


def object_report(objects: list[bpy.types.Object]) -> dict:
    bpy.context.view_layer.update()
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
    center = Vector(object_report(objects)["center_mm"])
    delta = target_center - center
    for obj in objects:
        obj.location += delta
    bpy.context.view_layer.update()


def create_mesh_object(
    name: str,
    vertices: list[tuple[float, float, float]],
    faces: list[tuple[int, ...]],
    mat: bpy.types.Material,
    collection: bpy.types.Collection,
    style: str,
    panel_scale: float | None = None,
) -> tuple[bpy.types.Object, dict]:
    mesh = bpy.data.meshes.new(f"{name}_mesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    for polygon in mesh.polygons:
        polygon.use_smooth = True
    mesh.materials.append(mat)
    obj = bpy.data.objects.new(name, mesh)
    if panel_scale is not None:
        obj.scale = (panel_scale, panel_scale, panel_scale)
    link_to_collection(obj, collection)
    obj["style"] = style
    report = object_report([obj])
    report.update({"object": obj.name, "style": style, "mesh_islands": mesh_island_count(vertices, faces)})
    return obj, report


def reduce_gap_surface(asset_id: str, component: str) -> dict:
    raw_path = GAP_RAW_DIR / asset_id / f"{asset_id}_surface_{component}.obj"
    reduced_path = GAP_REDUCED_DIR / asset_id / raw_path.name
    if not raw_path.exists() or raw_path.stat().st_size == 0:
        return {"asset_id": asset_id, "component": component, "status": "missing_raw_surface", "raw_obj": str(raw_path)}
    vertices, faces = parse_obj(raw_path, scale=1.0)
    target = min(len(faces), GAP_TARGETS[asset_id][component])
    mesh = bpy.data.meshes.new(f"{asset_id}_{component}_reduce_mesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    welded_vertices = weld_duplicate_vertices(mesh)
    obj = bpy.data.objects.new(f"{asset_id}_{component}_reduce", mesh)
    bpy.context.scene.collection.objects.link(obj)
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)
    ratio = 1.0
    if len(faces) > target and len(faces) > 0:
        ratio = max(0.001, target / len(faces))
        modifier = obj.modifiers.new("v2_face_target", "DECIMATE")
        modifier.decimate_type = "COLLAPSE"
        modifier.ratio = ratio
        bpy.ops.object.modifier_apply(modifier=modifier.name)
    for polygon in obj.data.polygons:
        polygon.use_smooth = True
    write_obj(reduced_path, obj.data)
    reduced_faces = len(obj.data.polygons)
    reduced_vertices = len(obj.data.vertices)
    reduced_vertices_coords = [(vertex.co.x, vertex.co.y, vertex.co.z) for vertex in obj.data.vertices]
    reduced_faces_list = [tuple(poly.vertices) for poly in obj.data.polygons]
    reduced_islands = mesh_island_count(reduced_vertices_coords, reduced_faces_list)
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
        "reduced_mesh_islands": reduced_islands,
        "target_met": reduced_faces <= max(target, int(target * 1.03)),
    }


def import_surface_component(
    path: Path,
    name: str,
    mat: bpy.types.Material,
    collection: bpy.types.Collection,
    style: str,
) -> tuple[bpy.types.Object | None, dict]:
    if not path.exists() or path.stat().st_size == 0:
        return None, {"object": name, "status": "missing", "path": str(path)}
    vertices, faces = parse_obj(path, scale=ANGSTROM_TO_MM)
    obj, report = create_mesh_object(name, vertices, faces, mat, collection, style, PANEL_DISPLAY_SCALE)
    report.update(
        {
            "status": "imported",
            "source_obj": str(path.relative_to(ROOT)).replace("\\", "/"),
            "coordinate_import_scale_mm_per_A": ANGSTROM_TO_MM,
            "panel_display_scale": PANEL_DISPLAY_SCALE,
        }
    )
    return obj, report


def import_route(route_name: str, asset_map: list[tuple[str, str, str]], base_dir: Path, materials: dict, collection: bpy.types.Collection, target_x: float) -> dict:
    objects_by_kind = {"dna": [], "rna": []}
    component_reports = []
    for asset_id, component, mat_name in asset_map:
        path = base_dir / asset_id / f"{asset_id}_surface_{component}.obj"
        obj, report = import_surface_component(path, f"{route_name} {asset_id} {component}", materials[mat_name], collection, route_name)
        report.update({"asset_id": asset_id, "component": component})
        component_reports.append(report)
        if obj is not None:
            kind = "dna" if "DNA" in asset_id else "rna"
            objects_by_kind[kind].append(obj)
    translate_objects(objects_by_kind["dna"], Vector((target_x, DNA_ROW_Y, 0.0)))
    translate_objects(objects_by_kind["rna"], Vector((target_x, RNA_ROW_Y, 0.0)))
    route_objects = objects_by_kind["dna"] + objects_by_kind["rna"]
    route_report = object_report(route_objects)
    route_report.update(
        {
            "route": route_name,
            "components": component_reports,
            "objects": [obj.name for obj in route_objects],
            "coordinate_import_scale_mm_per_A": ANGSTROM_TO_MM,
            "panel_display_scale": PANEL_DISPLAY_SCALE,
        }
    )
    return route_report


def build_gap_route(materials: dict, collection: bpy.types.Collection, target_x: float) -> dict:
    reductions = []
    for asset_id, components in [("DNA_PROXY_GAP", ["strand_A", "strand_B", "base_pairs"]), ("MRNA_PROXY_GAP", ["utr5", "coding", "utr3"])]:
        for component in components:
            reductions.append(reduce_gap_surface(asset_id, component))
    route = import_route(
        "legacy_pymol_gap_filled",
        [
            ("DNA_PROXY_GAP", "strand_A", "dna_a"),
            ("DNA_PROXY_GAP", "strand_B", "dna_b"),
            ("DNA_PROXY_GAP", "base_pairs", "dna_base"),
            ("MRNA_PROXY_GAP", "utr5", "rna_utr5"),
            ("MRNA_PROXY_GAP", "coding", "rna_coding"),
            ("MRNA_PROXY_GAP", "utr3", "rna_utr3"),
        ],
        GAP_REDUCED_DIR,
        materials,
        collection,
        target_x,
    )
    route["reductions"] = reductions
    return route


def sphere_template(segments: int = 8, rings: int = 5) -> tuple[list[Vector], list[tuple[int, ...]]]:
    vertices = [Vector((0.0, 0.0, 1.0))]
    for r in range(1, rings):
        phi = math.pi * r / rings
        z = math.cos(phi)
        radius = math.sin(phi)
        for s in range(segments):
            theta = 2.0 * math.pi * s / segments
            vertices.append(Vector((radius * math.cos(theta), radius * math.sin(theta), z)))
    vertices.append(Vector((0.0, 0.0, -1.0)))
    bottom = len(vertices) - 1
    faces: list[tuple[int, ...]] = []
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


def ellipsoid_mesh(centers: list[tuple[Vector, tuple[float, float, float]]]) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]]]:
    base_vertices, base_faces = sphere_template()
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, ...]] = []
    for center, radii in centers:
        offset = len(vertices)
        rx, ry, rz = radii
        for vertex in base_vertices:
            point = center + Vector((vertex.x * rx, vertex.y * ry, vertex.z * rz))
            vertices.append((point.x, point.y, point.z))
        faces.extend(tuple(offset + index for index in face) for face in base_faces)
    return vertices, faces


def add_box(vertices: list[tuple[float, float, float]], faces: list[tuple[int, ...]], center: Vector, u: Vector, v: Vector, w: Vector, half: tuple[float, float, float]) -> None:
    hx, hy, hz = half
    offset = len(vertices)
    corners = []
    for sx, sy, sz in [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1), (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]:
        p = center + u * (sx * hx) + v * (sy * hy) + w * (sz * hz)
        corners.append((p.x, p.y, p.z))
    vertices.extend(corners)
    faces.extend(
        [
            (offset + 0, offset + 1, offset + 2, offset + 3),
            (offset + 4, offset + 7, offset + 6, offset + 5),
            (offset + 0, offset + 4, offset + 5, offset + 1),
            (offset + 1, offset + 5, offset + 6, offset + 2),
            (offset + 2, offset + 6, offset + 7, offset + 3),
            (offset + 3, offset + 7, offset + 4, offset + 0),
        ]
    )


def as_vector(point) -> Vector:
    return Vector((point.x, point.y, point.z)) if hasattr(point, "x") else Vector(point)


def build_direct_realistic_route(manifest: dict, materials: dict, collection: bpy.types.Collection, target_x: float) -> dict:
    dna = geom.build_dna_model(manifest)
    mrna = geom.build_mrna_model(manifest)
    objects = {"dna": [], "rna": []}
    component_reports = {}

    dna_backbone_centers = []
    for strand in (dna["strand_a"], dna["strand_b"]):
        for point in strand:
            center = as_vector(point)
            dna_backbone_centers.append((center, (0.052, 0.044, 0.052)))
            axis_center = Vector((center.x, -45.0, 0.0))
            sugar = center.lerp(axis_center, 0.16)
            dna_backbone_centers.append((sugar, (0.042, 0.034, 0.042)))
    vertices, faces = ellipsoid_mesh(dna_backbone_centers)
    obj, report = create_mesh_object("direct Blender realistic DNA phosphate_sugar", vertices, faces, materials["dna_a"], collection, "direct_blender_realistic", PANEL_DISPLAY_SCALE)
    objects["dna"].append(obj)
    component_reports["dna_phosphate_sugar"] = dict(report, pseudoatom_like_count=len(dna_backbone_centers), detail_basis="phosphate and sugar ellipsoids")

    base_vertices: list[tuple[float, float, float]] = []
    base_faces: list[tuple[int, ...]] = []
    for i, (a_raw, b_raw) in enumerate(zip(dna["strand_a"], dna["strand_b"])):
        a = as_vector(a_raw)
        b = as_vector(b_raw)
        center = (a + b) * 0.5
        u = (b - a).normalized()
        tangent = (as_vector(dna["strand_a"][min(i + 1, len(dna["strand_a"]) - 1)]) - as_vector(dna["strand_a"][max(i - 1, 0)])).normalized()
        w = u.cross(tangent).normalized()
        if w.length < 1e-6:
            w = Vector((0.0, 0.0, 1.0))
        v = w.cross(u).normalized()
        add_box(base_vertices, base_faces, center, u, v, w, (0.235, 0.018, 0.045))
    obj, report = create_mesh_object("direct Blender realistic DNA stacked_base_pairs", base_vertices, base_faces, materials["dna_base"], collection, "direct_blender_realistic", PANEL_DISPLAY_SCALE)
    objects["dna"].append(obj)
    component_reports["dna_stacked_base_pairs"] = dict(report, base_pair_plates=len(dna["strand_a"]), detail_basis="planar stacked base-pair plates")

    nt_to_mm = manifest["units"]["mrna_nt_to_mm"]
    rna_component_reports = []
    for segment_model in mrna["segments"]:
        segment = segment_model["segment"]
        path = geom.SampledPath(segment_model["points"])
        centers = []
        base_centers = []
        for nt in range(segment["nt"]):
            distance = min(path.length, nt * nt_to_mm)
            center = as_vector(path.point_at_length(distance))
            tangent = as_vector(path.tangent_at_length(distance))
            normal = Vector((-tangent.y, tangent.x, 0.0))
            if normal.length < 1e-6:
                normal = Vector((0.0, 1.0, 0.0))
            normal.normalize()
            theta = 2.0 * math.pi * nt / 6.0
            side = normal if math.sin(theta) >= 0 else -normal
            centers.append((center, (0.046, 0.039, 0.046)))
            centers.append((center + side * 0.07, (0.038, 0.030, 0.038)))
            base_centers.append((center + side * 0.18, (0.090, 0.035, 0.060)))
        vertices, faces = ellipsoid_mesh(centers)
        obj, report = create_mesh_object(f"direct Blender realistic RNA {segment['name']} phosphate_sugar", vertices, faces, materials[segment["color"]], collection, "direct_blender_realistic", PANEL_DISPLAY_SCALE)
        objects["rna"].append(obj)
        rna_component_reports.append(dict(report, segment=segment["name"], detail_basis="RNA phosphate and ribose ellipsoids"))
        vertices, faces = ellipsoid_mesh(base_centers)
        obj, report = create_mesh_object(f"direct Blender realistic RNA {segment['name']} bases", vertices, faces, materials[segment["color"]], collection, "direct_blender_realistic", PANEL_DISPLAY_SCALE)
        objects["rna"].append(obj)
        rna_component_reports.append(dict(report, segment=segment["name"], detail_basis="offset RNA base ellipsoids"))
    translate_objects(objects["dna"], Vector((target_x, DNA_ROW_Y, 0.0)))
    translate_objects(objects["rna"], Vector((target_x, RNA_ROW_Y, 0.0)))
    all_objects = objects["dna"] + objects["rna"]
    route_report = object_report(all_objects)
    route_report.update(
        {
            "route": "direct_blender_realistic",
            "objects": [obj.name for obj in all_objects],
            "panel_display_scale": PANEL_DISPLAY_SCALE,
            "components": {**component_reports, "rna_components": rna_component_reports},
            "detail_basis": "B-DNA phosphate/sugar ellipsoids plus planar base-pair plates; ssRNA phosphate/sugar/base ellipsoids",
        }
    )
    return route_report


def import_calibrator(pdb_id: str, target_center: Vector, mat: bpy.types.Material, collection: bpy.types.Collection) -> dict:
    path = CALIBRATOR_DIR / pdb_id / f"{pdb_id}_nucleic_surface.obj"
    if not path.exists():
        return {"pdb_id": pdb_id, "status": "missing", "path": str(path)}
    vertices, faces = parse_obj(path, scale=ANGSTROM_TO_MM)
    obj, report = create_mesh_object(f"real PDB calibrator {pdb_id}", vertices, faces, mat, collection, "pymol_real_structure_calibrator")
    translate_objects([obj], target_center)
    report.update({"pdb_id": pdb_id, "status": "imported", "source_obj": str(path.relative_to(ROOT)).replace("\\", "/"), "coordinate_import_scale_mm_per_A": ANGSTROM_TO_MM})
    return report


def add_camera() -> bpy.types.Object:
    light_data = bpy.data.lights.new("softbox", "AREA")
    light_data.energy = 900.0
    light_data.size = 60.0
    light = bpy.data.objects.new("softbox", light_data)
    light.location = (0.0, -10.0, 65.0)
    bpy.context.scene.collection.objects.link(light)
    camera_data = bpy.data.cameras.new("Camera")
    camera = bpy.data.objects.new("Camera", camera_data)
    camera.location = (10.0, 0.0, 76.0)
    camera.rotation_euler = (0.0, 0.0, 0.0)
    camera_data.type = "ORTHO"
    camera_data.ortho_scale = 72.0
    bpy.context.scene.collection.objects.link(camera)
    bpy.context.scene.camera = camera
    return camera


def render(camera: bpy.types.Object, path: Path, location: tuple[float, float, float], ortho_scale: float) -> None:
    camera.location = location
    camera.data.ortho_scale = ortho_scale
    bpy.context.scene.render.filepath = str(path)
    bpy.ops.render.render(write_still=True)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    GAP_REDUCED_DIR.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    gap_report_path = ASSET_DIR / "legacy_gap_pymol_proxy" / "legacy_gap_proxy_generation_report.json"
    gap_generation_report = json.loads(gap_report_path.read_text(encoding="utf-8")) if gap_report_path.exists() else {}
    clean_scene()
    configure_scene()
    comparison = make_collection("DNA/RNA generation V2 comparison")
    labels = make_collection("Labels")
    materials = {
        "dna_a": material("dna_A_orange", (0.88, 0.48, 0.24, 1.0)),
        "dna_b": material("dna_B_brown", (0.52, 0.28, 0.15, 1.0)),
        "dna_base": material("dna_base_gold", (0.76, 0.55, 0.28, 1.0)),
        "rna_utr5": material("rna_utr5_olive", (0.60, 0.58, 0.25, 1.0)),
        "rna_coding": material("rna_coding_orange", (0.84, 0.45, 0.21, 1.0)),
        "rna_utr3": material("rna_utr3_yellow_olive", (0.72, 0.66, 0.32, 1.0)),
        "olive": material("olive", (0.60, 0.58, 0.25, 1.0)),
        "orange": material("orange", (0.84, 0.45, 0.21, 1.0)),
        "yellow_olive": material("yellow_olive", (0.72, 0.66, 0.32, 1.0)),
        "calibrator": material("calibrator_taupe", (0.62, 0.56, 0.48, 1.0)),
        "text": material("text_grey", (0.20, 0.20, 0.20, 1.0)),
    }

    route_positions = {
        "legacy_pymol_original": -15.0,
        "legacy_pymol_gap_filled": -5.0,
        "direct_blender_realistic": 5.0,
        "new_detail_reference": 15.0,
    }
    report = {
        "title": "DNA/RNA generation route comparison V2",
        "canonical_scene_changed": False,
        "units": {
            "angstrom_to_mm": ANGSTROM_TO_MM,
            "dna_bp_to_mm": manifest["units"]["dna_bp_to_mm"],
            "rna_nt_to_mm": manifest["units"]["mrna_nt_to_mm"],
            "panel_display_scale": PANEL_DISPLAY_SCALE,
            "layout": {
                "route_center_spacing_mm": 10.0,
                "dna_row_y_mm": DNA_ROW_Y,
                "rna_row_y_mm": RNA_ROW_Y,
                "calibrator_x_mm": CALIBRATOR_X,
            },
        },
        "routes": {},
        "calibrators": {},
        "gap_generation_report": str(gap_report_path.relative_to(ROOT)).replace("\\", "/") if gap_report_path.exists() else None,
        "outputs": {
            "blend": str(BLEND_PATH.relative_to(ROOT)).replace("\\", "/"),
            "preview": str(PREVIEW_PATH.relative_to(ROOT)).replace("\\", "/"),
            "dna_detail_preview": str(DNA_DETAIL_PREVIEW_PATH.relative_to(ROOT)).replace("\\", "/"),
            "rna_detail_preview": str(RNA_DETAIL_PREVIEW_PATH.relative_to(ROOT)).replace("\\", "/"),
        },
    }
    report["routes"]["legacy_pymol_original"] = import_route(
        "legacy_pymol_original",
        [
            ("DNA_PROXY", "strand_A", "dna_a"),
            ("DNA_PROXY", "strand_B", "dna_b"),
            ("DNA_PROXY", "base_pairs", "dna_base"),
            ("MRNA_PROXY", "utr5", "rna_utr5"),
            ("MRNA_PROXY", "coding", "rna_coding"),
            ("MRNA_PROXY", "utr3", "rna_utr3"),
        ],
        LEGACY_REDUCED_DIR,
        materials,
        comparison,
        route_positions["legacy_pymol_original"],
    )
    report["routes"]["legacy_pymol_original"]["route_note"] = (
        "Shared DNA_PROXY/MRNA_PROXY surfaces now use old-proxy geometry plus connector pseudoatoms "
        "to reduce PyMOL surface holes while preserving the older visual style."
    )
    report["routes"]["legacy_pymol_gap_filled"] = build_gap_route(materials, comparison, route_positions["legacy_pymol_gap_filled"])
    report["routes"]["direct_blender_realistic"] = build_direct_realistic_route(manifest, materials, comparison, route_positions["direct_blender_realistic"])
    report["routes"]["new_detail_reference"] = import_route(
        "new_detail_reference",
        [
            ("DNA_DETAIL_PROXY", "strand_A", "dna_a"),
            ("DNA_DETAIL_PROXY", "strand_B", "dna_b"),
            ("DNA_DETAIL_PROXY", "base_pairs", "dna_base"),
            ("RNA_DETAIL_PROXY", "utr5", "rna_utr5"),
            ("RNA_DETAIL_PROXY", "coding", "rna_coding"),
            ("RNA_DETAIL_PROXY", "utr3", "rna_utr3"),
        ],
        DETAIL_REDUCED_DIR,
        materials,
        comparison,
        route_positions["new_detail_reference"],
    )
    report["calibrators"]["dna_1BNA"] = import_calibrator("1BNA", Vector((CALIBRATOR_X, DNA_ROW_Y, 0.0)), materials["calibrator"], comparison)
    report["calibrators"]["rna_9IOB"] = import_calibrator("9IOB", Vector((CALIBRATOR_X, RNA_ROW_Y, 0.0)), materials["calibrator"], comparison)

    label_names = {
        "legacy_pymol_original": "old PyMOL repaired",
        "legacy_pymol_gap_filled": "old PyMOL gap-filled",
        "direct_blender_realistic": "direct Blender realistic",
        "new_detail_reference": "new detail reference",
    }
    for key, x in route_positions.items():
        add_text(f"label_{key}", label_names[key], (x, LABEL_Y, 0.5), 0.46, materials["text"], labels)
    add_text("label_calibrators", "real PDB calibrators", (CALIBRATOR_X, LABEL_Y, 0.5), 0.46, materials["text"], labels)
    add_text("label_scale", "Imported at 1 A = 0.04 mm; route panels displayed at 4.5% for side-by-side comparison", (PREVIEW_CENTER_X, SCALE_LABEL_Y, 0.5), 0.34, materials["text"], labels)

    old_dna_islands = sum(c.get("mesh_islands", 0) for c in report["routes"]["legacy_pymol_original"]["components"] if "DNA_PROXY" in c.get("asset_id", ""))
    gap_dna_islands = sum(c.get("mesh_islands", 0) for c in report["routes"]["legacy_pymol_gap_filled"]["components"] if "DNA_PROXY_GAP" in c.get("asset_id", ""))
    old_rna_islands = sum(c.get("mesh_islands", 0) for c in report["routes"]["legacy_pymol_original"]["components"] if "MRNA_PROXY" in c.get("asset_id", ""))
    gap_rna_islands = sum(c.get("mesh_islands", 0) for c in report["routes"]["legacy_pymol_gap_filled"]["components"] if "MRNA_PROXY_GAP" in c.get("asset_id", ""))
    full_mrna = gap_generation_report.get("full_actin_mrna_scale_validation", gap_geom.full_actin_mrna_scale_report(manifest))
    report["acceptance_checks"] = {
        "scale_constants_ok": ANGSTROM_TO_MM == 0.04 and manifest["units"]["dna_bp_to_mm"] == 0.136 and manifest["units"]["mrna_nt_to_mm"] == 0.12,
        "full_mrna_segments_mm": [segment["target_length_mm"] for segment in full_mrna["segments"]],
        "full_mrna_total_mm": full_mrna["total_length_mm"],
        "canonical_scene_changed": False,
        "gap_filled_dna_overlap_margin_A": gap_generation_report.get("assets", [{}])[0].get("geometry", {}).get("gap_filled_strand_overlap_margin_A"),
        "legacy_original_dna_overlap_margin_A": gap_generation_report.get("baseline_metrics", {}).get("legacy_original", {}).get("dna", {}).get("strand_overlap_margin_A"),
        "gap_filled_improves_dna_overlap": (
            gap_generation_report.get("assets", [{}])[0].get("geometry", {}).get("gap_filled_strand_overlap_margin_A", -999)
            > gap_generation_report.get("baseline_metrics", {}).get("legacy_original", {}).get("dna", {}).get("strand_overlap_margin_A", 999)
        ),
        "gap_filled_dna_mesh_islands": gap_dna_islands,
        "legacy_original_dna_mesh_islands": old_dna_islands,
        "gap_filled_rna_mesh_islands": gap_rna_islands,
        "legacy_original_rna_mesh_islands": old_rna_islands,
        "mesh_island_metric_note": (
            "PyMOL OBJ export plus decimation can duplicate or fragment surface topology; "
            "mesh-island counts are retained as diagnostics, while pseudoatom overlap margin "
            "is the primary gap-filling metric."
        ),
        "gap_reductions_target_met": all(item.get("target_met", False) for item in report["routes"]["legacy_pymol_gap_filled"].get("reductions", [])),
    }

    camera = add_camera()
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_PATH))
    render(camera, PREVIEW_PATH, (PREVIEW_CENTER_X, 0.0, 76.0), PREVIEW_ORTHO_SCALE)
    render(camera, DNA_DETAIL_PREVIEW_PATH, (route_positions["legacy_pymol_original"], DNA_ROW_Y, 35.0), 1.4)
    render(camera, RNA_DETAIL_PREVIEW_PATH, (route_positions["legacy_pymol_original"], RNA_ROW_Y, 35.0), 1.4)
    print(f"Wrote {BLEND_PATH}")
    print(f"Wrote {PREVIEW_PATH}")
    print(f"Wrote {DNA_DETAIL_PREVIEW_PATH}")
    print(f"Wrote {RNA_DETAIL_PREVIEW_PATH}")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
