#!/usr/bin/env python3
"""Compare Blender polish passes for the preferred repaired PyMOL DNA/RNA proxy."""

from __future__ import annotations

import json
import os
from collections import Counter
from pathlib import Path

import bpy
from mathutils import Vector
from mathutils.kdtree import KDTree


ROOT = Path(os.environ.get("GENE_SCENE_ROOT", Path(__file__).resolve().parents[3])).resolve()
REDUCED_DIR = ROOT / "assets" / "pymol_exports" / "surface_assets_reduced"
OUTPUT_DIR = ROOT / "experiments" / "procedural_nucleic_acids" / "outputs"
BLEND_PATH = OUTPUT_DIR / "pymol_proxy_polish_comparison.blend"
PREVIEW_PATH = OUTPUT_DIR / "preview_pymol_proxy_polish_comparison.png"
DNA_DETAIL_PATH = OUTPUT_DIR / "preview_pymol_proxy_polish_comparison_dna_detail.png"
RNA_DETAIL_PATH = OUTPUT_DIR / "preview_pymol_proxy_polish_comparison_rna_detail.png"
REPORT_PATH = OUTPUT_DIR / "pymol_proxy_polish_comparison_report.json"

ANGSTROM_TO_MM = 0.04
PANEL_DISPLAY_SCALE = 0.045
IMPORT_SCALE = ANGSTROM_TO_MM * PANEL_DISPLAY_SCALE
BACKGROUND_COLOR = (0.985, 0.985, 0.965, 1.0)
DNA_BASE_ASSIGNMENT_MARGIN = 0.006
DNA_MATERIAL_SMOOTHING_ITERATIONS = 2

DNA_COMPONENTS = [
    ("DNA_PROXY", "strand_A", "dna_a"),
    ("DNA_PROXY", "strand_B", "dna_b"),
    ("DNA_PROXY", "base_pairs", "dna_base"),
]
RNA_COMPONENTS = [
    ("MRNA_PROXY", "utr5", "rna_utr5"),
    ("MRNA_PROXY", "coding", "rna_coding"),
    ("MRNA_PROXY", "utr3", "rna_utr3"),
]


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
    scene.render.resolution_x = 2800
    scene.render.resolution_y = 1550
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
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, ...]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("v "):
                _, x, y, z = line.split()[:4]
                vertices.append((float(x) * scale, float(y) * scale, float(z) * scale))
            elif line.startswith("f "):
                face = []
                for token in line.split()[1:]:
                    index = token.split("/")[0]
                    if index:
                        face.append(int(index) - 1)
                if len(face) >= 3:
                    faces.append(tuple(face))
    return vertices, faces


def mesh_topology_stats(vertices: list[tuple[float, float, float]], faces: list[tuple[int, ...]]) -> dict:
    edge_counts: dict[tuple[int, int], int] = {}
    for face in faces:
        for i, a in enumerate(face):
            b = face[(i + 1) % len(face)]
            edge = (a, b) if a < b else (b, a)
            edge_counts[edge] = edge_counts.get(edge, 0) + 1
    boundary_edges = sum(1 for count in edge_counts.values() if count == 1)
    non_manifold_edges = sum(1 for count in edge_counts.values() if count > 2)
    return {
        "vertices": len(vertices),
        "faces": len(faces),
        "boundary_edges": boundary_edges,
        "non_manifold_edges": non_manifold_edges,
    }


def mesh_data_from_object(obj: bpy.types.Object) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]]]:
    vertices = [(vertex.co.x, vertex.co.y, vertex.co.z) for vertex in obj.data.vertices]
    faces = [tuple(poly.vertices) for poly in obj.data.polygons]
    return vertices, faces


def object_report(objects: list[bpy.types.Object]) -> dict:
    bpy.context.view_layer.update()
    coords = []
    faces = 0
    vertices = 0
    for obj in objects:
        if obj.type != "MESH":
            continue
        faces += len(obj.data.polygons)
        vertices += len(obj.data.vertices)
        for corner in obj.bound_box:
            coords.append(obj.matrix_world @ Vector(corner))
    if not coords:
        return {"center_mm": [0.0, 0.0, 0.0], "bbox_mm": [0.0, 0.0, 0.0], "faces": faces, "vertices": vertices}
    min_v = Vector((min(c.x for c in coords), min(c.y for c in coords), min(c.z for c in coords)))
    max_v = Vector((max(c.x for c in coords), max(c.y for c in coords), max(c.z for c in coords)))
    center = (min_v + max_v) * 0.5
    return {
        "center_mm": [center.x, center.y, center.z],
        "bbox_mm": [max_v.x - min_v.x, max_v.y - min_v.y, max_v.z - min_v.z],
        "faces": faces,
        "vertices": vertices,
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


def apply_voxel_polish(
    obj: bpy.types.Object,
    voxel_size: float,
    smooth_iterations: int,
    smooth_factor: float,
) -> dict:
    before_vertices, before_faces = mesh_data_from_object(obj)
    before = mesh_topology_stats(before_vertices, before_faces)
    bpy.ops.object.select_all(action="DESELECT")
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj
    remesh = obj.modifiers.new("voxel_close_small_holes", "REMESH")
    remesh.mode = "VOXEL"
    remesh.voxel_size = voxel_size
    remesh.adaptivity = 0.0
    remesh.use_smooth_shade = True
    bpy.ops.object.modifier_apply(modifier=remesh.name)
    if smooth_iterations > 0:
        smooth = obj.modifiers.new("surface_relax", "SMOOTH")
        smooth.factor = smooth_factor
        smooth.iterations = smooth_iterations
        bpy.ops.object.modifier_apply(modifier=smooth.name)
    bpy.ops.object.shade_smooth()
    after_vertices, after_faces = mesh_data_from_object(obj)
    after = mesh_topology_stats(after_vertices, after_faces)
    return {
        "voxel_size_display_mm": voxel_size,
        "voxel_size_scene_mm": voxel_size / PANEL_DISPLAY_SCALE,
        "smooth_iterations": smooth_iterations,
        "smooth_factor": smooth_factor,
        "before": before,
        "after": after,
    }


def load_component(asset_id: str, component: str) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]], Path]:
    path = REDUCED_DIR / asset_id / f"{asset_id}_surface_{component}.obj"
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    vertices, faces = parse_obj(path, IMPORT_SCALE)
    return vertices, faces, path


def merge_components(components: list[tuple[str, str, str]]) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]], list[dict]]:
    merged_vertices: list[tuple[float, float, float]] = []
    merged_faces: list[tuple[int, ...]] = []
    component_reports = []
    for asset_id, component, _material_key in components:
        vertices, faces, path = load_component(asset_id, component)
        offset = len(merged_vertices)
        merged_vertices.extend(vertices)
        merged_faces.extend(tuple(index + offset for index in face) for face in faces)
        component_reports.append(
            {
                "asset_id": asset_id,
                "component": component,
                "source_obj": str(path.relative_to(ROOT)).replace("\\", "/"),
                **mesh_topology_stats(vertices, faces),
            }
        )
    return merged_vertices, merged_faces, component_reports


def build_reference_panel(
    components: list[tuple[str, str, str]],
    materials: dict[str, bpy.types.Material],
    collection: bpy.types.Collection,
    target_center: Vector,
    label_prefix: str,
) -> tuple[list[bpy.types.Object], list[dict]]:
    objects = []
    reports = []
    for asset_id, component, material_key in components:
        vertices, faces, path = load_component(asset_id, component)
        obj = create_mesh_object(f"{label_prefix} {asset_id} {component}", vertices, faces, materials[material_key], collection, "repaired_pymol_surface_proxy")
        objects.append(obj)
        reports.append(
            {
                "object": obj.name,
                "asset_id": asset_id,
                "component": component,
                "source_obj": str(path.relative_to(ROOT)).replace("\\", "/"),
                **mesh_topology_stats(vertices, faces),
            }
        )
    translate_objects(objects, target_center)
    return objects, reports


def build_component_polish_panel(
    components: list[tuple[str, str, str]],
    materials: dict[str, bpy.types.Material],
    collection: bpy.types.Collection,
    target_center: Vector,
    label_prefix: str,
    voxel_size: float,
) -> tuple[list[bpy.types.Object], list[dict]]:
    objects = []
    reports = []
    for asset_id, component, material_key in components:
        vertices, faces, path = load_component(asset_id, component)
        obj = create_mesh_object(f"{label_prefix} {asset_id} {component}", vertices, faces, materials[material_key], collection, "component_voxel_polished_pymol_proxy")
        polish = apply_voxel_polish(obj, voxel_size, smooth_iterations=2, smooth_factor=0.22)
        objects.append(obj)
        reports.append(
            {
                "object": obj.name,
                "asset_id": asset_id,
                "component": component,
                "source_obj": str(path.relative_to(ROOT)).replace("\\", "/"),
                "polish": polish,
            }
        )
    translate_objects(objects, target_center)
    return objects, reports


def build_unified_polish_panel(
    components: list[tuple[str, str, str]],
    mat: bpy.types.Material,
    collection: bpy.types.Collection,
    target_center: Vector,
    name: str,
    voxel_size: float,
) -> tuple[bpy.types.Object, dict]:
    vertices, faces, component_reports = merge_components(components)
    obj = create_mesh_object(name, vertices, faces, mat, collection, "unified_voxel_polished_pymol_proxy")
    polish = apply_voxel_polish(obj, voxel_size, smooth_iterations=3, smooth_factor=0.18)
    translate_objects([obj], target_center)
    return obj, {
        "object": obj.name,
        "components": component_reports,
        "polish": polish,
    }


def build_component_kdtree(components: list[tuple[str, str, str]]) -> tuple[KDTree, list[str]]:
    """Index original component vertices so a fused remesh can recover source colors."""
    vertices_by_component = []
    total_vertices = 0
    material_keys = []
    for asset_id, component, material_key in components:
        vertices, _faces, _path = load_component(asset_id, component)
        vertices_by_component.append(vertices)
        total_vertices += len(vertices)
        material_keys.append(material_key)
    kd = KDTree(total_vertices)
    index = 0
    for component_index, vertices in enumerate(vertices_by_component):
        for vertex in vertices:
            kd.insert(Vector(vertex), component_index)
            index += 1
    kd.balance()
    return kd, material_keys


def build_component_kdtrees(components: list[tuple[str, str, str]]) -> dict[str, KDTree]:
    kdtrees: dict[str, KDTree] = {}
    for asset_id, component, material_key in components:
        vertices, _faces, _path = load_component(asset_id, component)
        kd = KDTree(len(vertices))
        for index, vertex in enumerate(vertices):
            kd.insert(Vector(vertex), index)
        kd.balance()
        kdtrees[material_key] = kd
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
        if len(users) < 2:
            continue
        for index in users:
            neighbors[index].update(other for other in users if other != index)
    return neighbors


def smooth_material_assignments(mesh: bpy.types.Mesh, assignments: list[int], iterations: int) -> list[int]:
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


def transfer_component_materials(
    obj: bpy.types.Object,
    components: list[tuple[str, str, str]],
    materials: dict[str, bpy.types.Material],
) -> dict:
    kd, material_keys = build_component_kdtree(components)
    obj.data.materials.clear()
    for material_key in material_keys:
        obj.data.materials.append(materials[material_key])
    face_counts = {material_key: 0 for material_key in material_keys}
    for polygon in obj.data.polygons:
        _co, component_index, _distance = kd.find(polygon.center)
        polygon.material_index = int(component_index)
        face_counts[material_keys[int(component_index)]] += 1
    obj.data.update()
    return {
        "method": "nearest_original_component_vertex",
        "component_materials": material_keys,
        "face_counts_by_material": face_counts,
    }


def transfer_dna_materials_strand_priority(
    obj: bpy.types.Object,
    components: list[tuple[str, str, str]],
    materials: dict[str, bpy.types.Material],
) -> dict:
    material_keys = [material_key for _asset_id, _component, material_key in components]
    material_index_by_key = {material_key: index for index, material_key in enumerate(material_keys)}
    kdtrees = build_component_kdtrees(components)
    obj.data.materials.clear()
    for material_key in material_keys:
        obj.data.materials.append(materials[material_key])

    assignments: list[int] = []
    raw_counts = {material_key: 0 for material_key in material_keys}
    for polygon in obj.data.polygons:
        center = polygon.center
        strand_a_distance = nearest_distance(kdtrees["dna_a"], center)
        strand_b_distance = nearest_distance(kdtrees["dna_b"], center)
        base_distance = nearest_distance(kdtrees["dna_base"], center)
        closest_strand_key = "dna_a" if strand_a_distance <= strand_b_distance else "dna_b"
        closest_strand_distance = min(strand_a_distance, strand_b_distance)
        if base_distance + DNA_BASE_ASSIGNMENT_MARGIN < closest_strand_distance:
            material_key = "dna_base"
        else:
            material_key = closest_strand_key
        raw_counts[material_key] += 1
        assignments.append(material_index_by_key[material_key])

    smoothed_assignments = smooth_material_assignments(obj.data, assignments, DNA_MATERIAL_SMOOTHING_ITERATIONS)
    smoothed_counts = {material_key: 0 for material_key in material_keys}
    for polygon, material_index in zip(obj.data.polygons, smoothed_assignments):
        polygon.material_index = material_index
        smoothed_counts[material_keys[material_index]] += 1
    obj.data.update()
    return {
        "method": "dna_strand_priority_nearest_component_vertex",
        "component_materials": material_keys,
        "base_assignment_margin_display_mm": DNA_BASE_ASSIGNMENT_MARGIN,
        "smoothing_iterations": DNA_MATERIAL_SMOOTHING_ITERATIONS,
        "raw_face_counts_by_material": raw_counts,
        "face_counts_by_material": smoothed_counts,
    }


def build_unified_transfer_polish_panel(
    components: list[tuple[str, str, str]],
    materials: dict[str, bpy.types.Material],
    collection: bpy.types.Collection,
    target_center: Vector,
    name: str,
    voxel_size: float,
) -> tuple[bpy.types.Object, dict]:
    vertices, faces, component_reports = merge_components(components)
    obj = create_mesh_object(name, vertices, faces, materials[components[0][2]], collection, "unified_voxel_polished_with_component_material_transfer")
    polish = apply_voxel_polish(obj, voxel_size, smooth_iterations=3, smooth_factor=0.18)
    if components == DNA_COMPONENTS:
        material_transfer = transfer_dna_materials_strand_priority(obj, components, materials)
    else:
        material_transfer = transfer_component_materials(obj, components, materials)
    translate_objects([obj], target_center)
    return obj, {
        "object": obj.name,
        "components": component_reports,
        "polish": polish,
        "material_transfer": material_transfer,
    }


def add_camera() -> bpy.types.Object:
    camera_data = bpy.data.cameras.new("Camera")
    camera = bpy.data.objects.new("Camera", camera_data)
    camera.location = (0.0, 0.0, 70.0)
    camera.rotation_euler = (0.0, 0.0, 0.0)
    camera_data.type = "ORTHO"
    camera_data.ortho_scale = 26.0
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
    comparison = make_collection("PyMOL proxy polish comparison")
    labels = make_collection("Labels")
    materials = {
        "dna_a": material("dna_A_orange", (0.88, 0.48, 0.24, 1.0)),
        "dna_b": material("dna_B_brown", (0.52, 0.28, 0.15, 1.0)),
        "dna_base": material("dna_base_gold", (0.76, 0.55, 0.28, 1.0)),
        "dna_unified": material("dna_unified_warm_brown", (0.72, 0.49, 0.31, 1.0)),
        "rna_utr5": material("rna_utr5_olive", (0.60, 0.58, 0.25, 1.0)),
        "rna_coding": material("rna_coding_orange", (0.84, 0.45, 0.21, 1.0)),
        "rna_utr3": material("rna_utr3_yellow_olive", (0.72, 0.66, 0.32, 1.0)),
        "rna_unified": material("rna_unified_olive", (0.62, 0.61, 0.38, 1.0)),
        "text": material("text_grey", (0.20, 0.20, 0.20, 1.0)),
    }
    x_positions = {
        "reference": -18.0,
        "component_polished": -6.0,
        "unified_polished": 6.0,
        "unified_transfer": 18.0,
    }
    dna_y = 3.6
    rna_y = -3.6
    report = {
        "title": "PyMOL proxy Blender polish comparison",
        "canonical_scene_changed": False,
        "units": {
            "angstrom_to_mm": ANGSTROM_TO_MM,
            "panel_display_scale": PANEL_DISPLAY_SCALE,
            "import_scale_mm_per_A": IMPORT_SCALE,
        },
        "routes": {},
        "outputs": {
            "blend": str(BLEND_PATH.relative_to(ROOT)).replace("\\", "/"),
            "preview": str(PREVIEW_PATH.relative_to(ROOT)).replace("\\", "/"),
            "dna_detail_preview": str(DNA_DETAIL_PATH.relative_to(ROOT)).replace("\\", "/"),
            "rna_detail_preview": str(RNA_DETAIL_PATH.relative_to(ROOT)).replace("\\", "/"),
        },
        "viability_note": (
            "Component polish preserves DNA/RNA subcomponent colors. Unified polish closes small holes "
            "and fuses transitions more cleanly. Unified polish with component material transfer keeps "
            "the fused/watertight surface while recoloring faces by nearest source component."
        ),
    }
    report["routes"]["repaired_reference"] = {
        "dna": build_reference_panel(DNA_COMPONENTS, materials, comparison, Vector((x_positions["reference"], dna_y, 0.0)), "reference")[1],
        "rna": build_reference_panel(RNA_COMPONENTS, materials, comparison, Vector((x_positions["reference"], rna_y, 0.0)), "reference")[1],
    }
    report["routes"]["component_voxel_polished"] = {
        "dna": build_component_polish_panel(DNA_COMPONENTS, materials, comparison, Vector((x_positions["component_polished"], dna_y, 0.0)), "component polished", voxel_size=0.0045)[1],
        "rna": build_component_polish_panel(RNA_COMPONENTS, materials, comparison, Vector((x_positions["component_polished"], rna_y, 0.0)), "component polished", voxel_size=0.0045)[1],
    }
    dna_unified_obj, dna_unified_report = build_unified_polish_panel(
        DNA_COMPONENTS,
        materials["dna_unified"],
        comparison,
        Vector((x_positions["unified_polished"], dna_y, 0.0)),
        "unified polished DNA_PROXY",
        voxel_size=0.0048,
    )
    rna_unified_obj, rna_unified_report = build_unified_polish_panel(
        RNA_COMPONENTS,
        materials["rna_unified"],
        comparison,
        Vector((x_positions["unified_polished"], rna_y, 0.0)),
        "unified polished MRNA_PROXY",
        voxel_size=0.0048,
    )
    report["routes"]["unified_voxel_polished"] = {"dna": dna_unified_report, "rna": rna_unified_report}
    dna_transfer_obj, dna_transfer_report = build_unified_transfer_polish_panel(
        DNA_COMPONENTS,
        materials,
        comparison,
        Vector((x_positions["unified_transfer"], dna_y, 0.0)),
        "unified material-transfer polished DNA_PROXY",
        voxel_size=0.0048,
    )
    rna_transfer_obj, rna_transfer_report = build_unified_transfer_polish_panel(
        RNA_COMPONENTS,
        materials,
        comparison,
        Vector((x_positions["unified_transfer"], rna_y, 0.0)),
        "unified material-transfer polished MRNA_PROXY",
        voxel_size=0.0048,
    )
    report["routes"]["unified_voxel_polished_material_transfer"] = {"dna": dna_transfer_report, "rna": rna_transfer_report}

    add_text("label_reference", "repaired PyMOL proxy", (x_positions["reference"], 8.6, 0.4), 0.45, materials["text"], labels)
    add_text("label_component", "component voxel polish", (x_positions["component_polished"], 8.6, 0.4), 0.45, materials["text"], labels)
    add_text("label_unified", "unified voxel polish", (x_positions["unified_polished"], 8.6, 0.4), 0.45, materials["text"], labels)
    add_text("label_unified_transfer", "unified + material transfer", (x_positions["unified_transfer"], 8.6, 0.4), 0.45, materials["text"], labels)
    add_text("label_scale", "All panels use repaired PyMOL proxy assets at 1 A = 0.04 mm, displayed at 4.5%", (0.0, -8.7, 0.4), 0.34, materials["text"], labels)

    camera = add_camera()
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_PATH))
    render(camera, PREVIEW_PATH, (0.0, 0.0, 72.0), 36.0)
    render(camera, DNA_DETAIL_PATH, (x_positions["unified_transfer"], dna_y, 35.0), 1.7)
    render(camera, RNA_DETAIL_PATH, (x_positions["unified_transfer"], rna_y, 35.0), 1.7)
    print(f"Wrote {BLEND_PATH}")
    print(f"Wrote {PREVIEW_PATH}")
    print(f"Wrote {DNA_DETAIL_PATH}")
    print(f"Wrote {RNA_DETAIL_PATH}")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
