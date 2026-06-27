#!/usr/bin/env python3
"""Build the V3 gene-expression scene with PyMOL PDB surfaces and polished procedural nucleic-acid proxies."""

from __future__ import annotations

import json
import math
import sys
from collections import Counter
from pathlib import Path

import bpy
from mathutils import Euler, Matrix, Vector
from mathutils.kdtree import KDTree


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import blender_nucleic_meshes as direct_nucleic_meshes  # noqa: E402
import build_gene_expression_scene as base  # noqa: E402


RAW_SURFACE_DIR = ROOT / "assets" / "pymol_exports" / "surface_assets"
REDUCED_SURFACE_DIR = ROOT / "assets" / "pymol_exports" / "surface_assets_reduced"
REDUCED_MANIFEST_PATH = REDUCED_SURFACE_DIR / "surface_assets_reduced_manifest.json"
OUTPUT_DIR = ROOT / "outputs" / "canonical"
BLEND_PATH = OUTPUT_DIR / "gene_expression_surface_style_v3.blend"
PREVIEW_PATH = OUTPUT_DIR / "preview_gene_expression_surface_style_v3.png"
REPORT_PATH = OUTPUT_DIR / "gene_expression_surface_scene_v3_report.json"
DETAIL_PREVIEWS = {
    "dna_transcription": OUTPUT_DIR / "preview_gene_expression_surface_style_v3_dna_transcription_detail.png",
    "cas9_binding": OUTPUT_DIR / "preview_gene_expression_surface_style_v3_cas9_binding_detail.png",
    "nucleosome_binding": OUTPUT_DIR / "preview_gene_expression_surface_style_v3_nucleosome_binding_detail.png",
    "translation": OUTPUT_DIR / "preview_gene_expression_surface_style_v3_translation_detail.png",
    "ribosome": OUTPUT_DIR / "preview_gene_expression_surface_style_v3_ribosome_detail.png",
}

BACKGROUND_COLOR = (0.985, 0.985, 0.965, 1.0)
RIBOSOME_PDBS = {"1J5E", "1JJ2"}
NUCLEIC_PROXY_STYLE = "procedural_pymol_surface_proxy_polished"
DIRECT_NUCLEIC_STYLE = "direct_blender_surface_proxy_polished"
NUCLEIC_PROXY_VOXEL_SIZE_MM = 0.075
NUCLEIC_PROXY_SMOOTH_ITERATIONS = 1
NUCLEIC_PROXY_SMOOTH_FACTOR = 0.10
DNA_BASE_STRICT_MARGIN_MM = 0.040
DNA_BASE_SOFT_MARGIN_MM = 0.060
DNA_STRAND_SYMMETRY_MM = 0.22
DNA_MATERIAL_SMOOTHING_ITERATIONS = 2
DNA_PROXY_COMPONENTS = [
    ("strand_A", "dna_orange"),
    ("strand_B", "dna_dark"),
    ("base_pairs", "rna_gold"),
]
MRNA_PROXY_COMPONENTS = [
    ("utr5", "olive"),
    ("coding", "orange"),
    ("utr3", "yellow_olive"),
]
DNA_BINDING_COLLECTIONS = {"DNA", "Transcription"}
DNA_ALIGNMENT_GUIDE_COMPONENTS = {"nucleic"}


class SampledPath:
    def __init__(self, points: list[tuple[float, float, float]]):
        self.points = [Vector(point) for point in points]
        self.cumulative = [0.0]
        for i in range(1, len(self.points)):
            self.cumulative.append(self.cumulative[-1] + (self.points[i] - self.points[i - 1]).length)
        self.length = self.cumulative[-1] if self.cumulative else 0.0

    def point_at_length(self, distance: float) -> Vector:
        if not self.points:
            return Vector((0.0, 0.0, 0.0))
        distance = min(max(distance, 0.0), self.length)
        for i in range(1, len(self.cumulative)):
            if self.cumulative[i] >= distance:
                start = self.points[i - 1]
                end = self.points[i]
                span = self.cumulative[i] - self.cumulative[i - 1]
                t = 0.0 if span == 0 else (distance - self.cumulative[i - 1]) / span
                return start.lerp(end, t)
        return self.points[-1]

    def tangent_at_length(self, distance: float) -> Vector:
        if len(self.points) < 2:
            return Vector((1.0, 0.0, 0.0))
        distance = min(max(distance, 0.0), self.length)
        for i in range(1, len(self.cumulative)):
            if self.cumulative[i] >= distance:
                tangent = self.points[i] - self.points[i - 1]
                return tangent.normalized() if tangent.length else Vector((1.0, 0.0, 0.0))
        tangent = self.points[-1] - self.points[-2]
        return tangent.normalized() if tangent.length else Vector((1.0, 0.0, 0.0))

    def closest_to_xy(self, location: tuple[float, float, float] | list[float]) -> tuple[Vector, Vector, float]:
        target = Vector((location[0], location[1], 0.0))
        best_i = 0
        best_distance = float("inf")
        for i, point in enumerate(self.points):
            distance = (Vector((point.x, point.y, 0.0)) - target).length
            if distance < best_distance:
                best_i = i
                best_distance = distance
        length = self.cumulative[best_i]
        return self.points[best_i], self.tangent_at_length(length), length

    def closest_to_x(self, x: float) -> tuple[Vector, Vector, float]:
        best_i = min(range(len(self.points)), key=lambda i: abs(self.points[i].x - x))
        length = self.cumulative[best_i]
        return self.points[best_i], self.tangent_at_length(length), length


def catmull_rom(control_points: list[tuple[float, float, float]], samples_per_segment: int = 24) -> list[tuple[float, float, float]]:
    points = [Vector(point) for point in control_points]
    if len(points) < 2:
        return control_points
    result = []
    for i in range(len(points) - 1):
        p0 = points[max(0, i - 1)]
        p1 = points[i]
        p2 = points[i + 1]
        p3 = points[min(len(points) - 1, i + 2)]
        for j in range(samples_per_segment):
            t = j / samples_per_segment
            t2 = t * t
            t3 = t2 * t
            value = 0.5 * (
                (2.0 * p1)
                + (-p0 + p2) * t
                + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2
                + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3
            )
            result.append((value.x, value.y, value.z))
    result.append(tuple(points[-1]))
    return result


def clean_scene() -> None:
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete()
    for collection in list(bpy.data.collections):
        bpy.data.collections.remove(collection)
    for mesh in list(bpy.data.meshes):
        if mesh.users == 0:
            bpy.data.meshes.remove(mesh)
    for curve in list(bpy.data.curves):
        if curve.users == 0:
            bpy.data.curves.remove(curve)


def configure_scene() -> None:
    scene = bpy.context.scene
    scene.unit_settings.system = "METRIC"
    scene.unit_settings.length_unit = "MILLIMETERS"
    scene.unit_settings.scale_length = 0.001
    scene.render.engine = "BLENDER_WORKBENCH"
    scene.render.resolution_x = 3000
    scene.render.resolution_y = 1800
    scene.display.shading.light = "STUDIO"
    scene.display.shading.color_type = "MATERIAL"
    scene.display.shading.background_type = "VIEWPORT"
    scene.display.shading.background_color = BACKGROUND_COLOR[:3]
    scene.view_settings.view_transform = "Standard"
    scene.view_settings.look = "None"
    scene.view_settings.exposure = 0.0
    scene.view_settings.gamma = 1.0

    world = scene.world or bpy.data.worlds.new("World")
    scene.world = world
    world.color = BACKGROUND_COLOR[:3]
    world.use_nodes = True
    background = world.node_tree.nodes.get("Background")
    if background:
        background.inputs["Color"].default_value = BACKGROUND_COLOR
        background.inputs["Strength"].default_value = 1.0


def soften_materials(materials: dict[str, bpy.types.Material]) -> None:
    for material in materials.values():
        material.diffuse_color = tuple(
            min(1.0, channel * 0.94 + 0.06) if i < 3 else channel
            for i, channel in enumerate(material.diffuse_color)
        )
        bsdf = material.node_tree.nodes.get("Principled BSDF") if material.use_nodes else None
        if bsdf:
            bsdf.inputs["Base Color"].default_value = material.diffuse_color
            bsdf.inputs["Roughness"].default_value = 0.9
            bsdf.inputs["Metallic"].default_value = 0.0


def path_frame(tangent: Vector) -> tuple[Vector, Vector]:
    tangent = tangent.normalized() if tangent.length else Vector((1.0, 0.0, 0.0))
    up = Vector((0.0, 0.0, 1.0))
    normal = up.cross(tangent)
    if normal.length < 1e-6:
        up = Vector((0.0, 1.0, 0.0))
        normal = up.cross(tangent)
    normal = normal.normalized()
    binormal = tangent.cross(normal).normalized()
    return normal, binormal


def is_dna_binding_asset(asset: dict) -> bool:
    return asset.get("anchor_kind") == "nucleic" and asset.get("collection") in DNA_BINDING_COLLECTIONS


def is_wrapped_dna_asset(asset: dict) -> bool:
    return asset.get("dna_binding_mode") == "wrapped_loop" or asset.get("name", "").lower() == "nucleosome"


def uses_surface_alignment(asset: dict) -> bool:
    path_anchor = asset.get("path_anchor") or {}
    binding_mode = str(path_anchor.get("binding_mode", ""))
    return (
        is_dna_binding_asset(asset)
        or binding_mode.startswith("co_crystal_rna")
        or binding_mode == "protein_surface_rna_contact"
    )


def project_perpendicular(vector: Vector, axis: Vector) -> Vector:
    axis = axis.normalized() if axis.length else Vector((1.0, 0.0, 0.0))
    return vector - axis * vector.dot(axis)


def normalized_or_none(vector: Vector) -> Vector | None:
    return vector.normalized() if vector.length > 1e-7 else None


def covariance_axes(vectors: list[Vector]) -> list[tuple[float, Vector]]:
    center = base.mean_vector(vectors)
    cov = [[0.0, 0.0, 0.0] for _ in range(3)]
    for vector in vectors:
        delta = vector - center
        values = (delta.x, delta.y, delta.z)
        for i in range(3):
            for j in range(3):
                cov[i][j] += values[i] * values[j]
    if vectors:
        inv_n = 1.0 / len(vectors)
        for i in range(3):
            for j in range(3):
                cov[i][j] *= inv_n

    eigenvectors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    for _ in range(36):
        p, q = 0, 1
        max_offdiag = abs(cov[0][1])
        for i, j in ((0, 2), (1, 2)):
            if abs(cov[i][j]) > max_offdiag:
                p, q = i, j
                max_offdiag = abs(cov[i][j])
        if max_offdiag < 1e-10:
            break
        app = cov[p][p]
        aqq = cov[q][q]
        apq = cov[p][q]
        angle = 0.5 * math.atan2(2.0 * apq, aqq - app)
        c = math.cos(angle)
        s = math.sin(angle)
        for k in range(3):
            if k in {p, q}:
                continue
            akp = cov[k][p]
            akq = cov[k][q]
            cov[k][p] = cov[p][k] = c * akp - s * akq
            cov[k][q] = cov[q][k] = s * akp + c * akq
        cov[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq
        cov[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq
        cov[p][q] = cov[q][p] = 0.0
        for k in range(3):
            vip = eigenvectors[k][p]
            viq = eigenvectors[k][q]
            eigenvectors[k][p] = c * vip - s * viq
            eigenvectors[k][q] = s * vip + c * viq

    axes = []
    for i in range(3):
        axis = Vector((eigenvectors[0][i], eigenvectors[1][i], eigenvectors[2][i]))
        axes.append((cov[i][i], axis.normalized() if axis.length else Vector((1.0, 0.0, 0.0))))
    return sorted(axes, key=lambda item: item[0], reverse=True)


def axis_roll_rotation(axis: Vector, current: Vector, target: Vector) -> tuple[Matrix, float | None]:
    axis = axis.normalized() if axis.length else Vector((1.0, 0.0, 0.0))
    current_perp = normalized_or_none(project_perpendicular(current, axis))
    target_perp = normalized_or_none(project_perpendicular(target, axis))
    if current_perp is None or target_perp is None:
        return Matrix.Identity(3), None
    cross = current_perp.cross(target_perp)
    angle = math.atan2(axis.dot(cross), max(-1.0, min(1.0, current_perp.dot(target_perp))))
    return Matrix.Rotation(angle, 3, axis), math.degrees(angle)


def target_binding_side(asset: dict, target_axis: Vector) -> Vector:
    preferred = Vector(asset.get("binding_side_mm", [0.0, 1.0, 0.0]))
    side = normalized_or_none(project_perpendicular(preferred, target_axis))
    if side is not None:
        return side
    normal, _binormal = path_frame(target_axis)
    return normal


def vector_sample(vectors: list[Vector], max_count: int) -> list[Vector]:
    if len(vectors) <= max_count:
        return vectors
    stride = max(1, math.ceil(len(vectors) / max_count))
    return vectors[::stride][:max_count]


def closest_surface_pair_A(protein_vectors: list[Vector], nucleic_vectors: list[Vector]) -> tuple[Vector | None, Vector | None, float | None]:
    if not protein_vectors or not nucleic_vectors:
        return None, None, None
    sampled_protein = vector_sample(protein_vectors, 30000)
    tree = KDTree(len(sampled_protein))
    for index, vertex in enumerate(sampled_protein):
        tree.insert(vertex, index)
    tree.balance()
    best_protein = None
    best_nucleic = None
    best_distance = float("inf")
    for nucleic_vertex in vector_sample(nucleic_vectors, 30000):
        protein_vertex, _index, distance = tree.find(nucleic_vertex)
        if distance < best_distance:
            best_distance = distance
            best_protein = Vector(protein_vertex)
            best_nucleic = Vector(nucleic_vertex)
    if best_protein is None or best_nucleic is None:
        return None, None, None
    return best_protein, best_nucleic, best_distance


def path_contact_radius_mm(path_name: str | None, path_context: dict) -> float:
    if path_name == "mrna":
        return float(path_context.get("mrna_contact_radius_mm", 0.30))
    if path_name == "dna":
        return float(path_context.get("dna_contact_radius_mm", 0.44))
    return 0.0


def dna_loop_center(path_context: dict) -> Vector | None:
    features = path_context.get("dna_features", {})
    loop = features.get("nucleosome_loop") if isinstance(features, dict) else None
    if not loop:
        return None
    return Vector(loop.get("center_mm", (72.0, -45.0, 0.0)))


def dna_loop_feature(path_context: dict) -> dict | None:
    features = path_context.get("dna_features", {})
    loop = features.get("nucleosome_loop") if isinstance(features, dict) else None
    return loop if isinstance(loop, dict) else None


def matrix_from_basis_columns(basis: list[list[float]]) -> Matrix:
    x_axis = Vector(basis[0])
    y_axis = Vector(basis[1])
    z_axis = Vector(basis[2])
    return Matrix(
        (
            (x_axis.x, y_axis.x, z_axis.x),
            (x_axis.y, y_axis.y, z_axis.y),
            (x_axis.z, y_axis.z, z_axis.z),
        )
    )


def rotation_from_basis_report(loop: dict) -> Matrix:
    source = matrix_from_basis_columns(loop["source_basis"])
    target = matrix_from_basis_columns(loop["target_basis"])
    return target @ source.transposed()


def surface_export_origin_A(source_vectors: list[Vector], exported_vectors: list[Vector]) -> Vector:
    if not source_vectors or not exported_vectors:
        return Vector((0.0, 0.0, 0.0))
    return base.bbox_center(source_vectors) - base.bbox_center(exported_vectors)


def build_procedural_dna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = direct_nucleic_meshes.build_dna_meshes(manifest, collections, materials)
    base.create_text("label_DNA", "DNA", (-108.0, -30.8, 0.2), 4.0, materials["black"], collections["Labels"])
    return {"path": path, "report": report, "features": {"nucleosome_loop": report.get("nucleosome_loop")}}


def build_procedural_mrna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = direct_nucleic_meshes.build_mrna_meshes(manifest, collections, materials)
    base.create_text("label_mRNA_total", "actin mRNA: 1852 nt", (-45.0, 62.5, 0.2), 4.0, materials["black"], collections["Labels"])
    for segment in manifest["mrna"]["segments"]:
        label_location = segment.get("label_location_mm", (-40.0, 52.0, 0.3))
        base.create_text(
            f"label_{segment['name']}",
            f"{segment['name']} ({segment['nt']} nt)",
            tuple(label_location),
            2.3,
            materials["label_grey"],
            collections["Labels"],
        )
    return {"path": path, "report": report}


def build_compact_mrna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = direct_nucleic_meshes.build_compact_mrna_meshes(manifest, collections, materials)
    compact_center = path_center(path)
    base.create_text(
        "label_mRNA_compact",
        "compact mRNP-like mRNA",
        (compact_center.x, compact_center.y + 8.0, compact_center.z + 0.4),
        2.3,
        materials["label_grey"],
        collections["Labels"],
    )
    return {"path": path, "report": report}


def mesh_bounds(vertices: list[tuple[float, float, float]]) -> dict[str, list[float]]:
    min_v = Vector((min(v[0] for v in vertices), min(v[1] for v in vertices), min(v[2] for v in vertices)))
    max_v = Vector((max(v[0] for v in vertices), max(v[1] for v in vertices), max(v[2] for v in vertices)))
    return {
        "min_mm": [min_v.x, min_v.y, min_v.z],
        "max_mm": [max_v.x, max_v.y, max_v.z],
        "bbox_mm": [max_v.x - min_v.x, max_v.y - min_v.y, max_v.z - min_v.z],
    }


def bounds_center(bounds: dict[str, list[float]]) -> Vector:
    min_v = bounds["min_mm"]
    max_v = bounds["max_mm"]
    return Vector(((min_v[0] + max_v[0]) * 0.5, (min_v[1] + max_v[1]) * 0.5, (min_v[2] + max_v[2]) * 0.5))


def path_center(path) -> Vector:
    points = path.points
    min_v = Vector((min(point.x for point in points), min(point.y for point in points), min(point.z for point in points)))
    max_v = Vector((max(point.x for point in points), max(point.y for point in points), max(point.z for point in points)))
    return (min_v + max_v) * 0.5


def translate_mesh_vertices(obj: bpy.types.Object, delta: Vector) -> None:
    for vertex in obj.data.vertices:
        vertex.co += delta
    obj.data.update()


def mesh_topology_stats(vertices: list[tuple[float, float, float]], faces: list[tuple[int, ...]]) -> dict:
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


def mesh_data_from_object(obj: bpy.types.Object) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]]]:
    vertices = [(vertex.co.x, vertex.co.y, vertex.co.z) for vertex in obj.data.vertices]
    faces = [tuple(poly.vertices) for poly in obj.data.polygons]
    return vertices, faces


def apply_voxel_polish(obj: bpy.types.Object) -> dict:
    before_vertices, before_faces = mesh_data_from_object(obj)
    before = mesh_topology_stats(before_vertices, before_faces)
    bpy.ops.object.select_all(action="DESELECT")
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj
    remesh = obj.modifiers.new("voxel_close_small_holes", "REMESH")
    remesh.mode = "VOXEL"
    remesh.voxel_size = NUCLEIC_PROXY_VOXEL_SIZE_MM
    remesh.adaptivity = 0.0
    remesh.use_smooth_shade = True
    bpy.ops.object.modifier_apply(modifier=remesh.name)
    smooth = obj.modifiers.new("surface_relax", "SMOOTH")
    smooth.factor = NUCLEIC_PROXY_SMOOTH_FACTOR
    smooth.iterations = NUCLEIC_PROXY_SMOOTH_ITERATIONS
    bpy.ops.object.modifier_apply(modifier=smooth.name)
    bpy.ops.object.shade_smooth()
    after_vertices, after_faces = mesh_data_from_object(obj)
    return {
        "voxel_size_mm": NUCLEIC_PROXY_VOXEL_SIZE_MM,
        "smooth_iterations": NUCLEIC_PROXY_SMOOTH_ITERATIONS,
        "smooth_factor": NUCLEIC_PROXY_SMOOTH_FACTOR,
        "before": before,
        "after": mesh_topology_stats(after_vertices, after_faces),
    }


def proxy_component_path(asset_id: str, component: str) -> Path:
    path = REDUCED_SURFACE_DIR / asset_id / f"{asset_id}_surface_{component}.obj"
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing reduced procedural PyMOL proxy surface: {path}")
    return path


def merge_proxy_components(
    asset_id: str,
    components: list[tuple[str, str]],
    scale: float,
) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]], list[dict], dict[str, list[tuple[float, float, float]]]]:
    merged_vertices: list[tuple[float, float, float]] = []
    merged_faces: list[tuple[int, ...]] = []
    reports = []
    component_vertices: dict[str, list[tuple[float, float, float]]] = {}
    transform = lambda vertex_A: vertex_A * scale
    for component, material_key in components:
        path = proxy_component_path(asset_id, component)
        vertices, faces = parse_obj_transformed(path, transform)
        offset = len(merged_vertices)
        merged_vertices.extend(vertices)
        merged_faces.extend(tuple(index + offset for index in face) for face in faces)
        component_vertices[material_key] = vertices
        reports.append(
            {
                "component": component,
                "material": material_key,
                "source_obj": str(path.relative_to(ROOT)).replace("\\", "/"),
                **mesh_topology_stats(vertices, faces),
            }
        )
    return merged_vertices, merged_faces, reports, component_vertices


def create_proxy_mesh_object(
    name: str,
    vertices: list[tuple[float, float, float]],
    faces: list[tuple[int, ...]],
    material: bpy.types.Material,
    collection: bpy.types.Collection,
    asset_id: str,
) -> bpy.types.Object:
    mesh = bpy.data.meshes.new(f"{name}_mesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    for polygon in mesh.polygons:
        polygon.use_smooth = True
    mesh.materials.append(material)
    obj = bpy.data.objects.new(name, mesh)
    base.link_to_collection(obj, collection)
    obj["style"] = NUCLEIC_PROXY_STYLE
    obj["asset_id"] = asset_id
    obj["source"] = "procedural_pseudoatom_pymol_surface"
    obj["coordinate_units"] = "millimeter"
    obj["source_coordinate_units"] = "angstrom"
    obj["angstrom_to_mm"] = 0.04
    return obj


def build_component_kdtrees(component_vertices: dict[str, list[tuple[float, float, float]]]) -> dict[str, KDTree]:
    kdtrees: dict[str, KDTree] = {}
    for material_key, vertices in component_vertices.items():
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


def transfer_rna_materials(
    obj: bpy.types.Object,
    components: list[tuple[str, str]],
    component_vertices: dict[str, list[tuple[float, float, float]]],
    materials: dict[str, bpy.types.Material],
) -> dict:
    material_keys = [material_key for _component, material_key in components]
    kdtrees = build_component_kdtrees(component_vertices)
    obj.data.materials.clear()
    for material_key in material_keys:
        obj.data.materials.append(materials[material_key])
    face_counts = {material_key: 0 for material_key in material_keys}
    for polygon in obj.data.polygons:
        distances = [(nearest_distance(kdtrees[material_key], polygon.center), material_key) for material_key in material_keys]
        _distance, material_key = min(distances, key=lambda item: item[0])
        polygon.material_index = material_keys.index(material_key)
        face_counts[material_key] += 1
    obj.data.update()
    return {"method": "nearest_original_component_vertex", "face_counts_by_material": face_counts}


def transfer_dna_materials(
    obj: bpy.types.Object,
    components: list[tuple[str, str]],
    component_vertices: dict[str, list[tuple[float, float, float]]],
    materials: dict[str, bpy.types.Material],
) -> dict:
    material_keys = [material_key for _component, material_key in components]
    material_index_by_key = {material_key: index for index, material_key in enumerate(material_keys)}
    kdtrees = build_component_kdtrees(component_vertices)
    obj.data.materials.clear()
    for material_key in material_keys:
        obj.data.materials.append(materials[material_key])

    assignments: list[int] = []
    raw_counts = {material_key: 0 for material_key in material_keys}
    for polygon in obj.data.polygons:
        center = polygon.center
        strand_a_distance = nearest_distance(kdtrees["dna_orange"], center)
        strand_b_distance = nearest_distance(kdtrees["dna_dark"], center)
        base_distance = nearest_distance(kdtrees["rna_gold"], center)
        closest_strand_key = "dna_orange" if strand_a_distance <= strand_b_distance else "dna_dark"
        closest_strand_distance = min(strand_a_distance, strand_b_distance)
        strand_symmetry = abs(strand_a_distance - strand_b_distance)
        base_between_strands = strand_symmetry <= DNA_STRAND_SYMMETRY_MM
        base_strict = base_distance + DNA_BASE_STRICT_MARGIN_MM < closest_strand_distance
        base_soft = base_distance < closest_strand_distance + DNA_BASE_SOFT_MARGIN_MM
        if base_between_strands and (base_strict or base_soft):
            material_key = "rna_gold"
        else:
            material_key = closest_strand_key
        raw_counts[material_key] += 1
        assignments.append(material_index_by_key[material_key])

    smoothed = smooth_material_assignments(obj.data, assignments, DNA_MATERIAL_SMOOTHING_ITERATIONS)
    face_counts = {material_key: 0 for material_key in material_keys}
    for polygon, material_index in zip(obj.data.polygons, smoothed):
        polygon.material_index = material_index
        face_counts[material_keys[material_index]] += 1
    obj.data.update()
    return {
        "method": "dna_strand_priority_nearest_component_vertex",
        "base_strict_margin_mm": DNA_BASE_STRICT_MARGIN_MM,
        "base_soft_margin_mm": DNA_BASE_SOFT_MARGIN_MM,
        "strand_symmetry_mm": DNA_STRAND_SYMMETRY_MM,
        "smoothing_iterations": DNA_MATERIAL_SMOOTHING_ITERATIONS,
        "raw_face_counts_by_material": raw_counts,
        "face_counts_by_material": face_counts,
    }


def build_polished_nucleic_proxy(
    asset_id: str,
    name: str,
    components: list[tuple[str, str]],
    materials: dict[str, bpy.types.Material],
    collection: bpy.types.Collection,
    transfer_kind: str,
    target_center: Vector,
) -> dict:
    scale = 0.04
    vertices, faces, component_reports, component_vertices = merge_proxy_components(asset_id, components, scale)
    obj = create_proxy_mesh_object(name, vertices, faces, materials[components[0][1]], collection, asset_id)
    polish_report = apply_voxel_polish(obj)
    if transfer_kind == "dna":
        material_report = transfer_dna_materials(obj, components, component_vertices, materials)
    else:
        material_report = transfer_rna_materials(obj, components, component_vertices, materials)
    pre_alignment_vertices, _pre_alignment_faces = mesh_data_from_object(obj)
    pre_alignment_bounds = mesh_bounds(pre_alignment_vertices)
    translate_mesh_vertices(obj, target_center - bounds_center(pre_alignment_bounds))
    final_vertices, final_faces = mesh_data_from_object(obj)
    final_bounds = mesh_bounds(final_vertices)
    return {
        "object": obj.name,
        "asset_id": asset_id,
        "style": NUCLEIC_PROXY_STYLE,
        "source_coordinate_units": "angstrom",
        "angstrom_to_mm": scale,
        "components": component_reports,
        "polish": polish_report,
        "material_transfer": material_report,
        "path_alignment": {
            "source_center_before_alignment_mm": list(bounds_center(pre_alignment_bounds)),
            "target_center_mm": [target_center.x, target_center.y, target_center.z],
        },
        "final_mesh": {**final_bounds, **mesh_topology_stats(final_vertices, final_faces)},
    }


def points_for_asset(asset: dict) -> tuple[list[dict], list[dict]]:
    pdb_id = asset["pdb_id"].upper()
    atoms = base.parse_atom_site(base.ASSET_DIR / f"{pdb_id}.cif")
    points = base.residue_points(atoms)
    return atoms, points


def parse_obj_vertices_A(path: Path) -> list[Vector]:
    vertices = []
    if not path.exists() or path.stat().st_size == 0:
        return vertices
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            if raw_line.startswith("v "):
                _, x, y, z = raw_line.split()[:4]
                vertices.append(Vector((float(x), float(y), float(z))))
    return vertices


def surface_component_path_for_alignment(pdb_id: str, component: str) -> Path | None:
    directory = surface_directory_for(pdb_id)
    path = directory / f"{pdb_id}_surface_{component}.obj"
    return path if path.exists() and path.stat().st_size > 0 else None


def surface_alignment_vectors(asset: dict) -> dict[str, list[Vector]]:
    if not uses_surface_alignment(asset):
        return {}
    pdb_id = asset["pdb_id"].upper()
    vectors = {}
    for component in ("nucleic", "protein"):
        path = surface_component_path_for_alignment(pdb_id, component)
        if path is not None:
            vertices = parse_obj_vertices_A(path)
            if vertices:
                vectors[component] = vertices
    return vectors


def path_adjusted_location_and_axis(asset: dict, path_context: dict[str, SampledPath]) -> tuple[Vector, Vector | None, float | None, str, dict | None]:
    original = Vector(asset["location_mm"])
    explicit_anchor = asset.get("path_anchor")
    if explicit_anchor:
        path_name = explicit_anchor.get("path")
        path = path_context.get(path_name)
        if path is None:
            raise KeyError(f"Unknown path anchor {path_name!r} for {asset.get('name')}")
        if "distance_mm" in explicit_anchor:
            length = max(0.0, min(path.length, float(explicit_anchor["distance_mm"])))
        else:
            length = max(0.0, min(path.length, float(explicit_anchor.get("fraction", 0.0)) * path.length))
        point = path.point_at_length(length)
        tangent = path.tangent_at_length(length)
        point = Vector((point.x, point.y, point.z))
        tangent = Vector((tangent.x, tangent.y, tangent.z))
        anchor_point = point
        if is_wrapped_dna_asset(asset) and path_name == "dna":
            loop_center = dna_loop_center(path_context)
            if loop_center is not None:
                anchor_point = Vector((loop_center.x, loop_center.y, loop_center.z))
        world_offset = Vector(explicit_anchor.get("offset_mm", (0.0, 0.0, 0.0)))
        local_offset = Vector(explicit_anchor.get("offset_local_mm", (0.0, 0.0, 0.0)))
        normal, binormal = path_frame(tangent)
        local_world_offset = normal * local_offset.x + binormal * local_offset.y + tangent * local_offset.z
        location = anchor_point + world_offset + local_world_offset
        attachment = {
            "path": path_name,
            "distance_mm": length,
            "path_length_mm": path.length,
            "fraction": length / path.length if path.length else None,
            "path_point_mm": [point.x, point.y, point.z],
            "anchor_point_mm": [anchor_point.x, anchor_point.y, anchor_point.z],
            "path_tangent": [tangent.x, tangent.y, tangent.z],
            "path_normal": [normal.x, normal.y, normal.z],
            "path_binormal": [binormal.x, binormal.y, binormal.z],
            "offset_mm": list(explicit_anchor.get("offset_mm", (0.0, 0.0, 0.0))),
            "offset_local_mm": list(explicit_anchor.get("offset_local_mm", (0.0, 0.0, 0.0))),
            "roll_deg": float(explicit_anchor.get("roll_deg", 0.0)),
            "binding_mode": explicit_anchor.get("binding_mode", "manual_tangent_frame"),
        }
        return location, tangent, length, f"{path_name}_path_explicit_anchor", attachment
    if is_wrapped_dna_asset(asset):
        loop_center = dna_loop_center(path_context)
        if loop_center is not None:
            return Vector((loop_center.x, loop_center.y, original.z)), Vector((0.0, 0.0, 1.0)), None, "dna_wrapped_loop_center", None
    if asset.get("anchor_kind") == "nucleic" and asset["collection"] in {"DNA", "Transcription"}:
        if path_context.get("dna_anchor_mode") == "closest_xy":
            point, tangent, length = path_context["dna"].closest_to_xy(original)
        else:
            point, tangent, length = path_context["dna"].closest_to_x(original.x)
        return Vector((point.x, point.y, original.z)), Vector((tangent.x, tangent.y, tangent.z)), length, "dna_path_closest_x", None
    if asset.get("anchor_kind") == "nucleic" and asset["collection"] == "RBPs":
        point, tangent, length = path_context["mrna"].closest_to_xy(original)
        return Vector((point.x, point.y, original.z)), Vector((tangent.x, tangent.y, tangent.z)), length, "mrna_path_closest_xy", None
    return original, None, None, "manifest_location", None


def transform_for_asset(asset: dict, scale: float, path_context: dict[str, SampledPath]) -> dict:
    atoms, points = points_for_asset(asset)
    all_vectors = [Vector(point["pos_A"]) for point in points]
    if not all_vectors:
        raise ValueError(f"No residue points parsed for {asset['pdb_id']}")
    min_v = Vector((min(v.x for v in all_vectors), min(v.y for v in all_vectors), min(v.z for v in all_vectors)))
    max_v = Vector((max(v.x for v in all_vectors), max(v.y for v in all_vectors), max(v.z for v in all_vectors)))
    bbox_A = max_v - min_v

    nucleic_vectors = [Vector(point["pos_A"]) for point in points if point["kind"] == "nucleic"]
    protein_vectors = [Vector(point["pos_A"]) for point in points if point["kind"] == "protein"]
    surface_vectors = surface_alignment_vectors(asset)
    nucleic_alignment_vectors = surface_vectors.get("nucleic", nucleic_vectors)
    protein_alignment_vectors = surface_vectors.get("protein", protein_vectors)
    anchor_kind = asset.get("anchor_kind", "center")
    anchor_vectors = nucleic_alignment_vectors if anchor_kind == "nucleic" and nucleic_alignment_vectors else all_vectors
    anchor_A = base.bbox_center(anchor_vectors)
    anchor_mode_detail = anchor_kind
    alignment_coordinate_source = "reduced_surface_obj_vertices" if surface_vectors else "mmcif_residue_points"
    if is_wrapped_dna_asset(asset) and protein_alignment_vectors:
        anchor_A = base.bbox_center(protein_alignment_vectors)
        anchor_mode_detail = "wrapped_protein_surface_center_with_nucleic_surface_plane_alignment"
    location, path_axis, path_length, anchor_mode, path_attachment = path_adjusted_location_and_axis(asset, path_context)
    loop_feature = dna_loop_feature(path_context) if is_wrapped_dna_asset(asset) else None

    align_rotation = Matrix.Identity(3)
    axis_before = None
    binding_alignment = None
    rotation_deg = asset.get("rotation_deg", [0.0, 0.0, 0.0])

    if is_wrapped_dna_asset(asset) and loop_feature and loop_feature.get("basis") == "pdb_guided_base_pair_centerline":
        source_center_A = Vector(loop_feature["source_center_A"])
        export_origin_A = surface_export_origin_A(nucleic_vectors, nucleic_alignment_vectors)
        anchor_A = source_center_A - export_origin_A
        anchor_mode_detail = "pdb_guided_nucleosome_dna_centerline"
        rotation = rotation_from_basis_report(loop_feature)
        source_z = Vector(loop_feature["source_basis"][2])
        target_z = Vector(loop_feature["target_basis"][2])
        axis_before = [source_z.x, source_z.y, source_z.z]
        binding_alignment = {
            "mode": "same_transform_as_procedural_nucleosome_dna_guide",
            "anchor_mode": anchor_mode,
            "guide_source_pdb_id": loop_feature.get("source_pdb_id"),
            "guide_source_bp_count": loop_feature.get("source_bp_count"),
            "guide_length_mm": loop_feature.get("approximate_centerline_length_mm"),
            "guide_bp_equivalent": loop_feature.get("approximate_bp_equivalent"),
            "source_center_A": loop_feature.get("source_center_A"),
            "surface_export_origin_A": [export_origin_A.x, export_origin_A.y, export_origin_A.z],
            "source_center_in_surface_obj_A": [anchor_A.x, anchor_A.y, anchor_A.z],
            "target_center_mm": loop_feature.get("center_mm"),
            "source_plane_normal_before": [source_z.x, source_z.y, source_z.z],
            "target_plane_normal": [target_z.x, target_z.y, target_z.z],
            "guide_roll_deg": loop_feature.get("roll_deg"),
            "manual_asset_rotation_deg_ignored": rotation_deg,
        }
    elif is_wrapped_dna_asset(asset) and len(nucleic_alignment_vectors) >= 3:
        axes = covariance_axes(nucleic_alignment_vectors)
        source_plane_normal = axes[-1][1]
        target_axis = Vector((0.0, 0.0, 1.0))
        if source_plane_normal.dot(target_axis) < 0:
            source_plane_normal = -source_plane_normal
        normal_alignment = source_plane_normal.rotation_difference(target_axis).to_matrix()
        source_in_plane = normal_alignment @ axes[0][1]
        roll_alignment, roll_deg = axis_roll_rotation(target_axis, source_in_plane, Vector((1.0, 0.0, 0.0)))
        extra_roll = Matrix.Rotation(math.radians(rotation_deg[2]), 3, target_axis)
        rotation = extra_roll @ roll_alignment @ normal_alignment
        axis_before = [source_plane_normal.x, source_plane_normal.y, source_plane_normal.z]
        binding_alignment = {
            "mode": "wrapped_loop_plane_to_scene",
            "anchor_mode": anchor_mode,
            "path_length_mm": path_length,
            "source_plane_normal_before": axis_before,
            "target_plane_normal": [0.0, 0.0, 1.0],
            "in_plane_roll_deg": roll_deg,
            "manual_roll_about_target_axis_deg": rotation_deg[2],
        }
    elif (
        path_attachment
        and path_attachment.get("path") == "mrna"
        and str(path_attachment.get("binding_mode", "")).startswith("co_crystal_rna")
        and len(nucleic_alignment_vectors) >= 2
    ):
        axes = covariance_axes(nucleic_alignment_vectors)
        source_axis = axes[0][1]
        target_axis = path_axis.normalized() if path_axis is not None else Vector((1.0, 0.0, 0.0))
        if source_axis.dot(target_axis) < 0:
            source_axis = -source_axis
        align_rotation = source_axis.rotation_difference(target_axis).to_matrix()
        protein_contact_A, nucleic_contact_A, source_contact_gap_A = closest_surface_pair_A(protein_alignment_vectors, nucleic_alignment_vectors)
        if protein_contact_A is not None and nucleic_contact_A is not None:
            anchor_A = protein_contact_A
            anchor_mode_detail = "source_protein_rna_interface_surface_point"
            source_binding = base.bbox_center(protein_alignment_vectors) - nucleic_contact_A
        elif protein_alignment_vectors:
            source_binding = base.mean_vector(protein_alignment_vectors) - anchor_A
        else:
            source_binding = axes[1][1]
        source_binding = project_perpendicular(source_binding, source_axis)
        if source_binding.length < 1e-7:
            source_binding = axes[1][1]
        target_side = target_binding_side(asset, target_axis)
        roll_alignment, roll_deg = axis_roll_rotation(target_axis, align_rotation @ source_binding, target_side)
        extra_roll = Matrix.Rotation(math.radians(rotation_deg[2]), 3, target_axis)
        rotation = extra_roll @ roll_alignment @ align_rotation
        target_radius = path_contact_radius_mm(path_attachment.get("path"), path_context)
        target_contact = location + target_side * target_radius
        location = target_contact
        path_attachment["target_contact_radius_mm"] = target_radius
        path_attachment["target_surface_contact_mm"] = [target_contact.x, target_contact.y, target_contact.z]
        axis_before = [source_axis.x, source_axis.y, source_axis.z]
        binding_alignment = {
            "mode": "co_crystal_rna_interface_surface_to_procedural_rna_surface",
            "anchor_mode": anchor_mode,
            "path_length_mm": path_length,
            "source_rna_axis_before": axis_before,
            "target_path_axis": [target_axis.x, target_axis.y, target_axis.z],
            "source_protein_contact_point_A": [protein_contact_A.x, protein_contact_A.y, protein_contact_A.z] if protein_contact_A else None,
            "source_nucleic_contact_point_A": [nucleic_contact_A.x, nucleic_contact_A.y, nucleic_contact_A.z] if nucleic_contact_A else None,
            "source_surface_gap_A": source_contact_gap_A,
            "source_binding_vector_before": [source_binding.x, source_binding.y, source_binding.z],
            "target_binding_side": [target_side.x, target_side.y, target_side.z],
            "target_contact_radius_mm": target_radius,
            "target_surface_contact_mm": [target_contact.x, target_contact.y, target_contact.z],
            "binding_roll_deg": roll_deg,
            "manual_roll_about_target_axis_deg": rotation_deg[2],
            "binding_mode": path_attachment.get("binding_mode"),
        }
    elif (
        path_attachment
        and path_attachment.get("path") == "mrna"
        and path_attachment.get("binding_mode") == "protein_surface_rna_contact"
        and len(protein_alignment_vectors) >= 3
    ):
        axes = covariance_axes(protein_alignment_vectors)
        source_axis = axes[0][1]
        target_axis = path_axis.normalized() if path_axis is not None else Vector((1.0, 0.0, 0.0))
        if source_axis.dot(target_axis) < 0:
            source_axis = -source_axis
        align_rotation = source_axis.rotation_difference(target_axis).to_matrix()

        source_center_A = base.bbox_center(protein_alignment_vectors)
        source_side = project_perpendicular(axes[1][1], source_axis)
        if source_side.length < 1e-7:
            source_side = project_perpendicular(axes[2][1], source_axis)
        source_side = source_side.normalized() if source_side.length else Vector((0.0, 1.0, 0.0))
        projected_vertices = [((vertex - source_center_A).dot(source_side), vertex) for vertex in protein_alignment_vectors]
        negative_extent, negative_vertex = min(projected_vertices, key=lambda item: item[0])
        positive_extent, positive_vertex = max(projected_vertices, key=lambda item: item[0])
        source_contact_A = negative_vertex if abs(negative_extent) >= abs(positive_extent) else positive_vertex
        source_binding = source_center_A - source_contact_A
        anchor_A = source_contact_A
        anchor_mode_detail = "protein_surface_contact_point_without_nucleic_guide"

        target_side = target_binding_side(asset, target_axis)
        roll_alignment, roll_deg = axis_roll_rotation(target_axis, align_rotation @ source_binding, target_side)
        extra_roll = Matrix.Rotation(math.radians(rotation_deg[2]), 3, target_axis)
        rotation = extra_roll @ roll_alignment @ align_rotation
        target_radius = path_contact_radius_mm(path_attachment.get("path"), path_context)
        target_contact = location + target_side * target_radius
        location = target_contact
        path_attachment["target_contact_radius_mm"] = target_radius
        path_attachment["target_surface_contact_mm"] = [target_contact.x, target_contact.y, target_contact.z]
        axis_before = [source_axis.x, source_axis.y, source_axis.z]
        binding_alignment = {
            "mode": "protein_surface_contact_to_procedural_rna_path_without_co_crystal_guide",
            "anchor_mode": anchor_mode,
            "path_length_mm": path_length,
            "co_crystal_guide_available": False,
            "source_protein_axis_before": axis_before,
            "target_path_axis": [target_axis.x, target_axis.y, target_axis.z],
            "source_center_A": [source_center_A.x, source_center_A.y, source_center_A.z],
            "source_contact_point_A": [source_contact_A.x, source_contact_A.y, source_contact_A.z],
            "source_binding_vector_before": [source_binding.x, source_binding.y, source_binding.z],
            "target_binding_side": [target_side.x, target_side.y, target_side.z],
            "target_contact_radius_mm": target_radius,
            "target_surface_contact_mm": [target_contact.x, target_contact.y, target_contact.z],
            "binding_roll_deg": roll_deg,
            "manual_roll_about_target_axis_deg": rotation_deg[2],
            "binding_mode": path_attachment.get("binding_mode"),
        }
    elif is_dna_binding_asset(asset) and len(nucleic_alignment_vectors) >= 2:
        axes = covariance_axes(nucleic_alignment_vectors)
        source_axis = axes[0][1]
        target_axis = path_axis.normalized() if path_axis is not None else Vector(asset.get("align_nucleic_axis_to", [1.0, 0.0, 0.0])).normalized()
        if source_axis.dot(target_axis) < 0:
            source_axis = -source_axis
        align_rotation = source_axis.rotation_difference(target_axis).to_matrix()
        protein_contact_A, nucleic_contact_A, source_contact_gap_A = closest_surface_pair_A(protein_alignment_vectors, nucleic_alignment_vectors)
        if protein_contact_A is not None and nucleic_contact_A is not None:
            anchor_A = protein_contact_A
            anchor_mode_detail = "source_protein_dna_interface_surface_point"
            source_binding = base.bbox_center(protein_alignment_vectors) - nucleic_contact_A
        else:
            source_binding = base.mean_vector(protein_alignment_vectors) - anchor_A if protein_alignment_vectors else Vector((0.0, 1.0, 0.0))
        source_binding = project_perpendicular(source_binding, source_axis)
        if source_binding.length < 1e-7:
            source_binding = axes[1][1]
        target_side = target_binding_side(asset, target_axis)
        roll_alignment, roll_deg = axis_roll_rotation(target_axis, align_rotation @ source_binding, target_side)
        extra_roll = Matrix.Rotation(math.radians(rotation_deg[2]), 3, target_axis)
        rotation = extra_roll @ roll_alignment @ align_rotation
        target_radius = path_contact_radius_mm(path_attachment.get("path") if path_attachment else "dna", path_context)
        target_contact = location + target_side * target_radius
        location = target_contact
        if path_attachment:
            path_attachment["target_contact_radius_mm"] = target_radius
            path_attachment["target_surface_contact_mm"] = [target_contact.x, target_contact.y, target_contact.z]
        axis_before = [source_axis.x, source_axis.y, source_axis.z]
        binding_alignment = {
            "mode": "co_crystal_dna_interface_surface_to_procedural_dna_surface",
            "anchor_mode": anchor_mode,
            "path_length_mm": path_length,
            "source_dna_axis_before": axis_before,
            "target_path_axis": [target_axis.x, target_axis.y, target_axis.z],
            "source_protein_contact_point_A": [protein_contact_A.x, protein_contact_A.y, protein_contact_A.z] if protein_contact_A else None,
            "source_nucleic_contact_point_A": [nucleic_contact_A.x, nucleic_contact_A.y, nucleic_contact_A.z] if nucleic_contact_A else None,
            "source_surface_gap_A": source_contact_gap_A,
            "source_binding_vector_before": [source_binding.x, source_binding.y, source_binding.z],
            "target_binding_side": [target_side.x, target_side.y, target_side.z],
            "target_contact_radius_mm": target_radius,
            "target_surface_contact_mm": [target_contact.x, target_contact.y, target_contact.z],
            "binding_roll_deg": roll_deg,
            "manual_roll_about_target_axis_deg": rotation_deg[2],
        }
    else:
        target_axis_data = path_axis if path_axis is not None else asset.get("align_nucleic_axis_to")
        if target_axis_data is not None and len(nucleic_vectors) >= 2:
            source_axis = base.principal_axis(nucleic_vectors)
            target_axis = Vector(target_axis_data).normalized()
            if source_axis.dot(target_axis) < 0:
                source_axis = -source_axis
            axis_before = [source_axis.x, source_axis.y, source_axis.z]
            align_rotation = source_axis.rotation_difference(target_axis).to_matrix()
        manual_rotation = Euler(tuple(math.radians(value) for value in rotation_deg), "XYZ").to_matrix()
        rotation = manual_rotation @ align_rotation

    def transform(vertex_A: Vector) -> Vector:
        return location + rotation @ ((vertex_A - anchor_A) * scale)

    return {
        "transform": transform,
        "atom_count": len(atoms),
        "residue_points": len(points),
        "bbox_A": [bbox_A.x, bbox_A.y, bbox_A.z],
        "bbox_mm": [bbox_A.x * scale, bbox_A.y * scale, bbox_A.z * scale],
        "anchor_A": [anchor_A.x, anchor_A.y, anchor_A.z],
        "anchor_kind": anchor_kind,
        "anchor_mode_detail": anchor_mode_detail,
        "alignment_coordinate_source": alignment_coordinate_source,
        "alignment_surface_vertex_counts": {component: len(vertices) for component, vertices in surface_vectors.items()},
        "location_mm": [location.x, location.y, location.z],
        "aligned_nucleic_axis_before": axis_before,
        "path_axis_after": [path_axis.x, path_axis.y, path_axis.z] if path_axis is not None else None,
        "path_anchor_mode": anchor_mode,
        "path_anchor_length_mm": path_length,
        "path_attachment": path_attachment,
        "binding_alignment": binding_alignment,
        "points": points,
        "rotation": rotation,
        "anchor_vector": anchor_A,
    }


def parse_obj_transformed(path: Path, transform) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]]]:
    vertices = []
    faces = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            if raw_line.startswith("v "):
                _, x, y, z = raw_line.split()[:4]
                world = transform(Vector((float(x), float(y), float(z))))
                vertices.append((world.x, world.y, world.z))
            elif raw_line.startswith("f "):
                face = []
                for token in raw_line.split()[1:]:
                    index = token.split("/")[0]
                    if index:
                        face.append(int(index) - 1)
                if len(face) >= 3:
                    faces.append(tuple(face))
    return vertices, faces


def import_surface_obj(
    path: Path,
    name: str,
    material: bpy.types.Material,
    collection: bpy.types.Collection,
    transform,
    style: str = "pymol_surface_reduced",
    coordinate_units: str = "angstrom",
) -> tuple[bpy.types.Object, dict]:
    vertices, faces = parse_obj_transformed(path, transform)
    if not vertices or not faces:
        raise ValueError(f"Surface OBJ has no usable geometry: {path}")
    mesh = bpy.data.meshes.new(f"{name}_mesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    for polygon in mesh.polygons:
        polygon.use_smooth = True
    mesh.materials.append(material)
    obj = bpy.data.objects.new(name, mesh)
    base.link_to_collection(obj, collection)
    obj["style"] = style
    obj["source_obj"] = str(path.relative_to(ROOT)).replace("\\", "/")
    obj["coordinate_units"] = coordinate_units
    if coordinate_units == "angstrom":
        obj["angstrom_to_mm"] = 0.04
    else:
        obj["coordinate_to_mm"] = 1.0
    report = mesh_bounds(vertices)
    report.update({"object": obj.name, "source_obj": str(path.relative_to(ROOT)).replace("\\", "/"), "vertices": len(vertices), "faces": len(faces)})
    return obj, report


def attachment_empty_name(asset: dict, path_name: str | None) -> str:
    safe = "".join(char if char.isalnum() else "_" for char in asset["name"]).strip("_")
    prefix = (path_name or "free").upper()
    return f"ATTACH_{prefix}_{safe}"


def create_attachment_empty(asset: dict, transform_data: dict, collections: dict) -> tuple[bpy.types.Object | None, dict | None]:
    attachment = transform_data.get("path_attachment")
    if not attachment:
        return None, None
    empty = bpy.data.objects.new(attachment_empty_name(asset, attachment.get("path")), None)
    empty.empty_display_type = "SPHERE"
    empty.empty_display_size = 2.0
    location = Vector(transform_data["location_mm"])
    empty.location = location
    empty.rotation_euler = transform_data["rotation"].to_euler()
    empty["attachment_path"] = attachment.get("path")
    empty["attachment_distance_mm"] = float(attachment["distance_mm"])
    empty["binding_mode"] = attachment.get("binding_mode")
    base.link_to_collection(empty, collections[asset["collection"]])
    report = {
        "empty": empty.name,
        "path": attachment.get("path"),
        "distance_mm": attachment["distance_mm"],
        "path_length_mm": attachment["path_length_mm"],
        "fraction": attachment["fraction"],
        "binding_mode": attachment.get("binding_mode"),
        "offset_mm": attachment.get("offset_mm"),
        "offset_local_mm": attachment.get("offset_local_mm"),
        "roll_deg": attachment.get("roll_deg"),
        "path_point_mm": attachment.get("path_point_mm"),
        "anchor_point_mm": attachment.get("anchor_point_mm"),
        "target_contact_radius_mm": attachment.get("target_contact_radius_mm"),
        "target_surface_contact_mm": attachment.get("target_surface_contact_mm"),
        "location_mm": [location.x, location.y, location.z],
    }
    return empty, report


def parent_to_attachment_empty(obj: bpy.types.Object, empty: bpy.types.Object) -> None:
    matrix_world = obj.matrix_world.copy()
    obj.parent = empty
    obj.matrix_world = matrix_world
    obj.matrix_parent_inverse = empty.matrix_world.inverted()


def surface_directory_for(pdb_id: str) -> Path:
    reduced = REDUCED_SURFACE_DIR / pdb_id
    return reduced if reduced.exists() else RAW_SURFACE_DIR / pdb_id


def surface_components_for(pdb_id: str) -> list[tuple[str, Path]]:
    directory = surface_directory_for(pdb_id)
    prefix = f"{pdb_id}_surface_"
    components = []
    for path in sorted(directory.glob(f"{prefix}*.obj")):
        if path.stat().st_size == 0:
            continue
        components.append((path.stem[len(prefix):], path))
    return components


def material_for_component(asset: dict, component: str, materials: dict[str, bpy.types.Material]) -> bpy.types.Material:
    nucleic_component = component == "nucleic" or component in {"rRNA", "mRNA"} or component.startswith("tRNA")
    return base.component_material(asset, "nucleic" if nucleic_component else "protein", materials)


def hide_component_used_for_dna_alignment(asset: dict, component: str) -> bool:
    if asset.get("show_alignment_nucleic", False):
        return False
    if asset.get("hide_alignment_nucleic", False) and component in DNA_ALIGNMENT_GUIDE_COMPONENTS:
        return True
    return is_dna_binding_asset(asset) and component in DNA_ALIGNMENT_GUIDE_COMPONENTS


def component_filter(point: dict, component: str) -> bool:
    if component == "protein" or component == "ribosomal_protein":
        return point["kind"] == "protein"
    if component == "nucleic":
        return point["kind"] == "nucleic"
    if component in {"rRNA", "mRNA"} or component.startswith("tRNA_"):
        return point["kind"] == "nucleic"
    return False


def build_component_bead_fallback(
    asset: dict,
    component: str,
    transform_data: dict,
    collections: dict,
    materials: dict,
    scale: float,
) -> dict | None:
    points = [point for point in transform_data["points"] if component_filter(point, component)]
    if not points:
        return None
    sampled = base.sample_points(points, min(int(asset.get("max_beads", 1500)), 3000))
    radius = float(asset["bead_radius_A"]) * scale
    centers = [transform_data["transform"](Vector(point["pos_A"])) for point in sampled]
    obj = base.create_bead_mesh(
        f"{asset['name']} ({asset['pdb_id']}) fallback {component}",
        centers,
        radius,
        material_for_component(asset, component, materials),
        collections[asset["collection"]],
    )
    if obj:
        obj["style"] = "component_residue_beads_fallback"
        obj["pdb_id"] = asset["pdb_id"]
        obj["component"] = component
        obj["angstrom_to_mm"] = scale
    return {"component": component, "style": "component_residue_beads_fallback", "sampled_beads": len(sampled)}


def add_asset_label(asset: dict, location: list[float], collections: dict, materials: dict) -> str:
    label_offset = asset.get("label_offset_mm", [0.0, 8.0, 0.4])
    label_size = float(asset.get("label_size_mm", 1.75))
    label_text = asset.get("label_text") or f"{asset['name']} ({asset['pdb_id'].upper()})"
    obj = base.create_text(
        f"label_{asset['name']}",
        label_text,
        (location[0] + label_offset[0], location[1] + label_offset[1], location[2] + label_offset[2]),
        label_size,
        materials["label_grey"],
        collections["Labels"],
    )
    obj["label_text"] = label_text
    obj["pdb_id"] = asset["pdb_id"].upper()
    return label_text


def import_asset(asset: dict, collections: dict, materials: dict, scale: float, path_context: dict[str, SampledPath]) -> dict:
    pdb_id = asset["pdb_id"].upper()
    transform_data = transform_for_asset(asset, scale, path_context)
    attachment_empty, attachment_report = create_attachment_empty(asset, transform_data, collections)
    components = surface_components_for(pdb_id)
    expected = tuple(component for component, _ in components)
    imported = []
    component_reports = []
    component_by_name = {component: path for component, path in components}
    for component in expected:
        path = component_by_name.get(component)
        if path is None:
            raise FileNotFoundError(f"Missing PyMOL surface for {pdb_id} component {component}")
        print(f"Importing reduced surface {pdb_id} {component}: {path}")
        obj, report = import_surface_obj(
            path,
            f"{asset['name']} ({pdb_id}) surface {component}",
            material_for_component(asset, component, materials),
            collections[asset["collection"]],
            transform_data["transform"],
        )
        obj["pdb_id"] = pdb_id
        obj["component"] = component
        obj["source"] = asset["source"]
        alignment_guide_only = hide_component_used_for_dna_alignment(asset, component)
        if alignment_guide_only:
            obj.hide_viewport = True
            obj.hide_render = True
            obj["alignment_guide_only"] = True
        if attachment_empty is not None:
            # Surface OBJ vertices are transformed into scene/world millimeter
            # coordinates before the mesh object is created. Parenting those
            # baked meshes to the attachment empty applies the empty transform
            # again after save/reload in Blender, visibly moving the mesh away
            # from the intended nucleic-acid contact. Keep the empty as a
            # marker/report object only.
            obj["attachment_empty"] = attachment_empty.name
        imported.append(obj)
        report["component"] = component
        report["alignment_guide_only"] = alignment_guide_only
        component_reports.append(report)
    if not imported:
        raise FileNotFoundError(f"No PyMOL surface components available for {pdb_id}")
    label_text = add_asset_label(asset, transform_data["location_mm"], collections, materials)
    return {
        "name": asset["name"],
        "pdb_id": pdb_id,
        "label_text": label_text,
        "style": "pymol_surface_reduced",
        "atom_count": transform_data["atom_count"],
        "residue_points": transform_data["residue_points"],
        "bbox_A": transform_data["bbox_A"],
        "bbox_mm": transform_data["bbox_mm"],
        "location_mm": transform_data["location_mm"],
        "anchor_kind": transform_data["anchor_kind"],
        "anchor_mode_detail": transform_data["anchor_mode_detail"],
        "alignment_coordinate_source": transform_data["alignment_coordinate_source"],
        "alignment_surface_vertex_counts": transform_data["alignment_surface_vertex_counts"],
        "anchor_A": transform_data["anchor_A"],
        "aligned_nucleic_axis_before": transform_data["aligned_nucleic_axis_before"],
        "path_axis_after": transform_data["path_axis_after"],
        "path_anchor_mode": transform_data["path_anchor_mode"],
        "path_anchor_length_mm": transform_data["path_anchor_length_mm"],
        "attachment_empty": attachment_report,
        "binding_alignment": transform_data["binding_alignment"],
        "components": component_reports,
    }


def build_assets(manifest: dict, collections: dict, materials: dict, path_context: dict[str, SampledPath]) -> list[dict]:
    scale = manifest["units"]["angstrom_to_mm"]
    return [import_asset(asset, collections, materials, scale, path_context) for asset in manifest["pdb_assets"]]


def add_scale_bars(manifest: dict, collections: dict, materials: dict) -> dict:
    return base.build_scale_bars(manifest, collections, materials)


def add_camera_views(collections: dict, materials: dict) -> None:
    base.add_lighting_and_camera(collections, materials)
    primary_camera = bpy.context.scene.camera
    if primary_camera:
        primary_camera.location = (0.0, 2.0, 285.0)
        primary_camera.data.ortho_scale = 245.0
    for name, location, ortho_scale in [
        ("Camera_DNA_transcription_detail", (-22.0, -44.0, 120.0), 82.0),
        ("Camera_Cas9_binding_detail", (-97.5, -47.0, 90.0), 18.0),
        ("Camera_Nucleosome_binding_detail", (72.0, -45.0, 90.0), 18.0),
        ("Camera_translation_detail", (12.0, 42.0, 120.0), 86.0),
        ("Camera_ribosome_surface_detail", (-13.5, 40.0, 90.0), 42.0),
    ]:
        camera_data = bpy.data.cameras.new(name)
        camera = bpy.data.objects.new(name, camera_data)
        camera.location = location
        camera.rotation_euler = (0.0, 0.0, 0.0)
        camera_data.type = "ORTHO"
        camera_data.ortho_scale = ortho_scale
        bpy.context.scene.collection.objects.link(camera)


def validate_scene(report: dict) -> None:
    expected_scale = report["units"]["angstrom_to_mm"]
    mismatches = []
    fallback_objects = []
    surface_style_objects = 0
    for obj in bpy.data.objects:
        if "angstrom_to_mm" in obj and abs(float(obj["angstrom_to_mm"]) - expected_scale) > 1e-9:
            mismatches.append(obj.name)
        style = str(obj.get("style", ""))
        if "fallback" in style:
            fallback_objects.append(obj.name)
        if style in {"pymol_surface_reduced", NUCLEIC_PROXY_STYLE, DIRECT_NUCLEIC_STYLE}:
            surface_style_objects += 1
    report["scale_validation"] = {
        "expected_angstrom_to_mm": expected_scale,
        "objects_with_scale_property": sum(1 for obj in bpy.data.objects if "angstrom_to_mm" in obj),
        "mismatches": mismatches,
    }
    report["surface_style_validation"] = {
        "surface_style_objects": surface_style_objects,
        "fallback_objects": fallback_objects,
        "allowed_styles": ["pymol_surface_reduced", NUCLEIC_PROXY_STYLE, DIRECT_NUCLEIC_STYLE],
    }
    report["object_count"] = len(bpy.data.objects)
    report["mesh_count"] = len(bpy.data.meshes)
    if mismatches:
        raise RuntimeError(f"Scale validation failed: {mismatches}")
    if fallback_objects:
        raise RuntimeError(f"Fallback geometry is not allowed in V3: {fallback_objects}")


def load_reduction_summary() -> dict:
    if not REDUCED_MANIFEST_PATH.exists():
        return {"available": False}
    data = json.loads(REDUCED_MANIFEST_PATH.read_text(encoding="utf-8"))
    ribosome = [entry for entry in data.get("entries", []) if entry.get("pdb_id") in RIBOSOME_PDBS]
    return {
        "available": True,
        "manifest": str(REDUCED_MANIFEST_PATH),
        "ribosome_entries": ribosome,
        "target_failures": [entry for entry in data.get("entries", []) if not entry.get("target_met", True)],
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = base.load_manifest()
    clean_scene()
    configure_scene()
    collections = base.make_collections(manifest["collections"])
    materials = base.make_materials()
    soften_materials(materials)

    dna = build_procedural_dna(manifest, collections, materials)
    mrna = build_procedural_mrna(manifest, collections, materials)
    compact_mrna = build_compact_mrna(manifest, collections, materials)
    path_context = {"dna": dna["path"], "dna_features": dna["features"], "mrna": mrna["path"]}

    report = {
        "title": f"{manifest['title']} - V3 direct Blender DNA/RNA",
        "units": manifest["units"],
        "source_manifest": str(base.MANIFEST_PATH),
        "surface_asset_dir": str(REDUCED_SURFACE_DIR),
        "nucleic_acid_pipeline": "direct_blender_meshes_for_procedural_dna_mrna_and_compact_mrna",
        "reduction": load_reduction_summary(),
        "outputs": {
            "blend": str(BLEND_PATH),
            "preview": str(PREVIEW_PATH),
            "report": str(REPORT_PATH),
            "detail_previews": {name: str(path) for name, path in DETAIL_PREVIEWS.items()},
        },
    }
    report["dna"] = dna["report"]
    report["mrna"] = mrna["report"]
    report["compact_mrna"] = compact_mrna["report"]
    report["scale_bars"] = add_scale_bars(manifest, collections, materials)
    report["pdb_assets"] = build_assets(manifest, collections, materials, path_context)
    add_camera_views(collections, materials)
    validate_scene(report)

    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_PATH))
    bpy.context.scene.render.filepath = str(PREVIEW_PATH)
    bpy.ops.render.render(write_still=True)
    primary_camera = bpy.context.scene.camera
    for key, camera_name in [
        ("dna_transcription", "Camera_DNA_transcription_detail"),
        ("cas9_binding", "Camera_Cas9_binding_detail"),
        ("nucleosome_binding", "Camera_Nucleosome_binding_detail"),
        ("translation", "Camera_translation_detail"),
        ("ribosome", "Camera_ribosome_surface_detail"),
    ]:
        camera = bpy.data.objects.get(camera_name)
        if camera:
            bpy.context.scene.camera = camera
            bpy.context.scene.render.filepath = str(DETAIL_PREVIEWS[key])
            bpy.ops.render.render(write_still=True)
    if primary_camera:
        bpy.context.scene.camera = primary_camera
    print(f"Wrote {BLEND_PATH}")
    print(f"Wrote {PREVIEW_PATH}")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
