#!/usr/bin/env python3
"""Direct Blender mesh builders for scale-correct procedural DNA/RNA."""

from __future__ import annotations

import math
from typing import Iterable

import bpy
from mathutils import Vector
from mathutils.kdtree import KDTree

import build_gene_expression_scene as base
import procedural_nucleic_geometry as geom


def as_vector(point) -> Vector:
    if isinstance(point, Vector):
        return point
    if hasattr(point, "x") and hasattr(point, "y") and hasattr(point, "z"):
        return Vector((point.x, point.y, point.z))
    return Vector(point)


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


def tangent_for(points: list[Vector], index: int) -> Vector:
    if len(points) < 2:
        return Vector((1.0, 0.0, 0.0))
    if index == 0:
        tangent = points[1] - points[0]
    elif index == len(points) - 1:
        tangent = points[-1] - points[-2]
    else:
        tangent = points[index + 1] - points[index - 1]
    return tangent.normalized() if tangent.length else Vector((1.0, 0.0, 0.0))


def bounds(vertices: list[tuple[float, float, float]]) -> dict:
    min_v = Vector((min(v[0] for v in vertices), min(v[1] for v in vertices), min(v[2] for v in vertices)))
    max_v = Vector((max(v[0] for v in vertices), max(v[1] for v in vertices), max(v[2] for v in vertices)))
    return {
        "min_mm": [min_v.x, min_v.y, min_v.z],
        "max_mm": [max_v.x, max_v.y, max_v.z],
        "bbox_mm": [max_v.x - min_v.x, max_v.y - min_v.y, max_v.z - min_v.z],
    }


def mesh_topology_stats(obj: bpy.types.Object) -> dict:
    edge_counts: dict[tuple[int, int], int] = {}
    for polygon in obj.data.polygons:
        vertices = list(polygon.vertices)
        for i, a in enumerate(vertices):
            b = vertices[(i + 1) % len(vertices)]
            edge = (a, b) if a < b else (b, a)
            edge_counts[edge] = edge_counts.get(edge, 0) + 1
    return {
        "vertices": len(obj.data.vertices),
        "faces": len(obj.data.polygons),
        "boundary_edges": sum(1 for count in edge_counts.values() if count == 1),
        "non_manifold_edges": sum(1 for count in edge_counts.values() if count > 2),
    }


def polish_object(obj: bpy.types.Object, voxel_size: float, smooth_factor: float, smooth_iterations: int) -> dict:
    before = mesh_topology_stats(obj)
    bpy.ops.object.select_all(action="DESELECT")
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj
    remesh = obj.modifiers.new("direct_voxel_surface_union", "REMESH")
    remesh.mode = "VOXEL"
    remesh.voxel_size = voxel_size
    remesh.adaptivity = 0.0
    remesh.use_smooth_shade = True
    bpy.ops.object.modifier_apply(modifier=remesh.name)
    smooth = obj.modifiers.new("direct_surface_relax", "SMOOTH")
    smooth.factor = smooth_factor
    smooth.iterations = smooth_iterations
    bpy.ops.object.modifier_apply(modifier=smooth.name)
    weighted = obj.modifiers.new("direct_weighted_normals", "WEIGHTED_NORMAL")
    bpy.ops.object.modifier_apply(modifier=weighted.name)
    bpy.ops.object.shade_smooth()
    return {
        "voxel_size_mm": voxel_size,
        "smooth_factor": smooth_factor,
        "smooth_iterations": smooth_iterations,
        "before": before,
        "after": mesh_topology_stats(obj),
    }


def create_mesh_object(
    name: str,
    vertices: list[tuple[float, float, float]],
    faces: list[tuple[int, ...]],
    material: bpy.types.Material,
    collection: bpy.types.Collection,
    style: str,
    polish: bool = True,
    voxel_size: float = 0.075,
    smooth_factor: float = 0.08,
    smooth_iterations: int = 1,
) -> tuple[bpy.types.Object, dict]:
    mesh = bpy.data.meshes.new(f"{name}_mesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    for polygon in mesh.polygons:
        polygon.use_smooth = True
    mesh.materials.append(material)
    obj = bpy.data.objects.new(name, mesh)
    base.link_to_collection(obj, collection)
    obj["style"] = style
    obj["coordinate_units"] = "millimeter"
    obj["coordinate_to_mm"] = 1.0
    report = bounds(vertices)
    report.update({"object": obj.name, "vertices": len(vertices), "faces": len(faces), "style": style})
    if polish:
        report["polish"] = polish_object(obj, voxel_size, smooth_factor, smooth_iterations)
        report.update(mesh_topology_stats(obj))
    return obj, report


def dna_annotation_segments(settings: dict, bp_spacing_mm: float) -> list[dict]:
    annotation = settings.get("gene_annotation") or {}
    segments = annotation.get("segments") or []
    resolved = []
    cursor = 0
    for segment in segments:
        length_bp = int(segment.get("length_bp", 0))
        start_bp = int(segment.get("start_bp", cursor))
        end_bp = int(segment.get("end_bp", start_bp + length_bp))
        cursor = end_bp
        if end_bp <= start_bp:
            continue
        resolved.append(
            {
                "name": segment.get("name", segment.get("kind", "segment")),
                "kind": segment.get("kind", "segment"),
                "start_bp": start_bp,
                "end_bp": end_bp,
                "start_mm": start_bp * bp_spacing_mm,
                "end_mm": end_bp * bp_spacing_mm,
                "material": segment.get("material", "dna_orange"),
            }
        )
    return resolved


def apply_dna_segment_materials(
    obj: bpy.types.Object,
    path: geom.SampledPath,
    settings: dict,
    materials: dict,
    bp_spacing_mm: float,
) -> dict | None:
    segments = dna_annotation_segments(settings, bp_spacing_mm)
    if not segments:
        return None
    mesh = obj.data
    mesh.materials.clear()
    material_keys = []
    for segment in segments:
        key = segment["material"] if segment["material"] in materials else "dna_orange"
        if key not in material_keys:
            material_keys.append(key)
            mesh.materials.append(materials[key])
    material_index = {key: index for index, key in enumerate(material_keys)}

    sample_spacing = max(bp_spacing_mm * 2.0, 0.25)
    sample_count = max(2, int(math.ceil(path.length / sample_spacing)) + 1)
    tree = KDTree(sample_count)
    cumulative_by_index = []
    for index in range(sample_count):
        distance = path.length * index / (sample_count - 1)
        point = as_vector(path.point_at_length(distance))
        tree.insert(point, index)
        cumulative_by_index.append(distance)
    tree.balance()

    face_counts = {segment["name"]: 0 for segment in segments}
    face_counts["unassigned"] = 0
    for polygon in mesh.polygons:
        _, sample_index, _ = tree.find(polygon.center)
        bp_position = cumulative_by_index[sample_index] / bp_spacing_mm
        selected = None
        for segment in segments:
            if segment["start_bp"] <= bp_position < segment["end_bp"]:
                selected = segment
                break
        if selected is None and segments and bp_position >= segments[-1]["end_bp"]:
            selected = segments[-1]
        if selected is None:
            face_counts["unassigned"] += 1
            continue
        key = selected["material"] if selected["material"] in materials else "dna_orange"
        polygon.material_index = material_index[key]
        face_counts[selected["name"]] += 1
    mesh.update()
    return {
        "basis": "nearest_centerline_distance_after_direct_mesh_generation",
        "material_keys": material_keys,
        "sample_spacing_mm": sample_spacing,
        "sample_count": sample_count,
        "segments": segments,
        "face_counts": face_counts,
    }


def append_mesh(
    vertices: list[tuple[float, float, float]],
    faces: list[tuple[int, ...]],
    new_vertices: list[tuple[float, float, float]],
    new_faces: list[tuple[int, ...]],
) -> None:
    offset = len(vertices)
    vertices.extend(new_vertices)
    faces.extend(tuple(offset + index for index in face) for face in new_faces)


def tube_mesh(
    points: Iterable,
    radius: float,
    sides: int,
    bump_amplitude: float = 0.0,
    cap: bool = True,
) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]]]:
    centers = [as_vector(point) for point in points]
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, ...]] = []
    if len(centers) < 2:
        return vertices, faces

    for i, center in enumerate(centers):
        normal, binormal = path_frame(tangent_for(centers, i))
        for j in range(sides):
            angle = 2.0 * math.pi * j / sides
            bump = 1.0 + bump_amplitude * (
                0.55 * math.sin(i * 0.73 + j * 1.91)
                + 0.45 * math.sin(i * 1.37 - j * 0.82)
            )
            ring_radius = max(radius * 0.5, radius * bump)
            point = center + normal * (math.cos(angle) * ring_radius) + binormal * (math.sin(angle) * ring_radius)
            vertices.append((point.x, point.y, point.z))

    for i in range(len(centers) - 1):
        ring = i * sides
        next_ring = (i + 1) * sides
        for j in range(sides):
            faces.append((ring + j, ring + (j + 1) % sides, next_ring + (j + 1) % sides, next_ring + j))

    if cap:
        start_center_index = len(vertices)
        vertices.append(tuple(centers[0]))
        end_center_index = len(vertices)
        vertices.append(tuple(centers[-1]))
        for j in range(sides):
            faces.append((start_center_index, (j + 1) % sides, j))
            last = (len(centers) - 1) * sides
            faces.append((end_center_index, last + j, last + (j + 1) % sides))
    return vertices, faces


def cylinder_segments_mesh(
    segments: Iterable[tuple[Vector, Vector]],
    radius: float,
    sides: int,
    bump_amplitude: float = 0.0,
) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]]]:
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, ...]] = []
    for segment_index, (start, end) in enumerate(segments):
        tangent = end - start
        if tangent.length < 1e-9:
            continue
        normal, binormal = path_frame(tangent)
        offset = len(vertices)
        for center_index, center in enumerate((start, end)):
            for j in range(sides):
                angle = 2.0 * math.pi * j / sides
                bump = 1.0 + bump_amplitude * math.sin(segment_index * 0.59 + j * 1.41 + center_index)
                ring_radius = max(radius * 0.5, radius * bump)
                point = center + normal * (math.cos(angle) * ring_radius) + binormal * (math.sin(angle) * ring_radius)
                vertices.append((point.x, point.y, point.z))
        for j in range(sides):
            faces.append((offset + j, offset + (j + 1) % sides, offset + sides + (j + 1) % sides, offset + sides + j))
        start_center_index = len(vertices)
        vertices.append((start.x, start.y, start.z))
        end_center_index = len(vertices)
        vertices.append((end.x, end.y, end.z))
        for j in range(sides):
            faces.append((start_center_index, offset + (j + 1) % sides, offset + j))
            faces.append((end_center_index, offset + sides + j, offset + sides + (j + 1) % sides))
    return vertices, faces


def append_ellipsoid_mesh(
    vertices: list[tuple[float, float, float]],
    faces: list[tuple[int, ...]],
    center: Vector,
    axis_u: Vector,
    axis_v: Vector,
    axis_w: Vector,
    radius_u: float,
    radius_v: float,
    radius_w: float,
    rings: int,
    sides: int,
    bump_amplitude: float = 0.0,
    phase: float = 0.0,
) -> None:
    """Append a closed low-poly ellipsoid oriented by a local molecular frame."""
    axis_u = axis_u.normalized() if axis_u.length else Vector((1.0, 0.0, 0.0))
    axis_v = axis_v.normalized() if axis_v.length else Vector((0.0, 1.0, 0.0))
    axis_w = axis_w.normalized() if axis_w.length else Vector((0.0, 0.0, 1.0))
    offset = len(vertices)
    vertices.append(tuple(center + axis_w * radius_w))
    for ring in range(1, rings):
        polar = math.pi * ring / rings
        sin_p = math.sin(polar)
        cos_p = math.cos(polar)
        for side in range(sides):
            azimuth = 2.0 * math.pi * side / sides
            bump = 1.0 + bump_amplitude * math.sin(phase + ring * 1.73 + side * 2.17)
            point = (
                center
                + axis_u * (math.cos(azimuth) * sin_p * radius_u * bump)
                + axis_v * (math.sin(azimuth) * sin_p * radius_v * bump)
                + axis_w * (cos_p * radius_w * bump)
            )
            vertices.append(tuple(point))
    bottom_index = len(vertices)
    vertices.append(tuple(center - axis_w * radius_w))

    if rings <= 1:
        return
    first_ring = offset + 1
    for side in range(sides):
        faces.append((offset, first_ring + side, first_ring + (side + 1) % sides))
    for ring in range(rings - 2):
        current = first_ring + ring * sides
        next_ring = current + sides
        for side in range(sides):
            faces.append((current + side, current + (side + 1) % sides, next_ring + (side + 1) % sides, next_ring + side))
    last_ring = first_ring + (rings - 2) * sides
    for side in range(sides):
        faces.append((bottom_index, last_ring + (side + 1) % sides, last_ring + side))


def sampled_path_points(path: geom.SampledPath, spacing_mm: float) -> list[Vector]:
    count = max(2, int(math.ceil(path.length / spacing_mm)) + 1)
    return [as_vector(path.point_at_length(path.length * i / (count - 1))) for i in range(count)]


def build_dna_meshes(manifest: dict, collections: dict, materials: dict) -> tuple[geom.SampledPath, dict]:
    settings = manifest.get("procedural_nucleic_acids", {}).get("dna", {})
    model = geom.build_dna_model(manifest)
    style = settings.get("style", "procedural_blender_surface_proxy")
    strand_radius = float(settings.get("strand_radius_mm", model["report"]["strand_radius_mm"]))
    strand_sides = int(settings.get("strand_sides", 9))
    base_radius = float(settings.get("base_pair_radius_mm", 0.034))
    base_sides = int(settings.get("base_pair_sides", 6))
    base_stack_radius = float(settings.get("base_stack_radius_mm", 0.16))
    base_stack_sides = int(settings.get("base_stack_sides", 10))
    every_bp = max(1, int(settings.get("base_pair_every_bp", 3)))
    bump = float(settings.get("surface_bump_amplitude", 0.08))
    voxel_size = float(settings.get("direct_voxel_size_mm", 0.075))
    smooth_factor = float(settings.get("direct_smooth_factor", 0.08))
    smooth_iterations = int(settings.get("direct_smooth_iterations", 1))
    unified = bool(settings.get("direct_unified_mesh", False))

    strand_a_vertices, strand_a_faces = tube_mesh(model["strand_a"], strand_radius, strand_sides, bump)
    strand_b_vertices, strand_b_faces = tube_mesh(model["strand_b"], strand_radius, strand_sides, bump)
    bridge_segments = [
        (as_vector(model["strand_a"][i]), as_vector(model["strand_b"][i]))
        for i in range(0, len(model["strand_a"]), every_bp)
    ]
    base_vertices, base_faces = cylinder_segments_mesh(bridge_segments, base_radius, base_sides, bump * 0.5)
    axis_points = sampled_path_points(model["path"], max(0.068, model["report"]["bp_spacing_mm"] * 0.5))
    core_vertices, core_faces = tube_mesh(axis_points, base_stack_radius, base_stack_sides, bump * 0.6)

    if unified:
        unified_vertices: list[tuple[float, float, float]] = []
        unified_faces: list[tuple[int, ...]] = []
        append_mesh(unified_vertices, unified_faces, core_vertices, core_faces)
        append_mesh(unified_vertices, unified_faces, base_vertices, base_faces)
        append_mesh(unified_vertices, unified_faces, strand_a_vertices, strand_a_faces)
        append_mesh(unified_vertices, unified_faces, strand_b_vertices, strand_b_faces)
        obj_unified, report_unified = create_mesh_object(
            "DNA B-form direct fused PyMOL-like surface",
            unified_vertices,
            unified_faces,
            materials["dna_orange"],
            collections["DNA"],
            style,
            voxel_size=voxel_size,
            smooth_factor=smooth_factor,
            smooth_iterations=smooth_iterations,
        )
        obj_unified["component"] = "fused_dna_surface"
        segment_material_report = apply_dna_segment_materials(
            obj_unified,
            model["path"],
            settings,
            materials,
            model["report"]["bp_spacing_mm"],
        )
        component_reports = [
            dict(report_unified, component="fused_dna_surface", geometry="strand_base_stack_voxel_union"),
            {
                "component": "strand_A_source",
                "vertices": len(strand_a_vertices),
                "faces": len(strand_a_faces),
                "geometry": "source_primitives",
            },
            {
                "component": "strand_B_source",
                "vertices": len(strand_b_vertices),
                "faces": len(strand_b_faces),
                "geometry": "source_primitives",
            },
            {
                "component": "base_pairs_source",
                "vertices": len(base_vertices),
                "faces": len(base_faces),
                "represented_every_bp": every_bp,
                "geometry": "source_primitives",
            },
            {
                "component": "base_stack_fill_source",
                "vertices": len(core_vertices),
                "faces": len(core_faces),
                "sample_spacing_mm": max(0.068, model["report"]["bp_spacing_mm"] * 0.5),
                "geometry": "source_primitives",
            },
        ]
    else:
        segment_material_report = None
        obj_a, report_a = create_mesh_object(
            "DNA B-form direct surface strand_A",
            strand_a_vertices,
            strand_a_faces,
            materials["dna_orange"],
            collections["DNA"],
            style,
            voxel_size=voxel_size,
            smooth_factor=smooth_factor,
            smooth_iterations=smooth_iterations,
        )
        obj_a["component"] = "strand_A"

        obj_b, report_b = create_mesh_object(
            "DNA B-form direct surface strand_B",
            strand_b_vertices,
            strand_b_faces,
            materials["dna_dark"],
            collections["DNA"],
            style,
            voxel_size=voxel_size,
            smooth_factor=smooth_factor,
            smooth_iterations=smooth_iterations,
        )
        obj_b["component"] = "strand_B"

        obj_base, report_base = create_mesh_object(
            "DNA B-form direct surface base_pairs",
            base_vertices,
            base_faces,
            materials["rna_gold"],
            collections["DNA"],
            style,
            voxel_size=voxel_size,
            smooth_factor=smooth_factor,
            smooth_iterations=smooth_iterations,
        )
        obj_base["component"] = "base_pairs"

        obj_core, report_core = create_mesh_object(
            "DNA B-form direct surface base_stack_fill",
            core_vertices,
            core_faces,
            materials["rna_gold"],
            collections["DNA"],
            style,
            voxel_size=voxel_size,
            smooth_factor=smooth_factor,
            smooth_iterations=smooth_iterations,
        )
        obj_core["component"] = "base_stack_fill"
        component_reports = [
            dict(report_a, component="strand_A"),
            dict(report_b, component="strand_B"),
            dict(report_base, component="base_pairs", represented_every_bp=every_bp),
            dict(report_core, component="base_stack_fill", sample_spacing_mm=max(0.068, model["report"]["bp_spacing_mm"] * 0.5)),
        ]

    report = dict(model["report"])
    report.update(
        {
            "style": style,
            "generation_pipeline": "direct_blender_mesh",
            "direct_unified_mesh": unified,
            "components": component_reports,
            "base_pair_every_bp": every_bp,
            "base_stack_radius_mm": base_stack_radius,
        }
    )
    if segment_material_report:
        report["gene_annotation"] = settings.get("gene_annotation", {})
        report["segment_materials"] = segment_material_report
    return model["path"], report


def build_mrna_meshes_from_model(
    manifest: dict,
    collections: dict,
    materials: dict,
    model: dict,
    object_prefix: str,
    style_suffix: str,
) -> tuple[geom.SampledPath, dict]:
    settings = manifest.get("procedural_nucleic_acids", {}).get("mrna", {})
    style = settings.get("style", "procedural_blender_surface_proxy")
    tube_radius = float(settings.get("tube_radius_mm", 0.22))
    tube_sides = int(settings.get("tube_sides", 9))
    bump = float(settings.get("surface_bump_amplitude", 0.12))
    lobe_radius = float(settings.get("base_lobe_radius_mm", 0.11))
    lobe_every = max(1, int(settings.get("base_lobe_every_nt", 8)))
    lobe_sides = int(settings.get("base_lobe_sides", 7))
    base_offset = float(settings.get("base_offset_mm", 0.16))
    detail_every = max(1, int(settings.get("nucleotide_detail_every_nt", lobe_every)))
    phosphate_radius = float(settings.get("phosphate_radius_mm", 0.095))
    sugar_radius = float(settings.get("sugar_radius_mm", 0.085))
    base_radii = settings.get("rna_base_ellipsoid_radii_mm", [0.145, 0.075, 0.05])
    base_connector_radius = float(settings.get("base_connector_radius_mm", 0.045))
    detail_rings = int(settings.get("nucleotide_detail_rings", 5))
    detail_sides = int(settings.get("nucleotide_detail_sides", 8))
    voxel_size = float(settings.get("direct_voxel_size_mm", 0.075))
    smooth_factor = float(settings.get("direct_smooth_factor", 0.08))
    smooth_iterations = int(settings.get("direct_smooth_iterations", 1))
    unified = bool(settings.get("direct_unified_mesh", False))

    component_reports = []
    nt_to_mm = manifest["units"]["mrna_nt_to_mm"]
    nt_cursor = 0
    for segment_model in model["segments"]:
        segment = segment_model["segment"]
        material = materials[segment["color"]]
        vertices, faces = tube_mesh(segment_model["points"], tube_radius, tube_sides, bump)
        segment_vertices: list[tuple[float, float, float]] = []
        segment_faces: list[tuple[int, ...]] = []
        if unified:
            append_mesh(segment_vertices, segment_faces, vertices, faces)
        else:
            obj, report = create_mesh_object(
                f"{object_prefix} direct surface {segment['name']}",
                vertices,
                faces,
                material,
                collections["mRNA"],
                style,
                voxel_size=voxel_size,
                smooth_factor=smooth_factor,
                smooth_iterations=smooth_iterations,
            )
            obj["component"] = segment["name"]
            obj["nt"] = segment["nt"]
            component_reports.append(dict(report, component=segment["name"], nt=segment["nt"], geometry="backbone_tube"))

        path = geom.SampledPath(segment_model["points"])
        lobe_centers = []
        for local_nt in range(0, segment["nt"], lobe_every):
            distance = min(path.length, local_nt * nt_to_mm)
            center = as_vector(path.point_at_length(distance))
            tangent = as_vector(path.tangent_at_length(distance))
            normal, binormal = path_frame(tangent)
            phase = (nt_cursor + local_nt) * 2.399963229728653
            radial = (math.cos(phase) * normal + math.sin(phase) * binormal).normalized()
            lobe_centers.append(center + radial * base_offset)
        lobe_vertices: list[tuple[float, float, float]] = []
        lobe_faces: list[tuple[int, ...]] = []
        for lobe_index, center in enumerate(lobe_centers):
            local_vertices, local_faces = tube_mesh(
                [
                    center + Vector((-0.001, 0.0, 0.0)),
                    center + Vector((0.001, 0.0, 0.0)),
                ],
                lobe_radius * (0.90 + 0.10 * math.sin(lobe_index)),
                lobe_sides,
                bump * 0.4,
            )
            offset = len(lobe_vertices)
            lobe_vertices.extend(local_vertices)
            lobe_faces.extend(tuple(offset + index for index in face) for face in local_faces)
        if lobe_vertices:
            if unified:
                append_mesh(segment_vertices, segment_faces, lobe_vertices, lobe_faces)
            else:
                lobe_obj, lobe_report = create_mesh_object(
                    f"{object_prefix} direct surface {segment['name']} base_lobes",
                    lobe_vertices,
                    lobe_faces,
                    material,
                    collections["mRNA"],
                    style,
                    voxel_size=voxel_size,
                    smooth_factor=smooth_factor,
                    smooth_iterations=smooth_iterations,
                )
                lobe_obj["component"] = f"{segment['name']}_base_lobes"
                component_reports.append(
                    dict(lobe_report, component=f"{segment['name']}_base_lobes", represented_every_nt=lobe_every, geometry="base_lobes")
                )

        detail_vertices: list[tuple[float, float, float]] = []
        detail_faces: list[tuple[int, ...]] = []
        connector_segments: list[tuple[Vector, Vector]] = []
        for local_nt in range(0, segment["nt"], detail_every):
            distance = min(path.length, local_nt * nt_to_mm)
            center = as_vector(path.point_at_length(distance))
            tangent = as_vector(path.tangent_at_length(distance))
            normal, binormal = path_frame(tangent)
            phase = nt_cursor + local_nt
            radial_phase = phase * 2.399963229728653
            radial = (math.cos(radial_phase) * normal + math.sin(radial_phase) * binormal).normalized()
            tangent_cross = tangent.cross(radial)
            circumferential = tangent_cross.normalized() if tangent_cross.length else binormal
            phosphate_center = center + radial * 0.035 + circumferential * (0.035 * math.sin(phase * 1.7))
            sugar_center = center + radial * 0.105 + circumferential * (0.030 * math.cos(phase * 1.1))
            base_center = center + radial * base_offset + circumferential * (0.040 * math.sin(phase * 0.83))
            append_ellipsoid_mesh(
                detail_vertices,
                detail_faces,
                phosphate_center,
                radial,
                tangent,
                circumferential,
                phosphate_radius * 1.08,
                phosphate_radius * 0.92,
                phosphate_radius * 0.86,
                detail_rings,
                detail_sides,
                bump * 0.35,
                phase,
            )
            append_ellipsoid_mesh(
                detail_vertices,
                detail_faces,
                sugar_center,
                tangent,
                radial,
                circumferential,
                sugar_radius * 1.15,
                sugar_radius * 0.95,
                sugar_radius * 0.82,
                detail_rings,
                detail_sides,
                bump * 0.30,
                phase + 0.41,
            )
            append_ellipsoid_mesh(
                detail_vertices,
                detail_faces,
                base_center,
                tangent,
                radial,
                circumferential,
                float(base_radii[0]),
                float(base_radii[1]),
                float(base_radii[2]),
                detail_rings,
                detail_sides,
                bump * 0.25,
                phase + 0.83,
            )
            connector_segments.append((sugar_center, base_center))
        connector_vertices, connector_faces = cylinder_segments_mesh(connector_segments, base_connector_radius, max(6, detail_sides - 2), bump * 0.2)
        connector_offset = len(detail_vertices)
        detail_vertices.extend(connector_vertices)
        detail_faces.extend(tuple(connector_offset + index for index in face) for face in connector_faces)
        if detail_vertices:
            if unified:
                append_mesh(segment_vertices, segment_faces, detail_vertices, detail_faces)
            else:
                detail_obj, detail_report = create_mesh_object(
                    f"{object_prefix} direct surface {segment['name']} nucleotide_detail",
                    detail_vertices,
                    detail_faces,
                    material,
                    collections["mRNA"],
                    style,
                    voxel_size=voxel_size,
                    smooth_factor=smooth_factor,
                    smooth_iterations=smooth_iterations,
                )
                detail_obj["component"] = f"{segment['name']}_nucleotide_detail"
                component_reports.append(
                    dict(
                        detail_report,
                        component=f"{segment['name']}_nucleotide_detail",
                        represented_every_nt=detail_every,
                        geometry="phosphate_sugar_base_ellipsoids",
                    )
                )
        if unified and segment_vertices:
            segment_obj, segment_report = create_mesh_object(
                f"{object_prefix} direct fused surface {segment['name']}",
                segment_vertices,
                segment_faces,
                material,
                collections["mRNA"],
                style,
                voxel_size=voxel_size,
                smooth_factor=smooth_factor,
                smooth_iterations=smooth_iterations,
            )
            segment_obj["component"] = f"{segment['name']}_fused_surface"
            segment_obj["nt"] = segment["nt"]
            component_reports.append(
                dict(
                    segment_report,
                    component=f"{segment['name']}_fused_surface",
                    nt=segment["nt"],
                    represented_every_nt=detail_every,
                    geometry="backbone_nucleotide_voxel_union",
                )
            )
        nt_cursor += segment["nt"]

    report = dict(model["report"])
    report.update(
        {
            "style": f"{style}_{style_suffix}",
            "generation_pipeline": "direct_blender_mesh",
            "direct_unified_mesh": unified,
            "tube_radius_mm": tube_radius,
            "base_lobe_radius_mm": lobe_radius,
            "base_lobe_every_nt": lobe_every,
            "nucleotide_detail_every_nt": detail_every,
            "phosphate_radius_mm": phosphate_radius,
            "sugar_radius_mm": sugar_radius,
            "rna_base_ellipsoid_radii_mm": base_radii,
            "components": component_reports,
        }
    )
    return model["path"], report


def build_mrna_meshes(manifest: dict, collections: dict, materials: dict) -> tuple[geom.SampledPath, dict]:
    return build_mrna_meshes_from_model(
        manifest,
        collections,
        materials,
        geom.build_mrna_model(manifest),
        "actin mRNA elongated",
        "elongated",
    )


def build_compact_mrna_meshes(manifest: dict, collections: dict, materials: dict) -> tuple[geom.SampledPath, dict]:
    return build_mrna_meshes_from_model(
        manifest,
        collections,
        materials,
        geom.build_compact_mrna_model(manifest),
        "actin mRNA compact",
        "compact",
    )
