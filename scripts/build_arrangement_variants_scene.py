#!/usr/bin/env python3
"""Build a DNA/RNA-only arrangement variants scene."""

from __future__ import annotations

import copy
import json
import math
import sys
from pathlib import Path

import bpy
from mathutils import Vector


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import blender_nucleic_meshes as direct_meshes  # noqa: E402
import build_gene_expression_scene as base  # noqa: E402
import build_gene_expression_surface_scene as surface_scene  # noqa: E402
import procedural_nucleic_geometry as geom  # noqa: E402


EXPERIMENT_DIR = ROOT / "experiments" / "arrangement_variants"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
BLEND_PATH = OUTPUT_DIR / "gene_expression_arrangement_variants.blend"
PREVIEW_PATH = OUTPUT_DIR / "preview_gene_expression_arrangement_variants.png"
REPORT_PATH = OUTPUT_DIR / "gene_expression_arrangement_variants_report.json"
VARIANT_PREVIEW_TEMPLATE = "preview_gene_expression_arrangement_variants_{key}.png"

TARGET_DNA_BP = 1900
CANONICAL_DNA_MM = 258.3432570566669
ACTB_CANONICAL_GENE_BP = 3454
UPSTREAM_PROMOTER_BP = 500
ACTB_WITH_PROMOTER_BP = ACTB_CANONICAL_GENE_BP + UPSTREAM_PROMOTER_BP
ACTB_WITH_PROMOTER_MM = ACTB_WITH_PROMOTER_BP * 0.136
TARGET_MRNA_MM = 222.24


def variant_collection(name: str, parent: bpy.types.Collection | None = None) -> bpy.types.Collection:
    collection = bpy.data.collections.new(name)
    if parent:
        parent.children.link(collection)
    else:
        bpy.context.scene.collection.children.link(collection)
    return collection


def path_length(points: list[tuple[float, float, float]]) -> float:
    if len(points) < 2:
        return 0.0
    return geom.SampledPath(geom.catmull_rom(points, 32)).length


def scale_path_to_length(
    points: list[tuple[float, float, float]],
    target_length: float = CANONICAL_DNA_MM,
) -> list[tuple[float, float, float]]:
    length = path_length(points)
    if length <= 0.0:
        return points
    anchor = Vector(points[0])
    scale = target_length / length
    return [(anchor + (Vector(point) - anchor) * scale).to_tuple() for point in points]


def sampled_point(points: list[tuple[float, float, float]], fraction: float) -> Vector:
    path = geom.SampledPath(geom.catmull_rom(points, 32))
    return path.point_at_length(path.length * max(0.0, min(1.0, fraction)))


def centered_sampled_point(points: list[tuple[float, float, float]], center: tuple[float, float, float]) -> tuple[Vector, float]:
    path = geom.SampledPath(geom.catmull_rom(points, 32))
    target = Vector((center[0], center[1], 0.0))
    best_index = 0
    best_distance = float("inf")
    for index, point in enumerate(path.points):
        distance = (Vector((point.x, point.y, 0.0)) - target).length
        if distance < best_distance:
            best_index = index
            best_distance = distance
    length = path.cumulative[best_index]
    fraction = length / path.length if path.length else 0.0
    return path.points[best_index], fraction


def organic_arc_controls(
    center: tuple[float, float, float],
    radius: float,
    start_deg: float,
    span_deg: float,
    samples: int,
    radial_amp: float,
    z_amp: float,
    phase: float,
    target_length: float = CANONICAL_DNA_MM,
) -> list[tuple[float, float, float]]:
    center_v = Vector(center)
    controls = []
    for index in range(samples + 1):
        t = index / samples
        angle = math.radians(start_deg + span_deg * t)
        envelope = math.sin(math.pi * t) ** 0.62
        radial_noise = envelope * radial_amp * (
            0.55 * math.sin(math.tau * (1.35 * t + phase))
            + 0.30 * math.sin(math.tau * (3.40 * t + phase * 0.37))
            + 0.15 * math.sin(math.tau * (7.10 * t + 0.19))
        )
        tangent_noise = envelope * radial_amp * 0.36 * math.sin(math.tau * (2.6 * t + phase * 0.71))
        radial = Vector((math.cos(angle), math.sin(angle), 0.0))
        tangent = Vector((-math.sin(angle), math.cos(angle), 0.0))
        z = z_amp * envelope * (
            0.7 * math.sin(math.tau * (2.2 * t + phase))
            + 0.3 * math.sin(math.tau * (5.4 * t + 0.33))
        )
        point = center_v + radial * (radius + radial_noise) + tangent * tangent_noise + Vector((0.0, 0.0, z))
        controls.append(point.to_tuple())
    return scale_path_to_length(controls, target_length)


def spiral_dna_controls(
    center: tuple[float, float, float],
    start_radius: float,
    end_radius: float,
    turns: float,
    start_deg: float,
    samples: int,
    radial_amp: float,
    z_amp: float,
    phase: float,
    target_length: float = CANONICAL_DNA_MM,
) -> list[tuple[float, float, float]]:
    center_v = Vector(center)
    controls = []
    for index in range(samples + 1):
        t = index / samples
        angle = math.radians(start_deg) + math.tau * turns * t
        ease = t * t * (3.0 - 2.0 * t)
        radius = start_radius * (1.0 - ease) + end_radius * ease
        envelope = math.sin(math.pi * t) ** 0.58
        r = radius + envelope * radial_amp * (
            0.62 * math.sin(math.tau * (2.25 * t + phase))
            + 0.38 * math.sin(math.tau * (5.85 * t + 0.21))
        )
        z = z_amp * envelope * (
            0.55 * math.sin(math.tau * (1.7 * t + phase * 0.5))
            + 0.45 * math.sin(math.tau * (4.8 * t + phase))
        )
        controls.append((center_v.x + r * math.cos(angle), center_v.y + r * math.sin(angle), center_v.z + z))
    return scale_path_to_length(controls, target_length)


def nested_coil_controls(
    center: tuple[float, float, float],
    radius: float,
    samples: int,
    phase: float,
    target_length: float = CANONICAL_DNA_MM,
) -> list[tuple[float, float, float]]:
    center_v = Vector(center)
    controls = []
    for index in range(samples + 1):
        t = index / samples
        angle = math.radians(210.0) - math.tau * 1.18 * t
        radius_drop = 0.36 * t + 0.08 * math.sin(math.tau * (2.0 * t + phase))
        local_radius = radius * (1.0 - radius_drop)
        local_radius += 4.6 * math.sin(math.tau * (3.3 * t + phase)) * (math.sin(math.pi * t) ** 0.55)
        z = 0.75 * math.sin(math.tau * (4.1 * t + phase)) * (math.sin(math.pi * t) ** 0.65)
        controls.append((center_v.x + local_radius * math.cos(angle), center_v.y + local_radius * math.sin(angle), center_v.z + z))
    return scale_path_to_length(controls, target_length)


def serpentine_dna_controls(
    center: tuple[float, float, float],
    width: float,
    height: float,
    lanes: int,
    samples_per_lane: int,
    irregularity: float,
    z_amp: float,
    phase: float,
    length_variation: float = 0.0,
    target_length: float = CANONICAL_DNA_MM,
) -> list[tuple[float, float, float]]:
    center_v = Vector(center)
    controls = []
    lanes = max(2, lanes)
    lane_pitch = height / (lanes - 1)
    turn_samples = max(10, samples_per_lane)
    lane_widths = []
    for lane in range(lanes):
        lane_phase = math.sin(math.tau * (0.37 * lane + phase)) + 0.45 * math.sin(math.tau * (0.19 * lane + phase * 1.7))
        width_scale = 1.0 + length_variation * lane_phase / 1.45
        lane_widths.append(width * max(0.62, min(1.22, width_scale)))

    def lane_y(lane_index: int, local_t: float) -> float:
        lane_t_inner = lane_index / (lanes - 1)
        y_base_inner = -height * 0.5 + height * lane_t_inner
        global_t_inner = (lane_index + local_t) / lanes
        return y_base_inner + irregularity * (
            0.42 * math.sin(math.tau * (1.3 * global_t_inner + phase))
            + 0.18 * math.sin(math.tau * (4.4 * global_t_inner + 0.18))
        )

    for lane in range(lanes):
        lane_t = lane / (lanes - 1)
        y_base = -height * 0.5 + height * lane_t
        lane_width = lane_widths[lane]
        x_start = -lane_width * 0.5 if lane % 2 == 0 else lane_width * 0.5
        x_end = lane_width * 0.5 if lane % 2 == 0 else -lane_width * 0.5
        for step in range(samples_per_lane):
            local_t = step / max(1, samples_per_lane - 1)
            x = x_start + (x_end - x_start) * local_t
            global_t = (lane + local_t) / lanes
            envelope = math.sin(math.pi * min(1.0, max(0.0, global_t))) ** 0.35
            y = lane_y(lane, local_t)
            z = z_amp * envelope * math.sin(math.tau * (5.2 * global_t + phase))
            controls.append((center_v.x + x, center_v.y + y, center_v.z + z))
        if lane < lanes - 1:
            next_width = lane_widths[lane + 1]
            side = 1.0 if lane % 2 == 0 else -1.0
            p0 = Vector(controls[-1]) - center_v
            x1 = side * next_width * 0.5
            y1 = lane_y(lane + 1, 0.0 if (lane + 1) % 2 == 1 else 1.0)
            p3 = Vector((x1, y1, 0.0))
            handle = max(8.0, lane_pitch * 0.62)
            p1 = p0 + Vector((side * handle, 0.0, 0.0))
            p2 = p3 + Vector((side * handle, 0.0, 0.0))
            for step in range(1, turn_samples):
                u = step / turn_samples
                omt = 1.0 - u
                point = p0 * (omt ** 3) + p1 * (3.0 * omt * omt * u) + p2 * (3.0 * omt * u * u) + p3 * (u ** 3)
                global_t = (lane + u) / lanes
                z = z_amp * 0.35 * math.sin(math.tau * (3.1 * global_t + phase))
                controls.append((center_v.x + point.x, center_v.y + point.y, center_v.z + z))
    return scale_path_to_length(controls, target_length)


def folded_ribbon_dna_controls(
    center: tuple[float, float, float],
    width: float,
    height: float,
    folds: float,
    samples: int,
    irregularity: float,
    z_amp: float,
    phase: float,
    target_length: float = CANONICAL_DNA_MM,
) -> list[tuple[float, float, float]]:
    center_v = Vector(center)
    controls = []
    bend_radius = width * 0.23
    for index in range(samples + 1):
        t = index / samples
        theta = math.tau * folds * t
        x = width * 0.36 * math.sin(theta) + bend_radius * 0.34 * math.sin(theta * 0.5 + phase)
        y = height * (t - 0.5)
        envelope = math.sin(math.pi * t) ** 0.50
        x += irregularity * envelope * (
            0.34 * math.sin(math.tau * (2.2 * t + phase))
            + 0.16 * math.sin(math.tau * (5.1 * t + 0.07))
        )
        y += irregularity * 0.40 * envelope * math.sin(math.tau * (3.4 * t + phase * 0.3))
        z = z_amp * envelope * math.sin(math.tau * (4.5 * t + phase))
        controls.append((center_v.x + x, center_v.y + y, center_v.z + z))
    return scale_path_to_length(controls, target_length)


def compact_loop_dna_controls(
    center: tuple[float, float, float],
    width: float,
    height: float,
    lobes: float,
    samples: int,
    irregularity: float,
    z_amp: float,
    phase: float,
    target_length: float = CANONICAL_DNA_MM,
) -> list[tuple[float, float, float]]:
    center_v = Vector(center)
    controls = []
    for index in range(samples + 1):
        t = index / samples
        angle = math.tau * lobes * t
        envelope = 0.72 + 0.16 * math.sin(math.tau * (3.0 * t + phase))
        x = width * 0.38 * envelope * math.sin(angle)
        y = height * 0.40 * math.sin(angle * 0.52 + 0.8) + height * 0.16 * math.sin(math.tau * (2.0 * t + phase))
        x += irregularity * math.sin(math.tau * (8.5 * t + phase)) * math.sin(math.pi * t)
        y += irregularity * 0.65 * math.sin(math.tau * (6.3 * t + 0.22)) * math.sin(math.pi * t)
        z = z_amp * math.sin(math.tau * (9.0 * t + phase)) * math.sin(math.pi * t)
        controls.append((center_v.x + x, center_v.y + y, center_v.z + z))
    return scale_path_to_length(controls, target_length)


def rna_spiral_points(
    start: Vector,
    vertical_rise: float,
    start_radius: float,
    end_radius: float,
    turns: float,
    samples: int,
    wiggle: float,
    phase: float,
) -> list[tuple[float, float, float]]:
    center0 = start - Vector((start_radius, 0.0, 0.0))
    points = []
    for index in range(samples + 1):
        t = index / samples
        ease = t * t * (3.0 - 2.0 * t)
        envelope = math.sin(math.pi * t) ** 0.72
        radius = start_radius * (1.0 - ease) + end_radius * ease
        angle = math.tau * turns * t
        radial_wiggle = wiggle * envelope * (
            0.50 * math.sin(math.tau * (3.7 * t + phase))
            + 0.30 * math.sin(math.tau * (9.4 * t + 0.27))
            + 0.20 * math.sin(math.tau * (17.0 * t + phase * 0.63))
        )
        lateral_x = wiggle * 0.90 * envelope * math.sin(math.tau * (1.3 * t + phase * 0.31))
        lateral_y = wiggle * 0.75 * envelope * math.sin(math.tau * (1.9 * t + phase))
        z_wiggle = wiggle * 0.55 * envelope * math.sin(math.tau * (5.5 * t + 0.18))
        center = center0 + Vector((lateral_x, lateral_y, vertical_rise * ease))
        radius_v = max(0.7, radius + radial_wiggle)
        point = center + Vector((math.cos(angle), math.sin(angle), 0.0)) * radius_v + Vector((0.0, 0.0, z_wiggle))
        points.append(point.to_tuple())
    points[0] = start.to_tuple()
    return points


def add_source_path_curve(
    name: str,
    path: geom.SampledPath,
    collection: bpy.types.Collection,
    material: bpy.types.Material,
) -> dict:
    points = [(point.x, point.y, point.z + 0.16) for point in path.points]
    obj = base.create_curve(name, points, 0.026, material, collection, resolution=1)
    obj.hide_render = True
    obj["source_path_for_rebake"] = True
    obj["curve_role"] = "editable_source_path"
    obj["sampled_points"] = len(points)
    return {"object": obj.name, "points": len(points), "hide_render": True}


def bounds_from_path(path: geom.SampledPath) -> dict:
    return geom.bounds(path.points)


def configure_direct_settings(manifest: dict) -> None:
    dna = manifest["procedural_nucleic_acids"]["dna"]
    dna.update(
        {
            "style": "direct_blender_surface_proxy_polished",
            "strand_sides": 12,
            "strand_radius_mm": 0.12,
            "base_pair_every_bp": 1,
            "base_pair_radius_mm": 0.082,
            "base_pair_sides": 10,
            "base_stack_radius_mm": 0.172,
            "base_stack_sides": 12,
            "surface_bump_amplitude": 0.05,
            "direct_voxel_size_mm": 0.052,
            "direct_smooth_factor": 0.035,
            "direct_smooth_iterations": 1,
            "direct_unified_mesh": True,
        }
    )
    # This experiment compares global DNA/RNA arrangements only. Protein-derived
    # nucleosome placement is intentionally left for the next attachment pass.
    dna.setdefault("nucleosome_loop", {})["enabled"] = False

    mrna = manifest["procedural_nucleic_acids"]["mrna"]
    mrna.update(
        {
            "style": "direct_blender_surface_proxy_polished",
            "tube_sides": 12,
            "tube_radius_mm": 0.165,
            "surface_bump_amplitude": 0.10,
            "base_lobe_every_nt": 1,
            "base_lobe_radius_mm": 0.095,
            "base_lobe_sides": 8,
            "base_offset_mm": 0.17,
            "nucleotide_detail_every_nt": 1,
            "phosphate_radius_mm": 0.09,
            "sugar_radius_mm": 0.085,
            "rna_base_ellipsoid_radii_mm": [0.135, 0.078, 0.055],
            "base_connector_radius_mm": 0.048,
            "nucleotide_detail_rings": 5,
            "nucleotide_detail_sides": 8,
            "direct_voxel_size_mm": 0.05,
            "direct_smooth_factor": 0.028,
            "direct_smooth_iterations": 1,
            "direct_unified_mesh": True,
        }
    )


def variant_specs() -> list[dict]:
    return [
        {
            "key": "compact_arc",
            "title": "Compact irregular 3/4 base",
            "panel_center": (-285.0, 315.0, 0.0),
            "dna_kind": "arc",
            "dna": {"radius": 44.0, "start_deg": 222.0, "span_deg": -268.0, "samples": 46, "radial_amp": 4.4, "z_amp": 0.45, "phase": 0.08},
            "rna_fraction": 0.96,
            "rna": {"vertical_rise": 86.0, "start_radius": 13.5, "end_radius": 5.8, "turns": 2.05, "samples": 1600, "wiggle": 2.8, "phase": 0.12},
            "intent": "more compact version of arrangement_v1 with DNA as an irregular partial circular base",
        },
        {
            "key": "dna_spiral_pedestal",
            "title": "Full-gene DNA spiral pedestal",
            "panel_center": (285.0, 315.0, 0.0),
            "dna_kind": "spiral",
            "dna_target_bp": ACTB_WITH_PROMOTER_BP,
            "dna_annotation": {
                "basis": "human ACTB canonical genomic transcript span plus upstream promoter",
                "actb_canonical_gene_bp": ACTB_CANONICAL_GENE_BP,
                "upstream_promoter_bp": UPSTREAM_PROMOTER_BP,
                "total_bp": ACTB_WITH_PROMOTER_BP,
            },
            "dna": {"start_radius": 48.0, "end_radius": 8.0, "turns": -3.25, "start_deg": 215.0, "samples": 96, "radial_amp": 3.2, "z_amp": 0.65, "phase": 0.28},
            "rna_anchor": "center",
            "rna": {"vertical_rise": 96.0, "start_radius": 8.0, "end_radius": 3.8, "turns": 2.05, "samples": 1600, "wiggle": 3.2, "phase": 0.39},
            "intent": "compact full-length ACTB-plus-promoter DNA spiral base with RNA emerging from the most central DNA point",
        },
        {
            "key": "nested_wiggle",
            "title": "Nested DNA coil, wiggly RNA",
            "panel_center": (-285.0, 95.0, 0.0),
            "dna_kind": "nested",
            "dna": {"radius": 58.0, "samples": 64, "phase": 0.18},
            "rna_fraction": 0.78,
            "rna": {"vertical_rise": 96.0, "start_radius": 14.0, "end_radius": 6.5, "turns": 1.85, "samples": 1700, "wiggle": 6.0, "phase": 0.55},
            "intent": "nested irregular DNA base with a visibly more wandering upward RNA path",
        },
        {
            "key": "broad_base_tight_rna",
            "title": "Broad base, tight RNA spiral",
            "panel_center": (285.0, 95.0, 0.0),
            "dna_kind": "arc",
            "dna": {"radius": 54.0, "start_deg": 236.0, "span_deg": -252.0, "samples": 52, "radial_amp": 5.8, "z_amp": 0.55, "phase": 0.42},
            "rna_fraction": 0.88,
            "rna": {"vertical_rise": 98.0, "start_radius": 8.5, "end_radius": 3.8, "turns": 3.2, "samples": 1750, "wiggle": 2.6, "phase": 0.71},
            "intent": "broader organic DNA footprint with a tighter vertical RNA spiral",
        },
        {
            "key": "full_gene_serpentine_base",
            "title": "Full-gene serpentine base",
            "panel_center": (-285.0, -125.0, 0.0),
            "dna_kind": "serpentine",
            "dna_target_bp": ACTB_WITH_PROMOTER_BP,
            "dna_annotation": {
                "basis": "human ACTB canonical genomic transcript span plus upstream promoter",
                "actb_canonical_gene_bp": ACTB_CANONICAL_GENE_BP,
                "upstream_promoter_bp": UPSTREAM_PROMOTER_BP,
                "total_bp": ACTB_WITH_PROMOTER_BP,
            },
            "dna": {"width": 104.0, "height": 66.0, "lanes": 6, "samples_per_lane": 14, "irregularity": 3.8, "z_amp": 0.42, "phase": 0.16},
            "rna_anchor": "center",
            "rna": {"vertical_rise": 94.0, "start_radius": 8.5, "end_radius": 4.0, "turns": 2.35, "samples": 1650, "wiggle": 3.5, "phase": 0.22},
            "intent": "full-length ACTB-plus-promoter DNA folded back and forth as a compact horizontal standing base",
        },
        {
            "key": "full_gene_irregular_ribbon",
            "title": "Full-gene irregular ribbon",
            "panel_center": (285.0, -125.0, 0.0),
            "dna_kind": "folded_ribbon",
            "dna_target_bp": ACTB_WITH_PROMOTER_BP,
            "dna_annotation": {
                "basis": "human ACTB canonical genomic transcript span plus upstream promoter",
                "actb_canonical_gene_bp": ACTB_CANONICAL_GENE_BP,
                "upstream_promoter_bp": UPSTREAM_PROMOTER_BP,
                "total_bp": ACTB_WITH_PROMOTER_BP,
            },
            "dna": {"width": 92.0, "height": 92.0, "folds": 4.35, "samples": 132, "irregularity": 4.4, "z_amp": 0.48, "phase": 0.44},
            "rna_anchor": "center",
            "rna": {"vertical_rise": 98.0, "start_radius": 9.0, "end_radius": 4.0, "turns": 2.15, "samples": 1650, "wiggle": 4.6, "phase": 0.58},
            "intent": "full-length ACTB-plus-promoter DNA as an irregular back-and-forth ribbon base",
        },
        {
            "key": "full_gene_compact_loop",
            "title": "Full-gene compact loop base",
            "panel_center": (-285.0, -345.0, 0.0),
            "dna_kind": "compact_loop",
            "dna_target_bp": ACTB_WITH_PROMOTER_BP,
            "dna_annotation": {
                "basis": "human ACTB canonical genomic transcript span plus upstream promoter",
                "actb_canonical_gene_bp": ACTB_CANONICAL_GENE_BP,
                "upstream_promoter_bp": UPSTREAM_PROMOTER_BP,
                "total_bp": ACTB_WITH_PROMOTER_BP,
            },
            "dna": {"width": 114.0, "height": 82.0, "lobes": 3.7, "samples": 148, "irregularity": 3.6, "z_amp": 0.50, "phase": 0.31},
            "rna_anchor": "center",
            "rna": {"vertical_rise": 92.0, "start_radius": 8.0, "end_radius": 3.8, "turns": 2.4, "samples": 1650, "wiggle": 4.0, "phase": 0.69},
            "intent": "full-length ACTB-plus-promoter DNA in a compact non-spiral looping base",
        },
        {
            "key": "full_gene_varied_serpentine",
            "title": "Full-gene varied serpentine",
            "panel_center": (285.0, -345.0, 0.0),
            "dna_kind": "serpentine",
            "dna_target_bp": ACTB_WITH_PROMOTER_BP,
            "dna_annotation": {
                "basis": "human ACTB canonical genomic transcript span plus upstream promoter",
                "actb_canonical_gene_bp": ACTB_CANONICAL_GENE_BP,
                "upstream_promoter_bp": UPSTREAM_PROMOTER_BP,
                "total_bp": ACTB_WITH_PROMOTER_BP,
            },
            "dna": {"width": 112.0, "height": 74.0, "lanes": 6, "samples_per_lane": 16, "irregularity": 3.2, "z_amp": 0.42, "phase": 0.53, "length_variation": 0.23},
            "rna_anchor": "center",
            "rna": {"vertical_rise": 96.0, "start_radius": 8.2, "end_radius": 3.9, "turns": 2.25, "samples": 1650, "wiggle": 3.9, "phase": 0.47},
            "intent": "full-length ACTB-plus-promoter DNA as a rounded serpentine base with unequal lane lengths",
        },
    ]


def make_dna_controls(spec: dict) -> list[tuple[float, float, float]]:
    center = spec["panel_center"]
    dna = spec["dna"]
    target_length = float(spec.get("dna_target_bp", TARGET_DNA_BP)) * 0.136
    if spec["dna_kind"] == "spiral":
        return spiral_dna_controls(center=center, target_length=target_length, **dna)
    if spec["dna_kind"] == "nested":
        return nested_coil_controls(center=center, target_length=target_length, **dna)
    if spec["dna_kind"] == "serpentine":
        return serpentine_dna_controls(center=center, target_length=target_length, **dna)
    if spec["dna_kind"] == "folded_ribbon":
        return folded_ribbon_dna_controls(center=center, target_length=target_length, **dna)
    if spec["dna_kind"] == "compact_loop":
        return compact_loop_dna_controls(center=center, target_length=target_length, **dna)
    return organic_arc_controls(center=center, target_length=target_length, **dna)


def build_variant(base_manifest: dict, spec: dict, materials: dict, labels: bpy.types.Collection) -> dict:
    master = variant_collection(f"Variant {spec['key']}")
    dna_collection = variant_collection(f"{spec['key']} DNA", master)
    mrna_collection = variant_collection(f"{spec['key']} mRNA", master)
    source_collection = variant_collection(f"{spec['key']} editable source paths", master)

    manifest = copy.deepcopy(base_manifest)
    configure_direct_settings(manifest)
    dna_controls = make_dna_controls(spec)
    if spec.get("rna_anchor") == "center":
        rna_anchor_point, rna_fraction = centered_sampled_point(dna_controls, spec["panel_center"])
    else:
        rna_fraction = spec["rna_fraction"]
        rna_anchor_point = sampled_point(dna_controls, rna_fraction)
    rna_start = rna_anchor_point + Vector((0.0, 0.0, 2.0))
    rna_points = rna_spiral_points(rna_start, **spec["rna"])

    manifest["title"] = f"{base_manifest['title']} - arrangement variant {spec['key']}"
    manifest["procedural_nucleic_acids"]["dna"]["path_mode"] = "custom_dna_controls"
    manifest["procedural_nucleic_acids"]["dna"]["custom_controls_mm"] = dna_controls
    manifest["procedural_nucleic_acids"]["mrna"]["path_mode"] = "custom_mrna_path"
    manifest["procedural_nucleic_acids"]["mrna"]["custom_points_mm"] = rna_points
    manifest["mrna"]["start_mm"] = list(rna_start)

    collections = {"DNA": dna_collection, "mRNA": mrna_collection}
    dna_path, dna_report = direct_meshes.build_dna_meshes(manifest, collections, materials)
    mrna_path, mrna_report = direct_meshes.build_mrna_meshes(manifest, collections, materials)
    dna_curve = add_source_path_curve(f"{spec['key']}_DNA_source_path", dna_path, source_collection, materials["black"])
    mrna_curve = add_source_path_curve(f"{spec['key']}_MRNA_source_path", mrna_path, source_collection, materials["black"])

    center = spec["panel_center"]
    label_y = center[1] - 74.0 if center[1] > 0.0 else center[1] + 74.0
    title_label = base.create_text(f"label_{spec['key']}_title", spec["title"], (center[0], label_y, 1.0), 4.0, materials["black"], labels)
    length_label = base.create_text(
        f"label_{spec['key']}_lengths",
        f"DNA {dna_report['represented_bp']} bp | RNA 1852 nt",
        (center[0], label_y - 8.0, 1.0),
        2.6,
        materials["label_grey"],
        labels,
    )

    return {
        "key": spec["key"],
        "title": spec["title"],
        "intent": spec["intent"],
        "panel_center_mm": list(center),
        "rna_start_fraction_on_dna": rna_fraction,
        "rna_anchor": spec.get("rna_anchor", "fraction"),
        "target_dna_bp": int(spec.get("dna_target_bp", TARGET_DNA_BP)),
        "dna_annotation": spec.get("dna_annotation"),
        "dna": dna_report,
        "mrna": mrna_report,
        "source_curves": {"dna": dna_curve, "mrna": mrna_curve},
        "path_bounds": {"dna": bounds_from_path(dna_path), "mrna": bounds_from_path(mrna_path)},
        "label_objects": [title_label.name, length_label.name],
        "collection": master.name,
    }


def look_at(camera: bpy.types.Object, target: Vector) -> None:
    direction = target - camera.location
    camera.rotation_euler = direction.to_track_quat("-Z", "Y").to_euler()


def add_camera(name: str, location: tuple[float, float, float], target: tuple[float, float, float], ortho_scale: float) -> str:
    camera_data = bpy.data.cameras.new(name)
    camera = bpy.data.objects.new(name, camera_data)
    camera.location = location
    look_at(camera, Vector(target))
    camera_data.type = "ORTHO"
    camera_data.ortho_scale = ortho_scale
    bpy.context.scene.collection.objects.link(camera)
    return name


def add_cameras(variants: list[dict]) -> dict[str, str]:
    cameras = {
        "overview": add_camera("Camera_arrangement_variants_overview", (0.0, -680.0, 590.0), (0.0, -25.0, 46.0), 900.0)
    }
    for variant in variants:
        center = variant["panel_center_mm"]
        key = variant["key"]
        cameras[key] = add_camera(
            f"Camera_arrangement_variants_{key}",
            (center[0] - 78.0, center[1] - 132.0, 142.0),
            (center[0], center[1], 36.0),
            245.0,
        )
    bpy.context.scene.camera = bpy.data.objects[cameras["overview"]]
    return cameras


def canonical_dna_scale_notes() -> dict:
    current_bp = TARGET_DNA_BP
    dna_bp_to_mm = 0.136
    actb_canonical_gene_bp = 3454
    actb_exon_bp = 1812
    actb_ensembl_broad_bp = 37494
    return {
        "current_canonical_dna": {
            "represented_bp": current_bp,
            "axis_length_mm": CANONICAL_DNA_MM,
            "bp_spacing_mm": dna_bp_to_mm,
            "interpretation": "short/exon-scale actin-region DNA",
        },
        "current_canonical_mrna": {"represented_nt": 1852, "contour_length_mm": TARGET_MRNA_MM},
        "human_actb_default_assumption": {
            "gene_symbol": "ACTB",
            "ncbi_gene_id": "60",
            "canonical_genomic_transcript_span_bp": actb_canonical_gene_bp,
            "ensembl_canonical_exon_total_bp": actb_exon_bp,
            "ensembl_broad_gene_span_bp": actb_ensembl_broad_bp,
        },
        "scale_comparison": {
            "current_dna_fraction_of_actb_canonical_span": current_bp / actb_canonical_gene_bp,
            "current_dna_fraction_of_actb_exon_total": current_bp / actb_exon_bp,
            "actb_canonical_span_full_scale_mm": actb_canonical_gene_bp * dna_bp_to_mm,
            "actb_canonical_span_1_to_10_mm": actb_canonical_gene_bp * dna_bp_to_mm / 10.0,
            "ensembl_broad_span_full_scale_mm": actb_ensembl_broad_bp * dna_bp_to_mm,
            "ensembl_broad_span_1_to_10_mm": actb_ensembl_broad_bp * dna_bp_to_mm / 10.0,
        },
        "sources": [
            "https://www.ncbi.nlm.nih.gov/gene/60",
            "https://rest.ensembl.org/lookup/symbol/homo_sapiens/ACTB?expand=1;content-type=application/json",
        ],
    }


def validate_report(report: dict) -> None:
    variants = report["variants"]
    if len(variants) < 4:
        raise RuntimeError("Arrangement variants experiment must contain at least four variants")
    z_extents = []
    dna_lengths = []
    for variant in variants:
        dna = variant["dna"]
        mrna = variant["mrna"]
        expected_dna_mm = float(variant["target_dna_bp"]) * 0.136
        if dna.get("generation_pipeline") != "direct_blender_mesh":
            raise RuntimeError(f"{variant['key']} DNA was not generated by the direct Blender pipeline")
        if mrna.get("generation_pipeline") != "direct_blender_mesh":
            raise RuntimeError(f"{variant['key']} RNA was not generated by the direct Blender pipeline")
        if abs(float(mrna["total_measured_mm"]) - TARGET_MRNA_MM) > 0.05:
            raise RuntimeError(f"{variant['key']} mRNA length is not {TARGET_MRNA_MM} mm")
        if not variant["source_curves"]["dna"].get("hide_render") or not variant["source_curves"]["mrna"].get("hide_render"):
            raise RuntimeError(f"{variant['key']} source curves are not marked render-hidden")
        z_extents.append(float(variant["path_bounds"]["mrna"]["bbox_mm"][2]))
        dna_lengths.append((variant["key"], float(dna["axis_length_mm"]), expected_dna_mm))
    if min(z_extents) < 68.0:
        raise RuntimeError("One or more RNA variants does not rise enough in Z")
    if max(z_extents) - min(z_extents) < 10.0:
        raise RuntimeError("RNA Z extents are too similar for arrangement comparison")
    for key, length, expected in dna_lengths:
        if abs(length - expected) > 1.25:
            raise RuntimeError(f"{key} DNA drifted too far from its target contour length")


def set_label_render_policy(variants: list[dict], current_variant_key: str | None = None) -> None:
    for variant in variants:
        visible = current_variant_key is None or variant["key"] == current_variant_key
        for name in variant.get("label_objects", []):
            obj = bpy.data.objects.get(name)
            if obj:
                obj.hide_render = not visible


def collection_objects_recursive(collection: bpy.types.Collection):
    for obj in collection.objects:
        yield obj
    for child in collection.children:
        yield from collection_objects_recursive(child)


def set_variant_render_policy(variants: list[dict], current_variant_key: str | None = None) -> None:
    for variant in variants:
        collection = bpy.data.collections.get(variant["collection"])
        if not collection:
            continue
        visible = current_variant_key is None or variant["key"] == current_variant_key
        for obj in collection_objects_recursive(collection):
            obj.hide_render = True if obj.get("source_path_for_rebake") else not visible


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    base_manifest = base.load_manifest()

    surface_scene.clean_scene()
    surface_scene.configure_scene()
    materials = base.make_materials()
    surface_scene.soften_materials(materials)
    labels = variant_collection("Labels")

    variants = [build_variant(base_manifest, spec, materials, labels) for spec in variant_specs()]
    cameras = add_cameras(variants)
    base.create_text("label_arrangement_variants_title", "DNA/RNA arrangement variants", (0.0, 475.0, 1.0), 5.2, materials["black"], labels)
    base.create_text(
        "label_arrangement_variants_scope",
        "direct Blender DNA/RNA only; protein attachment deferred",
        (0.0, 464.0, 1.0),
        3.0,
        materials["label_grey"],
        labels,
    )

    preview_paths = {
        key: OUTPUT_DIR / VARIANT_PREVIEW_TEMPLATE.format(key=key)
        for key in cameras
        if key != "overview"
    }
    report = {
        "title": "DNA/RNA arrangement variants experiment",
        "kind": "arrangement_variants_experiment",
        "units": base_manifest["units"],
        "outputs": {
            "blend": str(BLEND_PATH),
            "preview": str(PREVIEW_PATH),
            "report": str(REPORT_PATH),
            "variant_previews": {key: str(path) for key, path in preview_paths.items()},
        },
        "pipeline": {
            "nucleic_acid_generation": "direct_blender_meshes_for_dna_and_mrna",
            "pdb_asset_imports": False,
            "protein_attachment": "deferred",
            "variants_are_display_scaled": False,
        },
        "canonical_dna_scale_notes": canonical_dna_scale_notes(),
        "variants": variants,
    }
    validate_report(report)

    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_PATH))
    set_variant_render_policy(variants, None)
    set_label_render_policy(variants, None)
    bpy.context.scene.camera = bpy.data.objects[cameras["overview"]]
    bpy.context.scene.render.filepath = str(PREVIEW_PATH)
    bpy.ops.render.render(write_still=True)
    for key, path in preview_paths.items():
        set_variant_render_policy(variants, key)
        set_label_render_policy(variants, key)
        bpy.context.scene.camera = bpy.data.objects[cameras[key]]
        bpy.context.scene.render.filepath = str(path)
        bpy.ops.render.render(write_still=True)
    set_variant_render_policy(variants, None)
    set_label_render_policy(variants, None)
    print(f"Wrote {BLEND_PATH}")
    print(f"Wrote {PREVIEW_PATH}")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
