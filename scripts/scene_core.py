#!/usr/bin/env python3
"""Shared Blender scene primitives used by canonical and retained experiments."""

from __future__ import annotations

import json
import math
import os
import shlex
from collections import defaultdict
from pathlib import Path

import bpy
from mathutils import Euler, Matrix, Vector


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "config" / "scene_manifest.json"
ASSET_DIR = ROOT / "assets" / "rcsb"

NUCLEIC_COMPS = {
    "A",
    "C",
    "G",
    "U",
    "I",
    "DA",
    "DC",
    "DG",
    "DT",
    "DI",
    "ADE",
    "CYT",
    "GUA",
    "THY",
    "URA",
    # Common modified RNA bases, especially in tRNA structures such as 1EHZ.
    "1MA",
    "2MG",
    "5MC",
    "5MU",
    "7MG",
    "8AN",
    "H2U",
    "M2G",
    "OMC",
    "OMG",
    "PSU",
    "YYG",
}
WATER_COMPS = {"HOH", "WAT", "DOD"}

COLORS = {
    "dna_orange": (0.92, 0.45, 0.18, 1.0),
    "dna_dark": (0.58, 0.25, 0.08, 1.0),
    "dna_promoter_blue": (0.18, 0.42, 0.86, 1.0),
    "dna_exon_orange": (0.94, 0.48, 0.17, 1.0),
    "dna_intron_olive": (0.48, 0.55, 0.24, 1.0),
    "olive": (0.63, 0.58, 0.22, 1.0),
    "orange": (0.86, 0.43, 0.18, 1.0),
    "yellow_olive": (0.74, 0.68, 0.32, 1.0),
    "rna_gold": (0.80, 0.66, 0.30, 1.0),
    "rna_red": (0.88, 0.24, 0.22, 1.0),
    "ribosome_red": (0.88, 0.32, 0.38, 1.0),
    "ribosome_blue": (0.30, 0.43, 0.74, 1.0),
    "pol_green": (0.22, 0.64, 0.26, 1.0),
    "histone_green": (0.42, 0.70, 0.25, 1.0),
    "cas9_blue": (0.17, 0.58, 0.78, 1.0),
    "guide_orange": (0.95, 0.48, 0.22, 1.0),
    "actin_blue": (0.36, 0.48, 0.78, 1.0),
    "tf_cyan": (0.22, 0.72, 0.70, 1.0),
    "tf_lavender": (0.52, 0.50, 0.84, 1.0),
    "tf_purple": (0.38, 0.31, 0.78, 1.0),
    "ago_pink": (0.84, 0.34, 0.58, 1.0),
    "rbp_purple": (0.58, 0.36, 0.78, 1.0),
    "mcp_blue": (0.34, 0.54, 0.74, 1.0),
    "rfp_red": (0.90, 0.22, 0.22, 1.0),
    "label_grey": (0.43, 0.43, 0.43, 1.0),
    "scale_grey": (0.55, 0.55, 0.55, 1.0),
    "white": (0.92, 0.92, 0.90, 1.0),
    "black": (0.05, 0.05, 0.05, 1.0),
}


def load_manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def clean_scene() -> None:
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete()
    for block in list(bpy.data.meshes):
        if block.users == 0:
            bpy.data.meshes.remove(block)
    for block in list(bpy.data.curves):
        if block.users == 0:
            bpy.data.curves.remove(block)


def configure_scene() -> None:
    scene = bpy.context.scene
    scene.unit_settings.system = "METRIC"
    scene.unit_settings.length_unit = "MILLIMETERS"
    scene.unit_settings.scale_length = 0.001
    scene.render.engine = "BLENDER_EEVEE"
    scene.eevee.use_gtao = False
    scene.render.resolution_x = 2600
    scene.render.resolution_y = 1600
    scene.view_settings.view_transform = "Standard"
    scene.view_settings.look = "None"
    scene.view_settings.exposure = 0.0
    scene.view_settings.gamma = 1.0
    scene.world = bpy.data.worlds.new("World") if scene.world is None else scene.world
    scene.world.color = (1.0, 1.0, 1.0)
    scene.world.use_nodes = True
    background = scene.world.node_tree.nodes.get("Background")
    if background:
        background.inputs["Color"].default_value = (1.0, 1.0, 1.0, 1.0)
        background.inputs["Strength"].default_value = 1.0


def make_collections(names: list[str]) -> dict[str, bpy.types.Collection]:
    collections = {}
    for name in names:
        collection = bpy.data.collections.new(name)
        bpy.context.scene.collection.children.link(collection)
        collections[name] = collection
    return collections


def make_materials() -> dict[str, bpy.types.Material]:
    materials = {}
    for name, color in COLORS.items():
        mat = bpy.data.materials.new(name)
        mat.diffuse_color = color
        mat.use_nodes = True
        bsdf = mat.node_tree.nodes.get("Principled BSDF")
        if bsdf:
            bsdf.inputs["Base Color"].default_value = color
            bsdf.inputs["Roughness"].default_value = 0.72
            bsdf.inputs["Metallic"].default_value = 0.0
            if "Emission" in bsdf.inputs:
                bsdf.inputs["Emission"].default_value = color
            if "Emission Strength" in bsdf.inputs:
                bsdf.inputs["Emission Strength"].default_value = 0.12
        materials[name] = mat
    return materials


def link_to_collection(obj: bpy.types.Object, collection: bpy.types.Collection) -> None:
    collection.objects.link(obj)
    for coll in list(obj.users_collection):
        if coll != collection:
            coll.objects.unlink(obj)


def create_curve(
    name: str,
    points: list[tuple[float, float, float]],
    bevel_depth: float,
    material: bpy.types.Material,
    collection: bpy.types.Collection,
    resolution: int = 3,
) -> bpy.types.Object:
    curve = bpy.data.curves.new(name, "CURVE")
    curve.dimensions = "3D"
    curve.resolution_u = resolution
    curve.bevel_depth = bevel_depth
    curve.bevel_resolution = 4
    spline = curve.splines.new("POLY")
    spline.points.add(len(points) - 1)
    for point, co in zip(spline.points, points):
        point.co = (co[0], co[1], co[2], 1.0)
    curve.materials.append(material)
    obj = bpy.data.objects.new(name, curve)
    link_to_collection(obj, collection)
    return obj


def create_text(
    name: str,
    text: str,
    location: tuple[float, float, float],
    size: float,
    material: bpy.types.Material,
    collection: bpy.types.Collection,
    align: str = "CENTER",
) -> bpy.types.Object:
    curve = bpy.data.curves.new(name, "FONT")
    curve.body = text
    curve.align_x = align
    curve.align_y = "CENTER"
    curve.size = size
    curve.materials.append(material)
    obj = bpy.data.objects.new(name, curve)
    obj.location = location
    link_to_collection(obj, collection)
    return obj


def polyline_length(points: list[tuple[float, float, float]]) -> float:
    return sum(
        math.dist(points[i - 1], points[i])
        for i in range(1, len(points))
    )


def sine_arc_length(span: float, amplitude: float, period: float, phase: float, steps: int = 600) -> float:
    if span <= 0:
        return 0.0
    total = 0.0
    previous = (0.0, amplitude * math.sin(phase))
    for i in range(1, steps + 1):
        x = span * i / steps
        y = amplitude * math.sin((2 * math.pi * x / period) + phase)
        total += math.hypot(x - previous[0], y - previous[1])
        previous = (x, y)
    return total


def span_for_arc_length(target: float, amplitude: float, period: float, phase: float) -> float:
    low = 0.0
    high = target
    for _ in range(50):
        mid = (low + high) / 2.0
        if sine_arc_length(mid, amplitude, period, phase) < target:
            low = mid
        else:
            high = mid
    return high


def segment_points_for_span(
    start_point: tuple[float, float, float],
    base_y: float,
    amplitude: float,
    period: float,
    phase: float,
    span: float,
    point_count: int,
) -> list[tuple[float, float, float]]:
    points = [start_point]
    start_x = start_point[0]
    z = start_point[2]
    for i in range(1, point_count + 1):
        dx = span * i / point_count
        local_phase = (2 * math.pi * dx / period) + phase
        points.append((start_x + dx, base_y + amplitude * math.sin(local_phase), z))
    return points


def segment_length_for_span(
    start_point: tuple[float, float, float],
    base_y: float,
    amplitude: float,
    period: float,
    phase: float,
    span: float,
    steps: int = 600,
) -> float:
    points = segment_points_for_span(start_point, base_y, amplitude, period, phase, span, steps)
    return polyline_length(points)


def span_for_segment_length(
    target: float,
    start_point: tuple[float, float, float],
    base_y: float,
    amplitude: float,
    period: float,
    phase: float,
) -> float:
    low = 0.0
    high = max(target, 1.0)
    while segment_length_for_span(start_point, base_y, amplitude, period, phase, high) < target:
        high *= 2.0
    for _ in range(50):
        mid = (low + high) / 2.0
        if segment_length_for_span(start_point, base_y, amplitude, period, phase, mid) < target:
            low = mid
        else:
            high = mid
    return high


def build_mrna(manifest: dict, collections: dict, materials: dict) -> dict:
    units = manifest["units"]
    mrna = manifest["mrna"]
    nt_to_mm = units["mrna_nt_to_mm"]
    amp = mrna["wave_amplitude_mm"]
    period = mrna["wave_period_mm"]
    x, y0, z0 = mrna["start_mm"]
    phase = 0.7
    segment_reports = []
    segment_midpoints = []
    last_point = (x, y0 + amp * math.sin(phase), z0)

    for segment in mrna["segments"]:
        base_y = segment.get("track_y_mm", y0)
        target_len = segment["nt"] * nt_to_mm
        span = span_for_segment_length(target_len, last_point, base_y, amp, period, phase)
        point_count = max(30, int(target_len * 5))
        points = segment_points_for_span(last_point, base_y, amp, period, phase, span, point_count)
        obj = create_curve(
            f"mRNA_{segment['name']}",
            points,
            mrna["tube_radius_mm"],
            materials[segment["color"]],
            collections["mRNA"],
        )
        obj["nt"] = segment["nt"]
        obj["intended_contour_length_mm"] = target_len
        measured = polyline_length(points)
        midpoint = points[len(points) // 2]
        segment_midpoints.append((segment["name"], segment["nt"], midpoint))
        segment_reports.append(
            {
                "name": segment["name"],
                "nt": segment["nt"],
                "track_y_mm": base_y,
                "target_length_mm": target_len,
                "measured_length_mm": measured,
                "error_mm": measured - target_len,
                "start_mm": list(points[0]),
                "end_mm": list(points[-1]),
            }
        )
        x = points[-1][0]
        phase += 2 * math.pi * span / period
        last_point = points[-1]

    create_text(
        "label_mRNA_total",
        "actin mRNA: 1852 nt",
        (-45.0, 61.0, 0.2),
        4.0,
        materials["black"],
        collections["Labels"],
    )
    for name, nt, midpoint in segment_midpoints:
        segment = next(item for item in mrna["segments"] if item["name"] == name)
        label_location = segment.get("label_location_mm", (midpoint[0], midpoint[1] + 7.0, 0.3))
        create_text(
            f"label_{name}",
            f"{name} ({nt} nt)",
            tuple(label_location),
            2.3,
            materials["label_grey"],
            collections["Labels"],
        )
    return {"segments": segment_reports, "total_measured_mm": sum(item["measured_length_mm"] for item in segment_reports)}


def build_dna(manifest: dict, collections: dict, materials: dict) -> dict:
    dna = manifest["dna"]
    units = manifest["units"]
    start_x = dna["start_x_mm"]
    end_x = dna["end_x_mm"]
    y0 = dna["y_mm"]
    z0 = dna["center_z_mm"]
    span = end_x - start_x
    bp_to_mm = units["dna_bp_to_mm"]
    pitch = bp_to_mm * 10.5
    helix_radius = dna["helix_radius_mm"]
    step = 0.16
    n = int(span / step)
    strand_a = []
    strand_b = []
    for i in range(n + 1):
        x = start_x + span * i / n
        theta = 2 * math.pi * (x - start_x) / pitch
        strand_a.append((x, y0 + helix_radius * math.sin(theta), z0 + helix_radius * math.cos(theta)))
        strand_b.append((x, y0 + helix_radius * math.sin(theta + math.pi), z0 + helix_radius * math.cos(theta + math.pi)))
    create_curve("DNA_strand_A", strand_a, dna["strand_radius_mm"], materials["dna_orange"], collections["DNA"])
    create_curve("DNA_strand_B", strand_b, dna["strand_radius_mm"], materials["dna_dark"], collections["DNA"])

    rung_step = bp_to_mm * dna["rung_every_bp"]
    rung_count = int(span / rung_step)
    for i in range(rung_count + 1):
        x = start_x + i * rung_step
        theta = 2 * math.pi * (x - start_x) / pitch
        a = (x, y0 + helix_radius * math.sin(theta), z0 + helix_radius * math.cos(theta))
        b = (x, y0 + helix_radius * math.sin(theta + math.pi), z0 + helix_radius * math.cos(theta + math.pi))
        create_curve(
            f"DNA_basepair_rung_{i:03d}",
            [a, b],
            dna["rung_radius_mm"],
            materials["dna_orange"],
            collections["DNA"],
            resolution=1,
        )
    create_text("label_DNA", "DNA", (-100.0, -37.5, 0.2), 4.0, materials["black"], collections["Labels"])
    return {
        "axis_length_mm": span,
        "represented_bp": span / bp_to_mm,
        "helix_pitch_mm": pitch,
        "strand_a_length_mm": polyline_length(strand_a),
    }


def build_scale_bars(manifest: dict, collections: dict, materials: dict) -> dict:
    units = manifest["units"]
    reports = {}
    bp100 = 100 * units["dna_bp_to_mm"]
    nt100 = 100 * units["mrna_nt_to_mm"]
    for name, length, origin, label in [
        ("scale_100_bp", bp100, (-103.0, -56.5, 0.0), "100 bp"),
        ("scale_100_nt", nt100, (28.0, 20.0, 0.0), "100 nt"),
    ]:
        x, y, z = origin
        create_curve(name, [(x, y, z), (x + length, y, z)], 0.07, materials["scale_grey"], collections["Scale bars"])
        create_curve(f"{name}_left_tick", [(x, y - 0.8, z), (x, y + 0.8, z)], 0.055, materials["scale_grey"], collections["Scale bars"])
        create_curve(f"{name}_right_tick", [(x + length, y - 0.8, z), (x + length, y + 0.8, z)], 0.055, materials["scale_grey"], collections["Scale bars"])
        create_text(f"label_{name}", label, (x + length / 2.0, y - 2.3, z), 1.9, materials["label_grey"], collections["Labels"])
        reports[name] = {"length_mm": length}
    return reports


def parse_atom_site(path: Path) -> list[dict]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    atoms = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line != "loop_":
            i += 1
            continue
        i += 1
        cols = []
        while i < len(lines) and lines[i].strip().startswith("_"):
            cols.append(lines[i].strip())
            i += 1
        if not cols or not cols[0].startswith("_atom_site."):
            while i < len(lines) and lines[i].strip() != "#":
                i += 1
            continue
        index = {col: n for n, col in enumerate(cols)}
        required = ["_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z"]
        if not all(col in index for col in required):
            continue
        while i < len(lines):
            raw = lines[i].strip()
            if not raw:
                i += 1
                continue
            if raw == "#" or raw == "loop_" or raw.startswith("_"):
                break
            parts = raw.split()
            if len(parts) != len(cols):
                parts = shlex.split(raw)
            if len(parts) >= len(cols):
                def value(*names: str) -> str:
                    for name in names:
                        n = index.get(name)
                        if n is not None and n < len(parts):
                            return parts[n]
                    return ""

                model = value("_atom_site.pdbx_PDB_model_num")
                alt = value("_atom_site.label_alt_id")
                comp = value("_atom_site.auth_comp_id", "_atom_site.label_comp_id").strip("'\"").upper()
                if model and model not in {"1", ".", "?"}:
                    i += 1
                    continue
                if alt not in {"", ".", "?", "A", "1"}:
                    i += 1
                    continue
                if comp in WATER_COMPS:
                    i += 1
                    continue
                try:
                    atoms.append(
                        {
                            "group": value("_atom_site.group_PDB"),
                            "element": value("_atom_site.type_symbol").strip("'\"").upper(),
                            "atom": value("_atom_site.auth_atom_id", "_atom_site.label_atom_id").strip("'\""),
                            "comp": comp,
                            "chain": value("_atom_site.auth_asym_id", "_atom_site.label_asym_id").strip("'\""),
                            "seq": value("_atom_site.auth_seq_id", "_atom_site.label_seq_id").strip("'\""),
                            "x": float(value("_atom_site.Cartn_x")),
                            "y": float(value("_atom_site.Cartn_y")),
                            "z": float(value("_atom_site.Cartn_z")),
                        }
                    )
                except ValueError:
                    pass
            i += 1
    return atoms


def residue_points(atoms: list[dict]) -> list[dict]:
    grouped = defaultdict(list)
    ligand_counter = 0
    for atom in atoms:
        seq = atom["seq"]
        if seq in {"", ".", "?"}:
            ligand_counter += 1
            seq = f"ligand_{ligand_counter}"
        grouped[(atom["chain"], seq, atom["comp"])].append(atom)

    points = []
    for (chain, seq, comp), items in grouped.items():
        is_nucleic = comp in NUCLEIC_COMPS
        preferred = ["P", "C4'", "C3'", "C1'"] if is_nucleic else ["CA", "C4'", "P"]
        chosen = None
        atom_by_name = {item["atom"]: item for item in items}
        for atom_name in preferred:
            if atom_name in atom_by_name:
                chosen = atom_by_name[atom_name]
                break
        if chosen is None:
            x = sum(item["x"] for item in items) / len(items)
            y = sum(item["y"] for item in items) / len(items)
            z = sum(item["z"] for item in items) / len(items)
        else:
            x, y, z = chosen["x"], chosen["y"], chosen["z"]
        kind = "nucleic" if is_nucleic else "protein"
        points.append({"chain": chain, "seq": seq, "comp": comp, "kind": kind, "pos_A": (x, y, z)})
    return points


def sample_points(points: list[dict], max_beads: int) -> list[dict]:
    if len(points) <= max_beads:
        return points
    step = len(points) / float(max_beads)
    sampled = []
    seen = set()
    for i in range(max_beads):
        idx = min(len(points) - 1, int(round(i * step)))
        if idx not in seen:
            sampled.append(points[idx])
            seen.add(idx)
    return sampled


def icosahedron() -> tuple[list[Vector], list[tuple[int, int, int]]]:
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    verts = [
        Vector((-1, phi, 0)),
        Vector((1, phi, 0)),
        Vector((-1, -phi, 0)),
        Vector((1, -phi, 0)),
        Vector((0, -1, phi)),
        Vector((0, 1, phi)),
        Vector((0, -1, -phi)),
        Vector((0, 1, -phi)),
        Vector((phi, 0, -1)),
        Vector((phi, 0, 1)),
        Vector((-phi, 0, -1)),
        Vector((-phi, 0, 1)),
    ]
    verts = [v.normalized() for v in verts]
    faces = [
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
    ]
    return verts, faces


def create_bead_mesh(
    name: str,
    centers: list[Vector],
    radius: float,
    material: bpy.types.Material,
    collection: bpy.types.Collection,
) -> bpy.types.Object | None:
    if not centers:
        return None
    base_verts, base_faces = icosahedron()
    verts = []
    faces = []
    for center in centers:
        offset = len(verts)
        verts.extend([tuple(center + v * radius) for v in base_verts])
        faces.extend([(a + offset, b + offset, c + offset) for a, b, c in base_faces])
    mesh = bpy.data.meshes.new(f"{name}_mesh")
    mesh.from_pydata(verts, [], faces)
    mesh.update()
    obj = bpy.data.objects.new(name, mesh)
    mesh.materials.append(material)
    link_to_collection(obj, collection)
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)
    bpy.ops.object.shade_smooth()
    obj.select_set(False)
    return obj


def component_material(asset: dict, point_kind: str, materials: dict) -> bpy.types.Material:
    color_name = asset["nucleic_color"] if point_kind == "nucleic" else asset["protein_color"]
    return materials[color_name]


def mean_vector(vectors: list[Vector]) -> Vector:
    total = Vector((0.0, 0.0, 0.0))
    for vector in vectors:
        total += vector
    return total / max(1, len(vectors))


def bbox_center(vectors: list[Vector]) -> Vector:
    return Vector(
        (
            (min(v.x for v in vectors) + max(v.x for v in vectors)) * 0.5,
            (min(v.y for v in vectors) + max(v.y for v in vectors)) * 0.5,
            (min(v.z for v in vectors) + max(v.z for v in vectors)) * 0.5,
        )
    )


def principal_axis(vectors: list[Vector]) -> Vector:
    center = mean_vector(vectors)
    cov = [[0.0, 0.0, 0.0] for _ in range(3)]
    for vector in vectors:
        d = vector - center
        values = (d.x, d.y, d.z)
        for i in range(3):
            for j in range(3):
                cov[i][j] += values[i] * values[j]

    axis = Vector((1.0, 0.37, 0.19)).normalized()
    for _ in range(30):
        next_axis = Vector(
            (
                cov[0][0] * axis.x + cov[0][1] * axis.y + cov[0][2] * axis.z,
                cov[1][0] * axis.x + cov[1][1] * axis.y + cov[1][2] * axis.z,
                cov[2][0] * axis.x + cov[2][1] * axis.y + cov[2][2] * axis.z,
            )
        )
        if next_axis.length < 1e-9:
            return Vector((1.0, 0.0, 0.0))
        axis = next_axis.normalized()
    return axis


def build_pdb_component(asset: dict, collections: dict, materials: dict, angstrom_to_mm: float) -> dict:
    pdb_id = asset["pdb_id"].upper()
    path = ASSET_DIR / f"{pdb_id}.cif"
    if not path.exists():
        raise FileNotFoundError(f"Missing RCSB coordinate file: {path}")

    atoms = parse_atom_site(path)
    points = residue_points(atoms)
    sampled = sample_points(points, int(asset["max_beads"]))
    if not sampled:
        raise ValueError(f"No residue points parsed for {pdb_id}")

    all_vectors = [Vector(point["pos_A"]) for point in points]
    min_v = Vector((min(v.x for v in all_vectors), min(v.y for v in all_vectors), min(v.z for v in all_vectors)))
    max_v = Vector((max(v.x for v in all_vectors), max(v.y for v in all_vectors), max(v.z for v in all_vectors)))
    center_A = (min_v + max_v) * 0.5
    bbox_A = max_v - min_v

    nucleic_vectors = [Vector(point["pos_A"]) for point in points if point["kind"] == "nucleic"]
    anchor_kind = asset.get("anchor_kind", "center")
    anchor_vectors = nucleic_vectors if anchor_kind == "nucleic" and nucleic_vectors else all_vectors
    anchor_A = bbox_center(anchor_vectors)

    align_rotation = Matrix.Identity(3)
    axis_before = None
    if asset.get("align_nucleic_axis_to") and len(nucleic_vectors) >= 2:
        source_axis = principal_axis(nucleic_vectors)
        target_axis = Vector(asset["align_nucleic_axis_to"]).normalized()
        if source_axis.dot(target_axis) < 0:
            source_axis = -source_axis
        axis_before = [source_axis.x, source_axis.y, source_axis.z]
        align_rotation = source_axis.rotation_difference(target_axis).to_matrix()

    manual_rotation = Euler(tuple(math.radians(v) for v in asset["rotation_deg"]), "XYZ").to_matrix()
    rotation = manual_rotation @ align_rotation
    location = Vector(asset["location_mm"])
    radius_mm = float(asset["bead_radius_A"]) * angstrom_to_mm
    by_kind = defaultdict(list)
    for point in sampled:
        local = (Vector(point["pos_A"]) - anchor_A) * angstrom_to_mm
        by_kind[point["kind"]].append(location + rotation @ local)

    collection = collections[asset["collection"]]
    for kind, centers in by_kind.items():
        obj = create_bead_mesh(
            f"{asset['name']} ({pdb_id}) {kind}",
            centers,
            radius_mm,
            component_material(asset, kind, materials),
            collection,
        )
        if obj:
            obj["pdb_id"] = pdb_id
            obj["angstrom_to_mm"] = angstrom_to_mm
            obj["source"] = asset["source"]
            obj["source_residue_points"] = len(points)
            obj["sampled_beads"] = len(sampled)
            obj["bead_radius_A"] = asset["bead_radius_A"]
            obj["anchor_kind"] = anchor_kind

    label_offset = asset.get("label_offset_mm", [0.0, 8.0, 0.4])
    create_text(
        f"label_{asset['name']}",
        asset["name"],
        (
            asset["location_mm"][0] + label_offset[0],
            asset["location_mm"][1] + label_offset[1],
            asset["location_mm"][2] + label_offset[2],
        ),
        1.9,
        materials["label_grey"],
        collections["Labels"],
    )
    return {
        "name": asset["name"],
        "pdb_id": pdb_id,
        "atom_count": len(atoms),
        "residue_points": len(points),
        "sampled_beads": len(sampled),
        "bbox_A": [bbox_A.x, bbox_A.y, bbox_A.z],
        "bbox_mm": [bbox_A.x * angstrom_to_mm, bbox_A.y * angstrom_to_mm, bbox_A.z * angstrom_to_mm],
        "location_mm": asset["location_mm"],
        "anchor_kind": anchor_kind,
        "anchor_A": [anchor_A.x, anchor_A.y, anchor_A.z],
        "aligned_nucleic_axis_before": axis_before,
    }


def build_pdb_assets(manifest: dict, collections: dict, materials: dict) -> list[dict]:
    reports = []
    scale = manifest["units"]["angstrom_to_mm"]
    for asset in manifest["pdb_assets"]:
        reports.append(build_pdb_component(asset, collections, materials, scale))
    return reports


def add_lighting_and_camera(collections: dict, materials: dict) -> None:
    light_data = bpy.data.lights.new("large_softbox", "AREA")
    light_data.energy = 3000
    light_data.size = 120
    light = bpy.data.objects.new("large_softbox", light_data)
    light.location = (0.0, -30.0, 180.0)
    bpy.context.scene.collection.objects.link(light)

    camera_data = bpy.data.cameras.new("Camera")
    camera = bpy.data.objects.new("Camera", camera_data)
    camera.location = (0.0, 2.0, 285.0)
    camera.rotation_euler = (0.0, 0.0, 0.0)
    camera_data.type = "ORTHO"
    camera_data.ortho_scale = 175.0
    bpy.context.scene.collection.objects.link(camera)
    bpy.context.scene.camera = camera

def write_report(report: dict, output_path: Path) -> None:
    output_path.write_text(json.dumps(report, indent=2), encoding="utf-8")


def main() -> None:
    manifest = load_manifest()
    clean_scene()
    configure_scene()
    collections = make_collections(manifest["collections"])
    materials = make_materials()

    report = {
        "title": manifest["title"],
        "units": manifest["units"],
        "source_manifest": str(MANIFEST_PATH),
    }
    report["dna"] = build_dna(manifest, collections, materials)
    report["mrna"] = build_mrna(manifest, collections, materials)
    report["scale_bars"] = build_scale_bars(manifest, collections, materials)
    report["pdb_assets"] = build_pdb_assets(manifest, collections, materials)
    add_lighting_and_camera(collections, materials)

    outputs = manifest.get(
        "legacy_outputs",
        {
            "blend": "outputs/canonical/gene_expression_surface_style_v5.blend",
            "preview": "outputs/canonical/preview_gene_expression_surface_style_v5.png",
            "report": "outputs/canonical/gene_expression_surface_scene_v5_report.json",
        },
    )
    blend_path = ROOT / outputs["blend"]
    preview_path = ROOT / outputs["preview"]
    report_path = ROOT / outputs["report"]
    blend_path.parent.mkdir(parents=True, exist_ok=True)
    preview_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    write_report(report, report_path)

    bpy.ops.wm.save_as_mainfile(filepath=str(blend_path))
    bpy.context.scene.render.filepath = str(preview_path)
    bpy.ops.render.render(write_still=True)
    print(f"Wrote {blend_path}")
    print(f"Wrote {preview_path}")
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    raise RuntimeError("scene_core.py is a shared module; run run_canonical_v5_workflow.ps1")
