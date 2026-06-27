#!/usr/bin/env python3
"""Generate optimized V3 procedural DNA/RNA pseudoatom assets for comparison."""

from __future__ import annotations

import json
import math
import random
import shlex
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
import procedural_nucleic_geometry as geom  # noqa: E402


EXPERIMENT_DIR = ROOT / "experiments" / "procedural_nucleic_acids"
ASSET_DIR = EXPERIMENT_DIR / "assets" / "optimization_v3"
CIF_DIR = ASSET_DIR / "cif"
REPORT_PATH = ASSET_DIR / "optimization_v3_assets_report.json"
MANIFEST_PATH = ROOT / "config" / "scene_manifest.json"
RCSB_DIR = ROOT / "assets" / "rcsb"

ANGSTROM_TO_MM = 0.04
DNA_ID = "DNA_OPT_V3"
MRNA_ELONGATED_ID = "MRNA_ELONGATED_OPT_V3"
MRNA_COMPACT_ID = "MRNA_COMPACT_OPT_V3"
DNA_SURFACE_VDW_A = 2.15
RNA_SURFACE_VDW_A = 2.0
DNA_BP_PER_TURN = 10.5
DNA_BASE_ATOMS_PER_PAIR = 11
RNA_BACKBONE_CONNECTORS = 2
RNA_SUGAR_CONNECTORS = 1
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
}
WATER_COMPS = {"HOH", "WAT", "DOD"}


def parse_atom_site(path: Path) -> list[dict]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    atoms = []
    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
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
        if not all(col in index for col in ["_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z"]):
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
    for atom in atoms:
        grouped[(atom["chain"], atom["seq"], atom["comp"])].append(atom)
    points = []
    for (chain, seq, comp), items in grouped.items():
        atom_by_name = {item["atom"]: item for item in items}
        chosen = None
        for atom_name in ("P", "C4'", "C3'", "C1'"):
            if atom_name in atom_by_name:
                chosen = atom_by_name[atom_name]
                break
        if chosen is None:
            x = sum(item["x"] for item in items) / len(items)
            y = sum(item["y"] for item in items) / len(items)
            z = sum(item["z"] for item in items) / len(items)
        else:
            x, y, z = chosen["x"], chosen["y"], chosen["z"]
        points.append({"chain": chain, "seq": seq, "comp": comp, "kind": "nucleic" if comp in NUCLEIC_COMPS else "other", "pos_A": (x, y, z)})
    return points


def dot(a: geom.Vec, b: geom.Vec) -> float:
    return a.x * b.x + a.y * b.y + a.z * b.z


def cross(a: geom.Vec, b: geom.Vec) -> geom.Vec:
    return geom.Vec(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )


def as_vec(point) -> geom.Vec:
    if isinstance(point, geom.Vec):
        return point
    return geom.Vec(*point)


def project_perpendicular(vector: geom.Vec, axis: geom.Vec) -> geom.Vec:
    return vector - axis * dot(vector, axis)


def unwrap_angles(values: list[float]) -> list[float]:
    if not values:
        return values
    unwrapped = [values[0]]
    for value in values[1:]:
        previous = unwrapped[-1]
        while value - previous > math.pi:
            value -= 2.0 * math.pi
        while value - previous < -math.pi:
            value += 2.0 * math.pi
        unwrapped.append(value)
    return unwrapped


def chirality_metric_from_strands(strand_a: list[geom.Vec], strand_b: list[geom.Vec]) -> dict:
    count = min(len(strand_a), len(strand_b))
    if count < 4:
        return {"status": "insufficient_points", "sign": 0, "angle_slope_rad_per_unit": 0.0}
    centers = [(strand_a[i] + strand_b[i]) * 0.5 for i in range(count)]
    axis = (centers[-1] - centers[0]).normalized()
    radial0 = project_perpendicular(strand_a[0] - centers[0], axis).normalized()
    basis_v = cross(axis, radial0).normalized()
    distances = [dot(center - centers[0], axis) for center in centers]
    angles = []
    for a, center in zip(strand_a[:count], centers):
        radial = project_perpendicular(a - center, axis).normalized()
        angles.append(math.atan2(dot(radial, basis_v), dot(radial, radial0)))
    unwrapped = unwrap_angles(angles)
    span = distances[-1] - distances[0]
    slope = 0.0 if abs(span) < 1e-9 else (unwrapped[-1] - unwrapped[0]) / span
    return {
        "status": "ok",
        "sign": 1 if slope >= 0 else -1,
        "angle_slope_rad_per_unit": slope,
        "turns": (unwrapped[-1] - unwrapped[0]) / (2.0 * math.pi),
        "axis_span": span,
    }


def screw_handedness_from_strands(strand_a: list[geom.Vec], strand_b: list[geom.Vec]) -> dict:
    count = min(len(strand_a), len(strand_b))
    if count < 4:
        return {"status": "insufficient_points", "sign": 0, "right_handed": False, "mean_screw": 0.0}
    signed_steps = []
    for index in range(count - 1):
        center0 = (strand_a[index] + strand_b[index]) * 0.5
        center1 = (strand_a[index + 1] + strand_b[index + 1]) * 0.5
        tangent_raw = center1 - center0
        if tangent_raw.length < 1e-9:
            continue
        tangent = tangent_raw.normalized()
        radial0 = project_perpendicular(strand_a[index] - center0, tangent)
        radial1 = project_perpendicular(strand_a[index + 1] - center1, tangent)
        if radial0.length < 1e-9 or radial1.length < 1e-9:
            continue
        radial0 = radial0.normalized()
        radial1 = radial1.normalized()
        signed_steps.append(dot(tangent, cross(radial0, radial1)))
    if not signed_steps:
        return {"status": "degenerate", "sign": 0, "right_handed": False, "mean_screw": 0.0}
    mean_screw = sum(signed_steps) / len(signed_steps)
    sign = 1 if mean_screw >= 0 else -1
    return {
        "status": "ok",
        "sign": sign,
        "right_handed": sign > 0,
        "mean_screw": mean_screw,
        "sample_count": len(signed_steps),
        "definition": "sign(dot(local_axis_tangent, cross(radial_i, radial_i_plus_1))); positive is right-handed",
    }


def paired_dna_strands_from_1bna() -> tuple[list[geom.Vec], list[geom.Vec], list[int]]:
    atoms = parse_atom_site(RCSB_DIR / "1BNA.cif")
    points = [point for point in residue_points(atoms) if point["kind"] == "nucleic"]
    by_chain: dict[str, list[dict]] = {}
    for point in points:
        by_chain.setdefault(point["chain"], []).append(point)

    def seq_key(point: dict) -> int:
        try:
            return int(point["seq"])
        except ValueError:
            return 0

    chains = sorted(by_chain.values(), key=len, reverse=True)[:2]
    if len(chains) < 2:
        return [], [], []
    strand_a = [geom.Vec(*point["pos_A"]) for point in sorted(chains[0], key=seq_key)]
    strand_b = [geom.Vec(*point["pos_A"]) for point in sorted(chains[1], key=seq_key, reverse=True)]
    return strand_a, strand_b, [len(chains[0]), len(chains[1])]


def chirality_metric_from_1bna() -> dict:
    strand_a, strand_b, chain_lengths = paired_dna_strands_from_1bna()
    if len(strand_a) < 2 or len(strand_b) < 2:
        return {"status": "missing_two_dna_chains", "sign": 0}
    result = chirality_metric_from_strands(strand_a, strand_b)
    result["source"] = "1BNA"
    result["chain_lengths"] = chain_lengths
    return result


def screw_handedness_from_1bna() -> dict:
    strand_a, strand_b, chain_lengths = paired_dna_strands_from_1bna()
    if len(strand_a) < 2 or len(strand_b) < 2:
        return {"status": "missing_two_dna_chains", "sign": 0, "right_handed": False}
    result = screw_handedness_from_strands(strand_a, strand_b)
    result["source"] = "1BNA"
    result["chain_lengths"] = chain_lengths
    return result


def polyline_length(points: list[geom.Vec]) -> float:
    return sum((points[i] - points[i - 1]).length for i in range(1, len(points)))


def bounds(points: list[geom.Vec]) -> dict:
    min_x = min(point.x for point in points)
    min_y = min(point.y for point in points)
    min_z = min(point.z for point in points)
    max_x = max(point.x for point in points)
    max_y = max(point.y for point in points)
    max_z = max(point.z for point in points)
    return {
        "min_mm": [min_x, min_y, min_z],
        "max_mm": [max_x, max_y, max_z],
        "bbox_mm": [max_x - min_x, max_y - min_y, max_z - min_z],
    }


def path_frame(tangent: geom.Vec) -> tuple[geom.Vec, geom.Vec]:
    tangent = tangent.normalized()
    normal = geom.Vec(-tangent.y, tangent.x, 0.0)
    if normal.length < 1e-6:
        normal = geom.Vec(0.0, 1.0, 0.0)
    return normal.normalized(), geom.Vec(0.0, 0.0, 1.0)


def mm_to_angstrom(point: geom.Vec) -> tuple[float, float, float]:
    return (point.x / ANGSTROM_TO_MM, point.y / ANGSTROM_TO_MM, point.z / ANGSTROM_TO_MM)


def atom_record(name: str, resname: str, chain: str, resseq: int, point: geom.Vec, element: str, component: str) -> dict:
    return {
        "name": name,
        "resname": resname,
        "chain": chain,
        "resseq": resseq,
        "xyz_A": mm_to_angstrom(point),
        "element": element,
        "component": component,
    }


def add_connector_atoms(
    atoms: list[dict],
    points: list[geom.Vec],
    chain: str,
    resname: str,
    component: str,
    element: str,
    count: int,
) -> None:
    for index in range(len(points) - 1):
        for connector in range(1, count + 1):
            atoms.append(
                atom_record(
                    f"M{connector}",
                    resname,
                    chain,
                    index + 1,
                    points[index].lerp(points[index + 1], connector / (count + 1)),
                    element,
                    component,
                )
            )


def dna_strands_for_theta_sign(
    centerline: geom.SampledPath,
    represented_bp: int,
    bp_to_mm: float,
    strand_center_radius: float,
    theta_sign: int,
) -> tuple[list[geom.Vec], list[geom.Vec]]:
    strand_a: list[geom.Vec] = []
    strand_b: list[geom.Vec] = []
    for i in range(represented_bp + 1):
        distance = min(centerline.length, i * bp_to_mm)
        center = centerline.point_at_length(distance)
        tangent = centerline.tangent_at_length(distance)
        normal, binormal = path_frame(tangent)
        theta = theta_sign * 2.0 * math.pi * i / DNA_BP_PER_TURN
        radial = normal * math.sin(theta) + binormal * math.cos(theta)
        strand_a.append(center + radial * strand_center_radius)
        strand_b.append(center - radial * strand_center_radius)
    return strand_a, strand_b


def build_dna_asset(manifest: dict, calibrator: dict) -> dict:
    units = manifest["units"]
    bp_to_mm = units["dna_bp_to_mm"]
    target_envelope_mm = 0.8
    vdw_mm = DNA_SURFACE_VDW_A * ANGSTROM_TO_MM
    strand_center_radius = max(0.05, target_envelope_mm * 0.5 - vdw_mm)
    centerline = geom.SampledPath(geom.catmull_rom(geom.dna_controls(manifest), 32))
    represented_bp = int(round(centerline.length / bp_to_mm))
    calibrator_screw = screw_handedness_from_1bna()
    candidate_handedness = {}
    for candidate_sign in (-1, 1):
        candidate_a, candidate_b = dna_strands_for_theta_sign(
            centerline, represented_bp, bp_to_mm, strand_center_radius, candidate_sign
        )
        candidate_handedness[str(candidate_sign)] = {
            "chirality": chirality_metric_from_strands(candidate_a, candidate_b),
            "screw_handedness": screw_handedness_from_strands(candidate_a, candidate_b),
        }
    theta_sign = -1
    for candidate_sign, handedness in candidate_handedness.items():
        if handedness["screw_handedness"].get("sign") == calibrator_screw.get("sign", 1):
            theta_sign = int(candidate_sign)
            break
    strand_a: list[geom.Vec] = []
    strand_b: list[geom.Vec] = []
    sugar_a: list[geom.Vec] = []
    sugar_b: list[geom.Vec] = []
    base_rows: list[list[geom.Vec]] = []
    atoms: list[dict] = []
    for i in range(represented_bp + 1):
        distance = min(centerline.length, i * bp_to_mm)
        center = centerline.point_at_length(distance)
        tangent = centerline.tangent_at_length(distance)
        normal, binormal = path_frame(tangent)
        theta = theta_sign * 2.0 * math.pi * i / DNA_BP_PER_TURN
        radial = normal * math.sin(theta) + binormal * math.cos(theta)
        tangent_offset = tangent * (0.018 * math.sin(theta * 0.5))
        a = center + radial * strand_center_radius
        b = center - radial * strand_center_radius
        sa = a - radial * 0.045 + tangent_offset
        sb = b + radial * 0.045 - tangent_offset
        strand_a.append(a)
        strand_b.append(b)
        sugar_a.append(sa)
        sugar_b.append(sb)
        atoms.append(atom_record("P", "DNA", "A", i + 1, a, "P", "strand_A"))
        atoms.append(atom_record("S", "DNA", "A", i + 1, sa, "C", "strand_A"))
        atoms.append(atom_record("P", "DNA", "B", i + 1, b, "P", "strand_B"))
        atoms.append(atom_record("S", "DNA", "B", i + 1, sb, "C", "strand_B"))
        row = []
        for j in range(DNA_BASE_ATOMS_PER_PAIR):
            fraction = (j + 1) / (DNA_BASE_ATOMS_PER_PAIR + 1)
            base_center = a.lerp(b, fraction)
            groove = tangent * (0.016 * math.sin(j * 1.7 + i * 0.33))
            base_point = base_center + groove
            row.append(base_point)
            atoms.append(atom_record(f"B{j + 1}", "DBS", "C", i + 1, base_point, "N" if j % 2 else "C", "base_pairs"))
        base_rows.append(row)
    add_connector_atoms(atoms, strand_a, "A", "DNA", "strand_A", "P", 2)
    add_connector_atoms(atoms, sugar_a, "A", "DNA", "strand_A", "C", 1)
    add_connector_atoms(atoms, strand_b, "B", "DNA", "strand_B", "P", 2)
    add_connector_atoms(atoms, sugar_b, "B", "DNA", "strand_B", "C", 1)
    for j in range(DNA_BASE_ATOMS_PER_PAIR):
        add_connector_atoms(atoms, [row[j] for row in base_rows], "C", "DBS", "base_pairs", "C", 1)
    generated_chirality = chirality_metric_from_strands(strand_a, strand_b)
    generated_screw_handedness = screw_handedness_from_strands(strand_a, strand_b)
    all_points = strand_a + strand_b + [point for row in base_rows for point in row]
    return {
        "asset_id": DNA_ID,
        "atoms": atoms,
        "components": {"strand_A": "chain A", "strand_B": "chain B", "base_pairs": "chain C"},
        "surface_vdw_A": {"strand_A": DNA_SURFACE_VDW_A, "strand_B": DNA_SURFACE_VDW_A, "base_pairs": DNA_SURFACE_VDW_A},
        "geometry": {
            "style": "optimization_v3_right_handed_dna_proxy",
            "handedness": "right_handed",
            "handedness_reference": "1BNA",
            "axis_length_mm": centerline.length,
            "represented_bp": represented_bp,
            "bp_spacing_mm": bp_to_mm,
            "bp_per_turn": DNA_BP_PER_TURN,
            "theta_sign": theta_sign,
            "theta_sign_candidates": candidate_handedness,
            "calibrator_chirality": calibrator,
            "generated_chirality": generated_chirality,
            "chirality_matches_1BNA": generated_chirality.get("sign") == calibrator.get("sign"),
            "calibrator_screw_handedness": calibrator_screw,
            "generated_screw_handedness": generated_screw_handedness,
            "screw_handedness_matches_1BNA": generated_screw_handedness.get("sign") == calibrator_screw.get("sign"),
            "strand_center_radius_mm": strand_center_radius,
            "surface_vdw_A": DNA_SURFACE_VDW_A,
            "estimated_envelope_diameter_mm": 2.0 * (strand_center_radius + vdw_mm),
            "base_atoms_per_pair": DNA_BASE_ATOMS_PER_PAIR,
            "atom_count": len(atoms),
            **bounds(all_points),
        },
    }


def irregular_wave(t: float, spec: dict, rng_values: list[tuple[float, float, float]]) -> float:
    value = 0.0
    for frequency, phase, weight in rng_values:
        value += weight * math.sin(2.0 * math.pi * frequency * t + phase)
    return spec["amplitude"] * value


def smoothstep01(t: float) -> float:
    return t * t * (3.0 - 2.0 * t)


def solve_irregular_segment(
    start: tuple[float, float, float],
    target_length: float,
    end_y: float,
    amplitude: float,
    points_count: int,
    seed: int,
) -> list[geom.Vec]:
    rng = random.Random(seed)
    start_x, start_y, start_z = start
    control_count = max(7, min(30, int(target_length / 6.5)))
    raw_walk = [0.0]
    velocity = 0.0
    for index in range(1, control_count - 1):
        t = index / (control_count - 1)
        envelope = math.sin(math.pi * t) ** 0.65
        velocity = 0.64 * velocity + rng.uniform(-1.0, 1.0) * 0.58
        kink = rng.uniform(-0.36, 0.36)
        raw_walk.append(envelope * (velocity + kink))
    raw_walk.append(0.0)
    max_abs_walk = max(0.05, max(abs(value) for value in raw_walk))
    normalized_walk = [value / max_abs_walk for value in raw_walk]
    z_terms = [(rng.uniform(1.3, 5.6), rng.uniform(0.0, math.tau), rng.uniform(0.25, 0.85)) for _ in range(5)]
    z_norm = max(0.1, sum(abs(term[2]) for term in z_terms))
    z_terms = [(f, p, w / z_norm) for f, p, w in z_terms]

    def points_for_span(span: float) -> list[geom.Vec]:
        control_points = []
        for index, walk_value in enumerate(normalized_walk):
            t = index / (control_count - 1)
            smooth = smoothstep01(t)
            envelope = math.sin(math.pi * t)
            local_notch = 0.22 * amplitude * math.sin(2.0 * math.pi * (control_count * 0.17) * t + seed * 0.11)
            y = start_y + (end_y - start_y) * smooth + envelope * (amplitude * walk_value + local_notch)
            z = start_z + envelope * irregular_wave(t, {"amplitude": amplitude * 0.30}, z_terms)
            control_points.append((start_x + span * t, y, z))
        samples_per_segment = max(4, int(math.ceil(points_count / max(1, control_count - 1))))
        spline = [geom.Vec(*point) for point in geom.catmull_rom(control_points, samples_per_segment)]
        if len(spline) > points_count + 1:
            path = geom.SampledPath(spline)
            return [path.point_at_length(path.length * index / points_count) for index in range(points_count + 1)]
        return spline

    low = 0.0
    high = max(1.0, target_length)
    while polyline_length(points_for_span(high)) < target_length:
        high *= 1.35
    for _ in range(60):
        mid = (low + high) * 0.5
        if polyline_length(points_for_span(mid)) < target_length:
            low = mid
        else:
            high = mid
    points = points_for_span(high)
    points[0] = geom.Vec(start_x, start_y, start_z)
    points[-1] = geom.Vec(points[-1].x, end_y, start_z)
    return points


def compact_mrnp_path(
    center: tuple[float, float, float],
    target_length: float,
    width: float,
    height: float,
    depth: float,
    points_count: int,
    seed: int,
) -> list[geom.Vec]:
    center_vec = geom.Vec(*center)
    rng = random.Random(seed)
    phase_a = rng.uniform(0.0, math.tau)
    phase_b = rng.uniform(0.0, math.tau)
    phase_c = rng.uniform(0.0, math.tau)
    local_radius = min(width, height) * 0.18
    loop_count = max(18, int(target_length / max(3.6, 2.0 * math.pi * local_radius * 0.82)))
    dense_count = max(points_count * 5, loop_count * 36)

    def raw_points(loops: int) -> list[geom.Vec]:
        points = []
        for index in range(dense_count + 1):
            t = index / dense_count
            fold = 2.0 * math.pi * loops * t
            lobe = 2.0 * math.pi * (5.0 * t + 0.06 * math.sin(2.0 * math.pi * 3.0 * t + phase_a))
            drift = 2.0 * math.pi * (1.7 * t + 0.03 * math.sin(2.0 * math.pi * 7.0 * t + phase_b))
            envelope = 0.78 + 0.16 * math.sin(2.0 * math.pi * 11.0 * t + phase_c)
            macro_x = width * (0.22 * math.cos(lobe) + 0.10 * math.cos(drift * 2.3 + phase_a))
            macro_y = height * (0.20 * math.sin(lobe * 0.91 + phase_b) + 0.10 * math.sin(drift * 1.7))
            macro_z = depth * (0.18 * math.sin(lobe * 1.3 + phase_c) + 0.10 * math.cos(drift * 2.1))
            hairpin_x = local_radius * envelope * math.cos(fold)
            hairpin_y = local_radius * (0.76 + 0.10 * math.sin(lobe)) * math.sin(fold)
            hairpin_z = local_radius * 0.55 * math.sin(fold * 0.53 + phase_b)
            kink = 0.18 * local_radius * math.sin(2.0 * math.pi * 29.0 * t + phase_c)
            points.append(
                geom.Vec(
                    center_vec.x + macro_x + hairpin_x,
                    center_vec.y + macro_y + hairpin_y,
                    center_vec.z + macro_z + hairpin_z + kink,
                )
            )
        return points

    points = raw_points(loop_count)
    base_length = polyline_length(points)
    while base_length < target_length * 0.995:
        loop_count += 4
        points = raw_points(loop_count)
        base_length = polyline_length(points)
    scale = target_length / base_length if base_length else 1.0
    return [center_vec + (point - center_vec) * scale for point in points]


def split_path_by_lengths(points: list[geom.Vec], segments: list[dict], nt_to_mm: float) -> list[dict]:
    path = geom.SampledPath(points)
    cursor = 0.0
    segment_models = []
    for segment in segments:
        target = segment["nt"] * nt_to_mm
        count = max(480, int(segment["nt"] * 20))
        segment_points = [
            path.point_at_length(min(path.length, cursor + target * index / count))
            for index in range(count + 1)
        ]
        segment_models.append({"segment": segment, "points": [point.to_tuple() for point in segment_points]})
        cursor += target
    return segment_models


def sample_nt_points(points: list[geom.Vec], nt_count: int, nt_to_mm: float) -> list[tuple[geom.Vec, geom.Vec]]:
    path = geom.SampledPath(points)
    sampled = []
    for i in range(nt_count):
        distance = min(path.length, i * nt_to_mm)
        sampled.append((path.point_at_length(distance), path.tangent_at_length(distance)))
    return sampled


def rna_atoms_for_segments(manifest: dict, segment_models: list[dict], style: str) -> tuple[list[dict], dict]:
    nt_to_mm = manifest["units"]["mrna_nt_to_mm"]
    atoms: list[dict] = []
    chains = ["A", "B", "C"]
    resnames = ["U5R", "CDS", "U3R"]
    component_names = ["utr5", "coding", "utr3"]
    segment_reports = []
    all_path_points: list[geom.Vec] = []
    total_end_to_end = None
    for segment_index, segment_model in enumerate(segment_models):
        segment = segment_model["segment"]
        points = [as_vec(point) for point in segment_model["points"]]
        if all_path_points:
            all_path_points.extend(points[1:])
        else:
            all_path_points.extend(points)
        chain = chains[segment_index]
        resname = resnames[segment_index]
        component = component_names[segment_index]
        sampled = sample_nt_points(points, segment["nt"], nt_to_mm)
        segment_path_length = polyline_length(points)
        segment_end_to_end = (points[-1] - points[0]).length if len(points) > 1 else 0.0
        p_points: list[geom.Vec] = []
        sugar_points: list[geom.Vec] = []
        base_points: list[geom.Vec] = []
        for i, (center, tangent) in enumerate(sampled, start=1):
            normal, binormal = path_frame(tangent)
            theta = 2.0 * math.pi * (i - 1) / 5.5
            radial = normal * math.sin(theta) + binormal * math.cos(theta)
            sugar = center + radial * 0.065
            base_tip = center + radial * 0.16 + tangent * (0.018 * math.sin(theta * 0.7))
            base_mid = center + radial * 0.115
            atoms.append(atom_record("P", resname, chain, i, center, "P", component))
            atoms.append(atom_record("S", resname, chain, i, sugar, "C", component))
            atoms.append(atom_record("B", resname, chain, i, base_tip, "N", component))
            atoms.append(atom_record("BM", resname, chain, i, base_mid, "C", component))
            p_points.append(center)
            sugar_points.append(sugar)
            base_points.append(base_tip)
        add_connector_atoms(atoms, p_points, chain, resname, component, "P", RNA_BACKBONE_CONNECTORS)
        add_connector_atoms(atoms, sugar_points, chain, resname, component, "C", RNA_SUGAR_CONNECTORS)
        segment_reports.append(
            {
                "name": segment["name"],
                "nt": segment["nt"],
                "target_length_mm": segment["nt"] * nt_to_mm,
                "measured_length_mm": segment_path_length,
                "error_mm": segment_path_length - segment["nt"] * nt_to_mm,
                "end_to_end_mm": segment_end_to_end,
                "tortuosity": segment_path_length / segment_end_to_end if segment_end_to_end > 0 else None,
                "bbox_mm": bounds(points)["bbox_mm"],
            }
        )
    if all_path_points:
        total_end_to_end = (all_path_points[-1] - all_path_points[0]).length
    total_measured = sum(item["measured_length_mm"] for item in segment_reports)
    full_bounds = bounds(all_path_points)
    path_generator = "deterministic_correlated_random_control_path"
    if "compact" in style:
        path_generator = "continuous_mRNP_rosette_hairpin_globule"
    return atoms, {
        "style": style,
        "path_generator": path_generator,
        "segments": segment_reports,
        "total_measured_mm": total_measured,
        "total_end_to_end_mm": total_end_to_end,
        "contour_to_max_extent_ratio": total_measured / max(full_bounds["bbox_mm"]) if full_bounds["bbox_mm"] else None,
        "bbox": full_bounds,
        "atom_count": len(atoms),
        "surface_vdw_A": RNA_SURFACE_VDW_A,
    }


def build_mrna_asset(manifest: dict, compact: bool) -> dict:
    mrna = manifest["mrna"]
    nt_to_mm = manifest["units"]["mrna_nt_to_mm"]
    start = tuple(mrna["start_mm"])
    current = geom.Vec(start[0], start[1] + (1.3 if not compact else -18.0), start[2])
    segment_models = []
    if compact:
        total_target = sum(segment["nt"] * nt_to_mm for segment in mrna["segments"])
        full_points = compact_mrnp_path(
            center=(-16.0, -8.0, 0.0),
            target_length=total_target,
            width=10.5,
            height=8.0,
            depth=3.4,
            points_count=max(2200, int(sum(segment["nt"] for segment in mrna["segments"]) * 1.35)),
            seed=304,
        )
        segment_models = split_path_by_lengths(full_points, mrna["segments"], nt_to_mm)
    else:
        specs = [
            {"end_y": 40.6, "amplitude": 1.9, "seed": 201},
            {"end_y": 36.8, "amplitude": 7.4, "seed": 202},
            {"end_y": 31.2, "amplitude": 6.2, "seed": 203},
        ]
        for segment, spec in zip(mrna["segments"], specs):
            target = segment["nt"] * nt_to_mm
            count = max(80, int(segment["nt"] / 2))
            points = solve_irregular_segment(current.to_tuple(), target, spec["end_y"], spec["amplitude"], count, spec["seed"])
            segment_models.append({"segment": segment, "points": [point.to_tuple() for point in points]})
            current = points[-1]
    atoms, geometry = rna_atoms_for_segments(manifest, segment_models, "optimization_v3_compact_mrnp_rosette_hairpin_mrna" if compact else "optimization_v3_irregular_elongated_mrna")
    if compact:
        geometry["model_basis"] = {
            "biological_intent": "schematic compact mRNP-like fold with many local hairpin/rosette turns",
            "not_sequence_specific": True,
            "translation_state": "compact non-translating or storage/stress-associated mRNP, for visual contrast with elongated mRNA",
        }
    return {
        "asset_id": MRNA_COMPACT_ID if compact else MRNA_ELONGATED_ID,
        "atoms": atoms,
        "components": {"utr5": "chain A", "coding": "chain B", "utr3": "chain C"},
        "surface_vdw_A": {"utr5": RNA_SURFACE_VDW_A, "coding": RNA_SURFACE_VDW_A, "utr3": RNA_SURFACE_VDW_A},
        "geometry": geometry,
    }


def cif_value(value: str | int | float) -> str:
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def write_cif(path: Path, data_name: str, atoms: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    headers = [
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_atom_id",
        "_atom_site.pdbx_PDB_model_num",
    ]
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"data_{data_name}\n#\nloop_\n")
        for header in headers:
            handle.write(f"{header}\n")
        for serial, atom in enumerate(atoms, start=1):
            x, y, z = atom["xyz_A"]
            fields = [
                "ATOM",
                serial,
                atom["element"],
                atom["name"],
                ".",
                atom["resname"],
                atom["chain"],
                1,
                atom["resseq"],
                "?",
                x,
                y,
                z,
                1.0,
                0.0,
                atom["resseq"],
                atom["resname"],
                atom["chain"],
                atom["name"],
                1,
            ]
            handle.write(" ".join(cif_value(field) for field in fields) + "\n")
        handle.write("#\n")


def public_asset(asset: dict) -> dict:
    return {
        "asset_id": asset["asset_id"],
        "cif": f"assets/optimization_v3/cif/{asset['asset_id']}.cif",
        "atom_count": len(asset["atoms"]),
        "components": asset["components"],
        "surface_vdw_A": asset["surface_vdw_A"],
        "geometry": asset["geometry"],
    }


def main() -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    CIF_DIR.mkdir(parents=True, exist_ok=True)
    calibrator = chirality_metric_from_1bna()
    assets = [
        build_dna_asset(manifest, calibrator),
        build_mrna_asset(manifest, compact=False),
        build_mrna_asset(manifest, compact=True),
    ]
    for asset in assets:
        write_cif(CIF_DIR / f"{asset['asset_id']}.cif", asset["asset_id"], asset["atoms"])
    elongated = assets[1]["geometry"]
    compact = assets[2]["geometry"]
    report = {
        "title": "DNA/RNA optimization V3 procedural assets",
        "canonical_scene_changed": False,
        "coordinate_units": "angstrom",
        "angstrom_to_mm": ANGSTROM_TO_MM,
        "assets": [public_asset(asset) for asset in assets],
        "comparisons": {
            "compact_vs_elongated_mrna": {
                "elongated_bbox_mm": elongated["bbox"]["bbox_mm"],
                "compact_bbox_mm": compact["bbox"]["bbox_mm"],
                "bbox_x_compaction_ratio": compact["bbox"]["bbox_mm"][0] / elongated["bbox"]["bbox_mm"][0],
                "bbox_area_compaction_ratio": (
                    (compact["bbox"]["bbox_mm"][0] * compact["bbox"]["bbox_mm"][1])
                    / (elongated["bbox"]["bbox_mm"][0] * elongated["bbox"]["bbox_mm"][1])
                ),
            }
        },
    }
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
