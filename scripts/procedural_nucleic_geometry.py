#!/usr/bin/env python3
"""Shared scale-correct procedural DNA/RNA geometry helpers."""

from __future__ import annotations

import math
import random
import shlex
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NUCLEIC_COMPS = {"DA", "DC", "DG", "DT", "DI", "A", "C", "G", "U"}
WATER_COMPS = {"HOH", "WAT", "DOD"}


@dataclass(frozen=True)
class Vec:
    x: float
    y: float
    z: float

    def __iter__(self):
        yield self.x
        yield self.y
        yield self.z

    def __add__(self, other: "Vec") -> "Vec":
        return Vec(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other: "Vec") -> "Vec":
        return Vec(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, scalar: float) -> "Vec":
        return Vec(self.x * scalar, self.y * scalar, self.z * scalar)

    __rmul__ = __mul__

    def __neg__(self) -> "Vec":
        return Vec(-self.x, -self.y, -self.z)

    def dot(self, other: "Vec") -> float:
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other: "Vec") -> "Vec":
        return Vec(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )

    @property
    def length(self) -> float:
        return math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def normalized(self) -> "Vec":
        length = self.length
        if length < 1e-12:
            return Vec(1.0, 0.0, 0.0)
        return self * (1.0 / length)

    def lerp(self, other: "Vec", t: float) -> "Vec":
        return self + (other - self) * t

    def to_tuple(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)


class SampledPath:
    def __init__(self, points: list[tuple[float, float, float] | Vec]):
        self.points = [point if isinstance(point, Vec) else Vec(*point) for point in points]
        self.cumulative = [0.0]
        for i in range(1, len(self.points)):
            self.cumulative.append(self.cumulative[-1] + (self.points[i] - self.points[i - 1]).length)
        self.length = self.cumulative[-1] if self.cumulative else 0.0

    def point_at_length(self, distance: float) -> Vec:
        if not self.points:
            return Vec(0.0, 0.0, 0.0)
        distance = min(max(distance, 0.0), self.length)
        for i in range(1, len(self.cumulative)):
            if self.cumulative[i] >= distance:
                start = self.points[i - 1]
                end = self.points[i]
                span = self.cumulative[i] - self.cumulative[i - 1]
                t = 0.0 if span == 0 else (distance - self.cumulative[i - 1]) / span
                return start.lerp(end, t)
        return self.points[-1]

    def tangent_at_length(self, distance: float) -> Vec:
        if len(self.points) < 2:
            return Vec(1.0, 0.0, 0.0)
        distance = min(max(distance, 0.0), self.length)
        for i in range(1, len(self.cumulative)):
            if self.cumulative[i] >= distance:
                tangent = self.points[i] - self.points[i - 1]
                return tangent.normalized()
        return (self.points[-1] - self.points[-2]).normalized()

    def closest_to_xy(self, location: tuple[float, float, float] | list[float]) -> tuple[Vec, Vec, float]:
        target = Vec(location[0], location[1], 0.0)
        best_i = 0
        best_distance = float("inf")
        for i, point in enumerate(self.points):
            distance = (Vec(point.x, point.y, 0.0) - target).length
            if distance < best_distance:
                best_i = i
                best_distance = distance
        length = self.cumulative[best_i]
        return self.points[best_i], self.tangent_at_length(length), length

    def closest_to_x(self, x: float) -> tuple[Vec, Vec, float]:
        best_i = min(range(len(self.points)), key=lambda i: abs(self.points[i].x - x))
        length = self.cumulative[best_i]
        return self.points[best_i], self.tangent_at_length(length), length


def vec_tuple(point: Vec) -> tuple[float, float, float]:
    return (point.x, point.y, point.z)


def polyline_length(points: list[tuple[float, float, float] | Vec]) -> float:
    vectors = [point if isinstance(point, Vec) else Vec(*point) for point in points]
    return sum((vectors[i] - vectors[i - 1]).length for i in range(1, len(vectors)))


def bounds(points: list[tuple[float, float, float] | Vec]) -> dict:
    vectors = [point if isinstance(point, Vec) else Vec(*point) for point in points]
    min_x = min(point.x for point in vectors)
    min_y = min(point.y for point in vectors)
    min_z = min(point.z for point in vectors)
    max_x = max(point.x for point in vectors)
    max_y = max(point.y for point in vectors)
    max_z = max(point.z for point in vectors)
    return {
        "min_mm": [min_x, min_y, min_z],
        "max_mm": [max_x, max_y, max_z],
        "bbox_mm": [max_x - min_x, max_y - min_y, max_z - min_z],
    }


def catmull_rom(control_points: list[tuple[float, float, float]], samples_per_segment: int = 24) -> list[tuple[float, float, float]]:
    points = [Vec(*point) for point in control_points]
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
            value = (
                (p1 * 2.0)
                + (p2 - p0) * t
                + (p0 * 2.0 - p1 * 5.0 + p2 * 4.0 - p3) * t2
                + (-p0 + p1 * 3.0 - p2 * 3.0 + p3) * t3
            ) * 0.5
            result.append(value.to_tuple())
    result.append(points[-1].to_tuple())
    return result


def path_frame(tangent: Vec) -> tuple[Vec, Vec]:
    tangent = tangent.normalized()
    up = Vec(0.0, 0.0, 1.0)
    normal = up.cross(tangent)
    if normal.length < 1e-6:
        up = Vec(0.0, 1.0, 0.0)
        normal = up.cross(tangent)
    normal = normal.normalized()
    binormal = tangent.cross(normal).normalized()
    return normal, binormal


def mean_vec(points: list[Vec]) -> Vec:
    total = Vec(0.0, 0.0, 0.0)
    for point in points:
        total += point
    return total * (1.0 / max(1, len(points)))


def bbox_center(points: list[Vec]) -> Vec:
    return Vec(
        (min(point.x for point in points) + max(point.x for point in points)) * 0.5,
        (min(point.y for point in points) + max(point.y for point in points)) * 0.5,
        (min(point.z for point in points) + max(point.z for point in points)) * 0.5,
    )


def covariance_axes(points: list[Vec]) -> list[tuple[float, Vec]]:
    center = mean_vec(points)
    cov = [[0.0, 0.0, 0.0] for _ in range(3)]
    for point in points:
        delta = point - center
        values = (delta.x, delta.y, delta.z)
        for i in range(3):
            for j in range(3):
                cov[i][j] += values[i] * values[j]
    if points:
        inv_n = 1.0 / len(points)
        for i in range(3):
            for j in range(3):
                cov[i][j] *= inv_n

    eigenvectors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    for _ in range(40):
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
        axis = Vec(eigenvectors[0][i], eigenvectors[1][i], eigenvectors[2][i])
        axes.append((cov[i][i], axis.normalized()))
    return sorted(axes, key=lambda item: item[0], reverse=True)


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
            while i < len(lines) and lines[i].strip() not in {"#", "loop_"}:
                i += 1
            continue
        index = {col: n for n, col in enumerate(cols)}
        if "_atom_site.Cartn_x" not in index:
            continue

        def value(parts: list[str], *names: str) -> str:
            for name in names:
                n = index.get(name)
                if n is not None and n < len(parts):
                    return parts[n].strip("'\"")
            return ""

        while i < len(lines):
            raw = lines[i].strip()
            if raw == "#" or raw == "loop_" or raw.startswith("_"):
                break
            if not raw:
                i += 1
                continue
            parts = raw.split()
            if len(parts) != len(cols):
                parts = shlex.split(raw)
            if len(parts) >= len(cols):
                comp = value(parts, "_atom_site.auth_comp_id", "_atom_site.label_comp_id").upper()
                model = value(parts, "_atom_site.pdbx_PDB_model_num")
                alt = value(parts, "_atom_site.label_alt_id")
                if comp not in WATER_COMPS and (not model or model in {"1", ".", "?"}) and alt in {"", ".", "?", "A", "1"}:
                    try:
                        atoms.append(
                            {
                                "element": value(parts, "_atom_site.type_symbol").upper(),
                                "atom": value(parts, "_atom_site.auth_atom_id", "_atom_site.label_atom_id"),
                                "comp": comp,
                                "chain": value(parts, "_atom_site.auth_asym_id", "_atom_site.label_asym_id"),
                                "seq": value(parts, "_atom_site.auth_seq_id", "_atom_site.label_seq_id"),
                                "x": float(value(parts, "_atom_site.Cartn_x")),
                                "y": float(value(parts, "_atom_site.Cartn_y")),
                                "z": float(value(parts, "_atom_site.Cartn_z")),
                            }
                        )
                    except ValueError:
                        pass
            i += 1
    return atoms


def representative_residue_points(atoms: list[dict], chains: set[str]) -> dict[str, list[tuple[int, Vec]]]:
    grouped: dict[tuple[str, int, str], list[dict]] = {}
    for atom in atoms:
        if atom["comp"] not in NUCLEIC_COMPS or atom["chain"] not in chains:
            continue
        try:
            seq = int(atom["seq"])
        except ValueError:
            continue
        grouped.setdefault((atom["chain"], seq, atom["comp"]), []).append(atom)

    by_chain: dict[str, list[tuple[int, Vec]]] = {chain: [] for chain in chains}
    for (chain, seq, _comp), items in grouped.items():
        atom_by_name = {item["atom"]: item for item in items}
        chosen = None
        for atom_name in ("C1'", "C4'", "P", "C3'"):
            if atom_name in atom_by_name:
                chosen = atom_by_name[atom_name]
                break
        if chosen is None:
            chosen = {
                "x": sum(item["x"] for item in items) / len(items),
                "y": sum(item["y"] for item in items) / len(items),
                "z": sum(item["z"] for item in items) / len(items),
            }
        by_chain[chain].append((seq, Vec(chosen["x"], chosen["y"], chosen["z"])))
    for chain in by_chain:
        by_chain[chain].sort(key=lambda item: item[0])
    return by_chain


def unwrapped_turns(points: list[Vec], center: Vec, x_axis: Vec, y_axis: Vec) -> float:
    if len(points) < 2:
        return 0.0
    angles = [math.atan2((point - center).dot(y_axis), (point - center).dot(x_axis)) for point in points]
    unwrapped = [angles[0]]
    for angle in angles[1:]:
        while angle - unwrapped[-1] > math.pi:
            angle -= math.tau
        while angle - unwrapped[-1] < -math.pi:
            angle += math.tau
        unwrapped.append(angle)
    return (unwrapped[-1] - unwrapped[0]) / math.tau


def centerline_superhelix_metrics(points: list[Vec], z_axis: Vec | None = None, x_axis: Vec | None = None) -> dict:
    center = mean_vec(points)
    axes = covariance_axes(points)
    source_z = (z_axis or axes[-1][1]).normalized()
    endpoint = points[-1] - points[0]
    projected_endpoint = endpoint - source_z * endpoint.dot(source_z)
    source_x = (x_axis or (projected_endpoint if projected_endpoint.length > 5.0 else axes[0][1])).normalized()
    source_y = source_z.cross(source_x).normalized()
    source_x = source_y.cross(source_z).normalized()
    turns = unwrapped_turns(points, center, source_x, source_y)
    radii = []
    heights = []
    for point in points:
        local = point - center
        radii.append(math.sqrt(local.dot(source_x) ** 2 + local.dot(source_y) ** 2))
        heights.append(local.dot(source_z))
    return {
        "center_A": center,
        "axis": source_z,
        "x_axis": source_x,
        "y_axis": source_y,
        "superhelical_turns": turns,
        "superhelical_handedness": "left" if turns < 0 else "right",
        "mean_superhelical_radius_A": sum(radii) / len(radii),
        "min_superhelical_radius_A": min(radii),
        "max_superhelical_radius_A": max(radii),
        "height_span_A": max(heights) - min(heights),
        "path_length_A": polyline_length(points),
        "endpoint_distance_A": endpoint.length,
        "eigenvalues": [value for value, _axis in axes],
    }


def nucleosome_bp_centerline_from_pdb(loop: dict) -> dict:
    pdb_id = loop.get("source_pdb_id", "1AOI").upper()
    pdb_path = ROOT / "assets" / "rcsb" / f"{pdb_id}.cif"
    if not pdb_path.exists():
        raise FileNotFoundError(f"Missing nucleosome guide structure: {pdb_path}")
    chain_a = loop.get("chain_a", "I")
    chain_b = loop.get("chain_b", "J")
    by_chain = representative_residue_points(parse_atom_site(pdb_path), {chain_a, chain_b})
    strand_a = [point for _seq, point in by_chain.get(chain_a, [])]
    strand_b = [point for _seq, point in by_chain.get(chain_b, [])]
    count = min(len(strand_a), len(strand_b))
    if count < 24:
        raise ValueError(f"Too few paired nucleosome DNA residues in {pdb_id}: {count}")
    requested_pairing = loop.get("pairing", "same_order")
    candidates = {}
    for pairing, paired_b in (
        ("same_order", strand_b),
        ("reverse_order", list(reversed(strand_b))),
    ):
        distances = [(strand_a[index] - paired_b[index]).length for index in range(count)]
        centers = [(strand_a[index] + paired_b[index]) * 0.5 for index in range(count)]
        metrics = centerline_superhelix_metrics(centers)
        candidates[pairing] = {
            "points": centers,
            "mean_pair_distance_A": sum(distances) / len(distances),
            "min_pair_distance_A": min(distances),
            "max_pair_distance_A": max(distances),
            "metrics": metrics,
        }
    if requested_pairing in {"auto", "auto_antiparallel_min_distance"}:
        pairing = min(candidates, key=lambda key: candidates[key]["mean_pair_distance_A"])
    else:
        pairing = "reverse_order" if requested_pairing == "reverse_order" else "same_order"
    selected = candidates[pairing]
    return {
        "pdb_id": pdb_id,
        "chain_a": chain_a,
        "chain_b": chain_b,
        "pairing": pairing,
        "requested_pairing": requested_pairing,
        "candidate_pair_distances_A": {
            key: {
                "mean": value["mean_pair_distance_A"],
                "min": value["min_pair_distance_A"],
                "max": value["max_pair_distance_A"],
            }
            for key, value in candidates.items()
        },
        **selected,
    }


def basis_transform_points(
    source_points_A: list[Vec],
    source_center_A: Vec,
    source_basis: tuple[Vec, Vec, Vec],
    target_center_mm: Vec,
    target_basis: tuple[Vec, Vec, Vec],
    angstrom_to_mm: float,
) -> list[Vec]:
    source_x, source_y, source_z = source_basis
    target_x, target_y, target_z = target_basis
    transformed = []
    for point_A in source_points_A:
        local = (point_A - source_center_A) * angstrom_to_mm
        transformed.append(
            target_center_mm
            + target_x * local.dot(source_x)
            + target_y * local.dot(source_y)
            + target_z * local.dot(source_z)
        )
    return transformed


def nucleosome_loop_guide(manifest: dict, loop_override: dict | None = None) -> dict | None:
    proxy = manifest.get("procedural_nucleic_acids", {}).get(
        "dna",
        manifest.get("procedural_nucleic_surfaces", {}).get("dna", {}),
    )
    loop = dict(proxy.get("nucleosome_loop", {}))
    if loop_override:
        loop.update(loop_override)
    if not loop.get("enabled", False) or not loop.get("source_pdb_id"):
        return None
    guide_source = nucleosome_bp_centerline_from_pdb(loop)
    source_points_A = guide_source["points"]
    source_center_A = mean_vec(source_points_A)
    axes = covariance_axes(source_points_A)
    source_z = axes[-1][1]
    endpoint_axis = (source_points_A[-1] - source_points_A[0])
    endpoint_axis_projected = endpoint_axis - source_z * endpoint_axis.dot(source_z)
    source_x = endpoint_axis_projected.normalized() if endpoint_axis_projected.length > 5.0 else axes[0][1]
    source_y = source_z.cross(source_x).normalized()
    source_x = source_y.cross(source_z).normalized()
    turns = unwrapped_turns(source_points_A, source_center_A, source_x, source_y)
    if loop.get("expected_superhelical_handedness", "left") == "left" and turns > 0:
        source_z = -source_z
        source_y = source_z.cross(source_x).normalized()
        source_x = source_y.cross(source_z).normalized()
        turns = unwrapped_turns(source_points_A, source_center_A, source_x, source_y)
    source_metrics = centerline_superhelix_metrics(source_points_A, source_z, source_x)

    roll = math.radians(float(loop.get("roll_deg", 0.0)))
    target_z = Vec(0.0, 0.0, 1.0)
    target_x = Vec(math.cos(roll), math.sin(roll), 0.0)
    target_y = target_z.cross(target_x).normalized()
    target_center_mm = Vec(*loop.get("center_mm", (72.0, -45.0, 0.0)))
    transformed = basis_transform_points(
        source_points_A,
        source_center_A,
        (source_x, source_y, source_z),
        target_center_mm,
        (target_x, target_y, target_z),
        manifest["units"]["angstrom_to_mm"],
    )
    if endpoint_axis.length > 5.0 and transformed[0].x > transformed[-1].x:
        source_points_A = list(reversed(source_points_A))
        transformed = list(reversed(transformed))
        source_x = -source_x
        source_y = source_z.cross(source_x).normalized()
    guide_path = SampledPath(transformed)
    return {
        "source_pdb_id": guide_source["pdb_id"],
        "source_bp_count": len(source_points_A),
        "source_points_A": source_points_A,
        "scene_points_mm": transformed,
        "source_center_A": source_center_A,
        "source_basis": (source_x, source_y, source_z),
        "target_center_mm": target_center_mm,
        "target_basis": (target_x, target_y, target_z),
        "scene_length_mm": guide_path.length,
        "scene_bounds": bounds(transformed),
        "pairing": guide_source["pairing"],
        "requested_pairing": guide_source["requested_pairing"],
        "candidate_pair_distances_A": guide_source["candidate_pair_distances_A"],
        "mean_pair_distance_A": guide_source["mean_pair_distance_A"],
        "min_pair_distance_A": guide_source["min_pair_distance_A"],
        "max_pair_distance_A": guide_source["max_pair_distance_A"],
        "superhelical_turns": source_metrics["superhelical_turns"],
        "superhelical_handedness": source_metrics["superhelical_handedness"],
        "mean_superhelical_radius_A": source_metrics["mean_superhelical_radius_A"],
        "mean_superhelical_radius_mm": source_metrics["mean_superhelical_radius_A"] * manifest["units"]["angstrom_to_mm"],
        "height_span_A": source_metrics["height_span_A"],
        "height_span_mm": source_metrics["height_span_A"] * manifest["units"]["angstrom_to_mm"],
        "source_path_length_A": source_metrics["path_length_A"],
        "terminal_distance_A": endpoint_axis.length,
        "roll_deg": float(loop.get("roll_deg", 0.0)),
    }


def path_nucleosome_loop_override(manifest: dict, config: dict, arc_path: SampledPath) -> dict:
    loop = manifest.get("procedural_nucleic_acids", {}).get("dna", {}).get("nucleosome_loop", {})
    fraction = float(config.get("nucleosome_loop_fraction", 0.77))
    distance = max(0.0, min(arc_path.length, arc_path.length * fraction))
    center = arc_path.point_at_length(distance)
    tangent = arc_path.tangent_at_length(distance)
    roll_deg = math.degrees(math.atan2(tangent.y, tangent.x))
    override = {
        "enabled": True,
        "center_mm": center.to_tuple(),
        "roll_deg": roll_deg + float(config.get("nucleosome_loop_roll_offset_deg", 0.0)),
        "computed_from_arrangement_path": True,
        "arrangement_path_fraction": fraction,
        "arrangement_path_distance_mm": distance,
        "entry_bridge_mm": float(loop.get("entry_bridge_mm", config.get("nucleosome_entry_bridge_mm", 8.0))),
        "exit_bridge_mm": float(loop.get("exit_bridge_mm", config.get("nucleosome_exit_bridge_mm", 8.0))),
    }
    return override


def rounded_serpentine_controls(config: dict, xy_scale: float = 1.0) -> list[tuple[float, float, float]]:
    center = Vec(*config.get("center_mm", (0.0, 0.0, 0.0)))
    width = float(config.get("width_mm", 112.0)) * xy_scale
    height = float(config.get("height_mm", 74.0)) * xy_scale
    lanes = max(2, int(config.get("lanes", 6)))
    samples_per_lane = max(4, int(config.get("samples_per_lane", 16)))
    irregularity = float(config.get("irregularity_mm", 3.2)) * xy_scale
    z_amp = float(config.get("z_irregularity_mm", 0.42))
    phase = float(config.get("phase", 0.53))
    length_variation = float(config.get("length_variation", 0.23))
    controls: list[tuple[float, float, float]] = []
    lane_pitch = height / (lanes - 1)
    turn_samples = max(10, samples_per_lane)
    lane_widths = []
    for lane in range(lanes):
        lane_phase = math.sin(math.tau * (0.37 * lane + phase)) + 0.45 * math.sin(math.tau * (0.19 * lane + phase * 1.7))
        width_scale = 1.0 + length_variation * lane_phase / 1.45
        lane_widths.append(width * max(0.62, min(1.22, width_scale)))

    def lane_y(lane_index: int, local_t: float) -> float:
        lane_t = lane_index / (lanes - 1)
        y_base = -height * 0.5 + height * lane_t
        global_t = (lane_index + local_t) / lanes
        return y_base + irregularity * (
            0.42 * math.sin(math.tau * (1.3 * global_t + phase))
            + 0.18 * math.sin(math.tau * (4.4 * global_t + 0.18))
        )

    for lane in range(lanes):
        lane_width = lane_widths[lane]
        x_start = -lane_width * 0.5 if lane % 2 == 0 else lane_width * 0.5
        x_end = lane_width * 0.5 if lane % 2 == 0 else -lane_width * 0.5
        for step in range(samples_per_lane):
            local_t = step / max(1, samples_per_lane - 1)
            global_t = (lane + local_t) / lanes
            envelope = math.sin(math.pi * min(1.0, max(0.0, global_t))) ** 0.35
            x = x_start + (x_end - x_start) * local_t
            y = lane_y(lane, local_t)
            z = z_amp * envelope * math.sin(math.tau * (5.2 * global_t + phase))
            controls.append((center.x + x, center.y + y, center.z + z))
        if lane < lanes - 1:
            next_width = lane_widths[lane + 1]
            side = 1.0 if lane % 2 == 0 else -1.0
            p0 = Vec(*controls[-1]) - center
            x1 = side * next_width * 0.5
            y1 = lane_y(lane + 1, 0.0 if (lane + 1) % 2 == 1 else 1.0)
            p3 = Vec(x1, y1, 0.0)
            handle = max(8.0 * xy_scale, lane_pitch * float(config.get("turn_handle_fraction", 0.62)))
            p1 = p0 + Vec(side * handle, 0.0, 0.0)
            p2 = p3 + Vec(side * handle, 0.0, 0.0)
            for step in range(1, turn_samples):
                u = step / turn_samples
                omt = 1.0 - u
                point = p0 * (omt ** 3) + p1 * (3.0 * omt * omt * u) + p2 * (3.0 * omt * u * u) + p3 * (u ** 3)
                global_t = (lane + u) / lanes
                z = z_amp * 0.35 * math.sin(math.tau * (3.1 * global_t + phase))
                controls.append((center.x + point.x, center.y + point.y, center.z + z))
    return controls


def reader_order_serpentine_controls(config: dict, xy_scale: float = 1.0) -> list[tuple[float, float, float]]:
    center = Vec(*config.get("center_mm", (0.0, 0.0, 0.0)))
    width = float(config.get("width_mm", 112.0)) * xy_scale
    height = float(config.get("height_mm", 74.0)) * xy_scale
    lanes = max(2, int(config.get("lanes", 5)))
    samples_per_lane = max(4, int(config.get("samples_per_lane", 16)))
    irregularity = float(config.get("irregularity_mm", 3.2)) * xy_scale
    z_amp = float(config.get("z_irregularity_mm", 0.42))
    phase = float(config.get("phase", 0.53))
    length_variation = float(config.get("length_variation", 0.23))
    controls: list[tuple[float, float, float]] = []
    lane_pitch = height / (lanes - 1)
    turn_samples = max(10, samples_per_lane)
    lane_widths = []
    for lane in range(lanes):
        lane_phase = math.sin(math.tau * (0.37 * lane + phase)) + 0.45 * math.sin(math.tau * (0.19 * lane + phase * 1.7))
        width_scale = 1.0 + length_variation * lane_phase / 1.45
        lane_widths.append(width * max(0.62, min(1.22, width_scale)))

    def lane_y(lane_index: int, local_t: float) -> float:
        lane_t = lane_index / (lanes - 1)
        y_base = height * 0.5 - height * lane_t
        global_t = (lane_index + local_t) / lanes
        return y_base + irregularity * (
            0.42 * math.sin(math.tau * (1.3 * global_t + phase))
            + 0.18 * math.sin(math.tau * (4.4 * global_t + 0.18))
        )

    for lane in range(lanes):
        lane_width = lane_widths[lane]
        x_start = -lane_width * 0.5 if lane % 2 == 0 else lane_width * 0.5
        x_end = lane_width * 0.5 if lane % 2 == 0 else -lane_width * 0.5
        for step in range(samples_per_lane):
            local_t = step / max(1, samples_per_lane - 1)
            global_t = (lane + local_t) / lanes
            envelope = math.sin(math.pi * min(1.0, max(0.0, global_t))) ** 0.35
            x = x_start + (x_end - x_start) * local_t
            y = lane_y(lane, local_t)
            z = z_amp * envelope * math.sin(math.tau * (5.2 * global_t + phase))
            controls.append((center.x + x, center.y + y, center.z + z))
        if lane < lanes - 1:
            next_width = lane_widths[lane + 1]
            side = 1.0 if lane % 2 == 0 else -1.0
            p0 = Vec(*controls[-1]) - center
            x1 = side * next_width * 0.5
            y1 = lane_y(lane + 1, 0.0 if (lane + 1) % 2 == 1 else 1.0)
            p3 = Vec(x1, y1, 0.0)
            handle = max(8.0 * xy_scale, lane_pitch * float(config.get("turn_handle_fraction", 0.62)))
            p1 = p0 + Vec(side * handle, 0.0, 0.0)
            p2 = p3 + Vec(side * handle, 0.0, 0.0)
            for step in range(1, turn_samples):
                u = step / turn_samples
                omt = 1.0 - u
                point = p0 * (omt ** 3) + p1 * (3.0 * omt * omt * u) + p2 * (3.0 * omt * u * u) + p3 * (u ** 3)
                global_t = (lane + u) / lanes
                z = z_amp * 0.35 * math.sin(math.tau * (3.1 * global_t + phase))
                controls.append((center.x + point.x, center.y + point.y, center.z + z))
    return controls


def full_gene_serpentine_nucleosome_layout(
    manifest: dict,
    config: dict,
    controls_builder=rounded_serpentine_controls,
) -> dict:
    bp_to_mm = manifest["units"]["dna_bp_to_mm"]
    target_bp = int(config.get("target_bp", 3954))
    target_length = target_bp * bp_to_mm
    entry_bridge = float(config.get("nucleosome_entry_bridge_mm", 8.0))
    exit_bridge = float(config.get("nucleosome_exit_bridge_mm", 8.0))

    def build_for_scale(scale: float) -> dict:
        base_controls = controls_builder(config, scale)
        base_path = SampledPath(catmull_rom(base_controls, 32))
        loop_override = path_nucleosome_loop_override(manifest, config, base_path)
        loop_override["entry_bridge_mm"] = entry_bridge
        loop_override["exit_bridge_mm"] = exit_bridge
        guide = nucleosome_loop_guide(manifest, loop_override)
        if guide is None:
            return {"controls": base_controls, "path": SampledPath(catmull_rom(base_controls, 32)), "loop_override": loop_override, "guide": None}
        insert_distance = loop_override["arrangement_path_distance_mm"]
        before_end = max(0.0, insert_distance - entry_bridge)
        after_start = min(base_path.length, insert_distance + exit_bridge)
        before_count = max(8, int(before_end / 5.0))
        after_count = max(8, int((base_path.length - after_start) / 5.0))
        guide_points = guide["scene_points_mm"]
        start = guide_points[0]
        end = guide_points[-1]
        start_tangent = (guide_points[1] - start).normalized()
        end_tangent = (end - guide_points[-2]).normalized()
        controls: list[Vec] = []
        for index in range(before_count + 1):
            controls.append(base_path.point_at_length(before_end * index / before_count))
        controls.append(start - start_tangent * min(3.5, entry_bridge * 0.45))
        controls.extend(guide_points)
        controls.append(end + end_tangent * min(3.5, exit_bridge * 0.45))
        for index in range(1, after_count + 1):
            distance = after_start + (base_path.length - after_start) * index / after_count
            controls.append(base_path.point_at_length(distance))
        path = SampledPath(catmull_rom([point.to_tuple() for point in controls], 32))
        return {"controls": [point.to_tuple() for point in controls], "path": path, "loop_override": loop_override, "guide": guide}

    low = float(config.get("min_xy_scale", 0.45))
    high = float(config.get("max_xy_scale", 2.2))
    while build_for_scale(high)["path"].length < target_length:
        high *= 1.25
    best = build_for_scale(high)
    for _ in range(40):
        mid = (low + high) * 0.5
        trial = build_for_scale(mid)
        if trial["path"].length < target_length:
            low = mid
        else:
            high = mid
            best = trial
    best["target_bp"] = target_bp
    best["target_length_mm"] = target_length
    best["xy_scale"] = high
    return best


def full_gene_serpentine_with_nucleosome_loop_controls(manifest: dict, config: dict) -> list[tuple[float, float, float]]:
    return full_gene_serpentine_nucleosome_layout(manifest, config)["controls"]


def v5_reader_order_serpentine_nucleosome_layout(manifest: dict, config: dict) -> dict:
    return full_gene_serpentine_nucleosome_layout(manifest, config, reader_order_serpentine_controls)


def v5_reader_order_serpentine_with_nucleosome_loop_controls(manifest: dict, config: dict) -> list[tuple[float, float, float]]:
    return v5_reader_order_serpentine_nucleosome_layout(manifest, config)["controls"]


def organic_dna_arc_controls(manifest: dict, config: dict) -> list[tuple[float, float, float]]:
    center = Vec(*config.get("center_mm", (0.0, -8.0, 0.0)))
    radius = float(config.get("radius_mm", 64.0))
    start_angle = math.radians(float(config.get("start_angle_deg", 215.0)))
    span = math.radians(float(config.get("span_deg", -240.0)))
    samples = int(config.get("samples", 32))
    radial_amp = float(config.get("radial_irregularity_mm", 4.2))
    tangential_amp = float(config.get("tangential_irregularity_mm", 1.4))
    z_amp = float(config.get("z_irregularity_mm", 0.75))
    controls: list[tuple[float, float, float]] = []
    for index in range(samples + 1):
        t = index / samples
        angle = start_angle + span * t
        envelope = math.sin(math.pi * t) ** 0.55
        radial_noise = envelope * radial_amp * (
            0.58 * math.sin(math.tau * (1.35 * t + 0.08))
            + 0.31 * math.sin(math.tau * (3.10 * t + 0.31))
            + 0.18 * math.sin(math.tau * (6.70 * t + 0.17))
        )
        tangent_noise = envelope * tangential_amp * (
            0.70 * math.sin(math.tau * (2.20 * t + 0.46))
            + 0.25 * math.sin(math.tau * (5.30 * t + 0.12))
        )
        r = radius + radial_noise
        radial = Vec(math.cos(angle), math.sin(angle), 0.0)
        tangent = Vec(-math.sin(angle), math.cos(angle), 0.0)
        z = center.z + envelope * z_amp * (
            0.65 * math.sin(math.tau * (2.65 * t + 0.21))
            + 0.35 * math.sin(math.tau * (5.15 * t + 0.48))
        )
        point = center + radial * r + tangent * tangent_noise + Vec(0.0, 0.0, z - center.z)
        controls.append(point.to_tuple())
    return controls


def dna_controls(manifest: dict) -> list[tuple[float, float, float]]:
    dna = manifest["dna"]
    proxy = manifest.get("procedural_nucleic_acids", {}).get(
        "dna",
        manifest.get("procedural_nucleic_surfaces", {}).get("dna", {}),
    )
    if proxy.get("path_mode") == "custom_dna_controls":
        return [tuple(point) for point in proxy.get("custom_controls_mm", [])]
    if proxy.get("path_mode") == "full_gene_serpentine_with_nucleosome_loop":
        return full_gene_serpentine_with_nucleosome_loop_controls(manifest, proxy.get("full_gene_serpentine", {}))
    if proxy.get("path_mode") in {
        "v5_reader_order_serpentine_with_nucleosome_loop",
        "v6_reader_order_serpentine_with_nucleosome_loop",
    }:
        return v5_reader_order_serpentine_with_nucleosome_loop_controls(manifest, proxy.get("full_gene_serpentine", {}))
    loop = proxy.get("nucleosome_loop", {})
    if loop.get("enabled", False):
        guide = nucleosome_loop_guide(manifest)
        if guide is not None:
            loop_points_vec = guide["scene_points_mm"]
            start = loop_points_vec[0]
            end = loop_points_vec[-1]
            start_tangent = (loop_points_vec[1] - start).normalized()
            end_tangent = (end - loop_points_vec[-2]).normalized()
            entry = start - start_tangent * float(loop.get("entry_bridge_mm", 6.0))
            exit = end + end_tangent * float(loop.get("exit_bridge_mm", 6.0))
            return [
                (dna["start_x_mm"], -49.0, dna["center_z_mm"]),
                (-98.0, -46.8, 0.0),
                (-62.0, -43.8, 0.0),
                (-35.0, -45.2, 0.0),
                (-20.0, -44.3, 0.0),
                (-6.0, -45.0, 0.0),
                (22.0, -43.3, 0.0),
                (52.0, -44.5, 0.0),
                entry.to_tuple(),
                *[point.to_tuple() for point in loop_points_vec],
                exit.to_tuple(),
                (96.0, -44.5, 0.0),
                (dna["end_x_mm"], -42.8, 0.0),
            ]
        center = Vec(*loop.get("center_mm", (72.0, -45.0, 0.0)))
        radius = float(loop.get("radius_mm", 1.93))
        turns = float(loop.get("turns", 1.65))
        samples = int(loop.get("samples", 72))
        direction = -1.0 if float(loop.get("direction", -1.0)) < 0 else 1.0
        start_angle = math.radians(float(loop.get("start_angle_deg", 205.0)))
        arc = direction * turns * math.tau
        z_rise = float(loop.get("z_rise_mm", 1.35))
        loop_points = []
        for index in range(samples + 1):
            t = index / samples
            angle = start_angle + arc * t
            loop_points.append(
                (
                    center.x + radius * math.cos(angle),
                    center.y + radius * math.sin(angle),
                    center.z + z_rise * (t - 0.5),
                )
            )
        return [
            (dna["start_x_mm"], -49.0, dna["center_z_mm"]),
            (-98.0, -46.8, 0.0),
            (-62.0, -43.8, 0.0),
            (-35.0, -45.2, 0.0),
            (-20.0, -44.3, 0.0),
            (-6.0, -45.0, 0.0),
            (22.0, -43.3, 0.0),
            (52.0, -44.5, 0.0),
            (64.0, -45.0, -0.15),
            *loop_points,
            (82.0, -45.5, 0.15),
            (96.0, -44.5, 0.0),
            (dna["end_x_mm"], -42.8, 0.0),
        ]
    return [
        (dna["start_x_mm"], -49.0, dna["center_z_mm"]),
        (-98.0, -46.8, 0.0),
        (-62.0, -43.8, 0.0),
        (-35.0, -45.2, 0.0),
        (-20.0, -44.3, 0.0),
        (-6.0, -45.0, 0.0),
        (22.0, -43.3, 0.0),
        (72.0, -46.0, 0.0),
        (dna["end_x_mm"], -42.8, 0.0),
    ]


def dna_nucleosome_loop_report(manifest: dict) -> dict | None:
    proxy = manifest.get("procedural_nucleic_acids", {}).get(
        "dna",
        manifest.get("procedural_nucleic_surfaces", {}).get("dna", {}),
    )
    loop = proxy.get("nucleosome_loop", {})
    if not loop.get("enabled", False):
        return None
    units = manifest["units"]
    guide_override = None
    if proxy.get("path_mode") == "full_gene_serpentine_with_nucleosome_loop":
        layout = full_gene_serpentine_nucleosome_layout(manifest, proxy.get("full_gene_serpentine", {}))
        guide_override = layout["loop_override"]
    if proxy.get("path_mode") in {
        "v5_reader_order_serpentine_with_nucleosome_loop",
        "v6_reader_order_serpentine_with_nucleosome_loop",
    }:
        layout = v5_reader_order_serpentine_nucleosome_layout(manifest, proxy.get("full_gene_serpentine", {}))
        guide_override = layout["loop_override"]
    guide = nucleosome_loop_guide(manifest, guide_override)
    if guide is not None:
        source_center = guide["source_center_A"]
        source_basis = guide["source_basis"]
        target_center = guide["target_center_mm"]
        target_basis = guide["target_basis"]
        return {
            "enabled": True,
            "basis": "pdb_guided_base_pair_centerline",
            "source_pdb_id": guide["source_pdb_id"],
            "source_bp_count": guide["source_bp_count"],
            "center_mm": target_center.to_tuple(),
            "roll_deg": guide["roll_deg"],
            "pairing": guide["pairing"],
            "requested_pairing": guide["requested_pairing"],
            "candidate_pair_distances_A": guide["candidate_pair_distances_A"],
            "mean_pair_distance_A": guide["mean_pair_distance_A"],
            "min_pair_distance_A": guide["min_pair_distance_A"],
            "max_pair_distance_A": guide["max_pair_distance_A"],
            "superhelical_handedness": guide["superhelical_handedness"],
            "superhelical_turns": guide["superhelical_turns"],
            "mean_superhelical_radius_A": guide["mean_superhelical_radius_A"],
            "mean_superhelical_radius_mm": guide["mean_superhelical_radius_mm"],
            "height_span_A": guide["height_span_A"],
            "height_span_mm": guide["height_span_mm"],
            "source_path_length_A": guide["source_path_length_A"],
            "terminal_distance_A": guide["terminal_distance_A"],
            "approximate_centerline_length_mm": guide["scene_length_mm"],
            "approximate_bp_equivalent": guide["scene_length_mm"] / units["dna_bp_to_mm"],
            "scene_bounds": guide["scene_bounds"],
            "source_center_A": source_center.to_tuple(),
            "source_basis": [axis.to_tuple() for axis in source_basis],
            "target_basis": [axis.to_tuple() for axis in target_basis],
            "computed_from_arrangement_path": bool(guide_override),
            "arrangement_path_fraction": guide_override.get("arrangement_path_fraction") if guide_override else None,
            "arrangement_path_distance_mm": guide_override.get("arrangement_path_distance_mm") if guide_override else None,
            "entry_bridge_mm": guide_override.get("entry_bridge_mm") if guide_override else loop.get("entry_bridge_mm"),
            "exit_bridge_mm": guide_override.get("exit_bridge_mm") if guide_override else loop.get("exit_bridge_mm"),
        }
    radius = float(loop.get("radius_mm", 1.93))
    turns = float(loop.get("turns", 1.65))
    z_rise = float(loop.get("z_rise_mm", 1.35))
    planar_arc = radius * turns * math.tau
    approximate_3d_arc = math.sqrt(planar_arc * planar_arc + z_rise * z_rise)
    return {
        "enabled": True,
        "center_mm": loop.get("center_mm", [72.0, -45.0, 0.0]),
        "radius_mm": radius,
        "turns": turns,
        "direction": loop.get("direction", -1),
        "start_angle_deg": loop.get("start_angle_deg", 205.0),
        "z_rise_mm": z_rise,
        "approximate_centerline_length_mm": approximate_3d_arc,
        "approximate_bp_equivalent": approximate_3d_arc / units["dna_bp_to_mm"],
    }


def build_dna_model(manifest: dict) -> dict:
    units = manifest["units"]
    proxy = manifest.get("procedural_nucleic_acids", {}).get(
        "dna",
        manifest.get("procedural_nucleic_surfaces", {}).get("dna", {}),
    )
    bp_to_mm = units["dna_bp_to_mm"]
    pitch_bp = proxy.get("bp_per_turn", 10.5)
    target_envelope_diameter_mm = proxy.get("envelope_diameter_nm", 2.0) * units["nm_to_mm"]
    strand_radius_mm = proxy.get("strand_radius_mm")
    surface_vdw_A = proxy.get("surface_vdw_A", 2.4)
    surface_vdw_mm = (strand_radius_mm if strand_radius_mm is not None else surface_vdw_A * units["angstrom_to_mm"])
    strand_center_radius_mm = proxy.get(
        "strand_center_radius_mm",
        max(0.05, target_envelope_diameter_mm * 0.5 - surface_vdw_mm),
    )
    base_atoms_per_pair = int(proxy.get("base_atoms_per_pair", 11))
    theta_sign = int(proxy.get("theta_sign", -1))
    centerline = SampledPath(catmull_rom(dna_controls(manifest), 32))
    represented_bp = max(1, int(round(centerline.length / bp_to_mm)))
    strand_a: list[Vec] = []
    strand_b: list[Vec] = []
    sugar_a: list[Vec] = []
    sugar_b: list[Vec] = []
    base_pair_atoms: list[Vec] = []
    base_pair_rows: list[list[Vec]] = []
    for i in range(represented_bp + 1):
        distance = min(centerline.length, i * bp_to_mm)
        center = centerline.point_at_length(distance)
        tangent = centerline.tangent_at_length(distance)
        normal, binormal = path_frame(tangent)
        theta = theta_sign * 2.0 * math.pi * i / pitch_bp
        radial = normal * math.sin(theta) + binormal * math.cos(theta)
        tangent_offset = tangent * (0.018 * math.sin(theta * 0.5))
        a = center + radial * strand_center_radius_mm
        b = center - radial * strand_center_radius_mm
        sa = a - radial * 0.045 + tangent_offset
        sb = b + radial * 0.045 - tangent_offset
        strand_a.append(a)
        strand_b.append(b)
        sugar_a.append(sa)
        sugar_b.append(sb)
        row = []
        for j in range(base_atoms_per_pair):
            fraction = (j + 1) / (base_atoms_per_pair + 1)
            groove = tangent * (0.016 * math.sin(j * 1.7 + i * 0.33))
            point = a.lerp(b, fraction) + groove
            row.append(point)
            base_pair_atoms.append(point)
        base_pair_rows.append(row)
    return {
        "path": centerline,
        "strand_a": strand_a,
        "strand_b": strand_b,
        "sugar_a": sugar_a,
        "sugar_b": sugar_b,
        "base_pair_atoms": base_pair_atoms,
        "base_pair_rows": base_pair_rows,
        "report": {
            "style": proxy.get("style", "procedural_blender_surface_proxy"),
            "path_mode": proxy.get("path_mode", "default_scene_dna"),
            "handedness": "right_handed",
            "theta_sign": theta_sign,
            "axis_length_mm": centerline.length,
            "represented_bp": represented_bp,
            "bp_spacing_mm": bp_to_mm,
            "pitch_bp": pitch_bp,
            "strand_center_radius_mm": strand_center_radius_mm,
            "strand_radius_mm": surface_vdw_mm,
            "surface_vdw_A": surface_vdw_A if strand_radius_mm is None else None,
            "target_envelope_diameter_mm": target_envelope_diameter_mm,
            "estimated_envelope_diameter_mm": 2.0 * (strand_center_radius_mm + surface_vdw_mm),
            "strand_a_length_mm": polyline_length(strand_a),
            "base_atoms_per_pair": base_atoms_per_pair,
            "nucleosome_loop": dna_nucleosome_loop_report(manifest),
        },
    }


def smoothstep(t: float) -> float:
    return t * t * (3.0 - 2.0 * t)


def irregular_wave(t: float, amplitude: float, terms: list[tuple[float, float, float]]) -> float:
    return amplitude * sum(weight * math.sin(2.0 * math.pi * frequency * t + phase) for frequency, phase, weight in terms)


def solve_irregular_segment(
    start: tuple[float, float, float],
    target_length: float,
    end_y: float,
    amplitude: float,
    points_count: int,
    seed: int,
) -> list[Vec]:
    rng = random.Random(seed)
    start_x, start_y, start_z = start
    control_count = max(7, min(30, int(target_length / 6.5)))
    raw_walk = [0.0]
    velocity = 0.0
    for index in range(1, control_count - 1):
        t = index / (control_count - 1)
        envelope = math.sin(math.pi * t) ** 0.65
        velocity = 0.64 * velocity + rng.uniform(-1.0, 1.0) * 0.58
        raw_walk.append(envelope * (velocity + rng.uniform(-0.36, 0.36)))
    raw_walk.append(0.0)
    max_abs_walk = max(0.05, max(abs(value) for value in raw_walk))
    normalized_walk = [value / max_abs_walk for value in raw_walk]
    z_terms = [(rng.uniform(1.3, 5.6), rng.uniform(0.0, math.tau), rng.uniform(0.25, 0.85)) for _ in range(5)]
    z_norm = max(0.1, sum(abs(term[2]) for term in z_terms))
    z_terms = [(frequency, phase, weight / z_norm) for frequency, phase, weight in z_terms]

    def points_for_span(span: float) -> list[Vec]:
        control_points = []
        for index, walk_value in enumerate(normalized_walk):
            t = index / (control_count - 1)
            envelope = math.sin(math.pi * t)
            local_notch = 0.22 * amplitude * math.sin(2.0 * math.pi * (control_count * 0.17) * t + seed * 0.11)
            y = start_y + (end_y - start_y) * smoothstep(t) + envelope * (amplitude * walk_value + local_notch)
            z = start_z + envelope * irregular_wave(t, amplitude * 0.30, z_terms)
            control_points.append((start_x + span * t, y, z))
        samples_per_segment = max(4, int(math.ceil(points_count / max(1, control_count - 1))))
        spline = [Vec(*point) for point in catmull_rom(control_points, samples_per_segment)]
        if len(spline) > points_count + 1:
            path = SampledPath(spline)
            return [path.point_at_length(path.length * index / points_count) for index in range(points_count + 1)]
        return spline

    low = 0.0
    high = max(target_length, 1.0)
    while polyline_length(points_for_span(high)) < target_length:
        high *= 1.35
    for _ in range(60):
        mid = (low + high) * 0.5
        if polyline_length(points_for_span(mid)) < target_length:
            low = mid
        else:
            high = mid
    points = points_for_span(high)
    points[0] = Vec(start_x, start_y, start_z)
    points[-1] = Vec(points[-1].x, end_y, start_z)
    return points


def compact_mrnp_path(
    center: tuple[float, float, float],
    target_length: float,
    width: float,
    height: float,
    depth: float,
    points_count: int,
    seed: int,
) -> list[Vec]:
    center_vec = Vec(*center)
    rng = random.Random(seed)
    phase_a = rng.uniform(0.0, math.tau)
    phase_b = rng.uniform(0.0, math.tau)
    phase_c = rng.uniform(0.0, math.tau)
    local_radius = min(width, height) * 0.18
    loop_count = max(18, int(target_length / max(3.6, 2.0 * math.pi * local_radius * 0.82)))
    dense_count = max(points_count * 5, loop_count * 36)

    def raw_points(loops: int) -> list[Vec]:
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
            points.append(Vec(center_vec.x + macro_x + hairpin_x, center_vec.y + macro_y + hairpin_y, center_vec.z + macro_z + hairpin_z + kink))
        return points

    points = raw_points(loop_count)
    base_length = polyline_length(points)
    while base_length < target_length * 0.995:
        loop_count += 4
        points = raw_points(loop_count)
        base_length = polyline_length(points)
    scale = target_length / base_length if base_length else 1.0
    return [center_vec + (point - center_vec) * scale for point in points]


def _bounded_integer_allocation(
    total: int,
    count: int,
    minimum: int,
    maximum: int,
    rng: random.Random,
) -> list[int]:
    """Distribute an integer total deterministically within inclusive bounds."""
    if count <= 0:
        return []
    total = max(count * minimum, min(count * maximum, total))
    values = [minimum] * count
    remaining = total - count * minimum
    order = list(range(count))
    rng.shuffle(order)
    while remaining:
        changed = False
        for index in order:
            if remaining <= 0:
                break
            if values[index] < maximum:
                values[index] += 1
                remaining -= 1
                changed = True
        if not changed:
            break
        order = order[1:] + order[:1]
    return values


def _sampled_self_clearance(points: list[Vec], sample_count: int = 320) -> float | None:
    """Estimate non-local backbone clearance without an expensive all-pairs scan."""
    if len(points) < 8:
        return None
    stride = max(1, len(points) // sample_count)
    sampled = [(index, points[index]) for index in range(0, len(points), stride)]
    exclusion = stride * 3
    best = float("inf")
    for left in range(len(sampled)):
        index_a, point_a = sampled[left]
        for right in range(left + 1, len(sampled)):
            index_b, point_b = sampled[right]
            if abs(index_b - index_a) <= exclusion:
                continue
            best = min(best, (point_b - point_a).length)
    return best if math.isfinite(best) else None


def structured_stem_loop_path(
    guide_points: list[Vec],
    total_nt: int,
    nt_to_mm: float,
    settings: dict,
    *,
    compact: bool,
) -> dict:
    """Fold a guide into a continuous, schematic secondary-structure-rich RNA.

    Each motif follows the primary chain up one stem arm, around a hairpin loop,
    and down its antiparallel partner.  The final uniform scale is deliberately
    small and enforces the canonical contour-length invariant exactly.
    """
    secondary = settings.get("secondary_structure", {})
    stem_count = int(
        secondary.get("compact_stem_count" if compact else "elongated_stem_count", 40 if compact else 24)
    )
    paired_target = float(
        secondary.get(
            "compact_paired_fraction_target" if compact else "elongated_paired_fraction_target",
            0.60 if compact else 0.40,
        )
    )
    seed = int(secondary.get("compact_seed" if compact else "elongated_seed", 51852 if compact else 41852))
    rng = random.Random(seed)
    stem_min = int(secondary.get("stem_bp_min", 6))
    stem_max = int(secondary.get("stem_bp_max", 18))
    loop_min = int(secondary.get("hairpin_loop_nt_min", 4))
    loop_max = int(secondary.get("hairpin_loop_nt_max", 12))
    requested_bp = round(total_nt * paired_target * 0.5)
    stem_bp = _bounded_integer_allocation(requested_bp, stem_count, stem_min, stem_max, rng)
    loop_nt = [rng.randint(loop_min, loop_max) for _ in range(stem_count)]
    allocated = 2 * sum(stem_bp) + sum(loop_nt)
    if allocated >= total_nt:
        raise ValueError(f"Secondary-structure allocation consumes {allocated}/{total_nt} nt")
    linker_nt = total_nt - allocated

    guide = SampledPath(guide_points)
    if guide.length <= 0:
        raise ValueError("Structured RNA guide has zero length")
    anchor_distances = [guide.length * index / (stem_count + 1) for index in range(stem_count + 2)]
    anchors = [guide.point_at_length(distance) for distance in anchor_distances]
    tangents = [guide.tangent_at_length(distance) for distance in anchor_distances]

    points: list[Vec] = [anchors[0]]
    pair_indices: list[tuple[int, int]] = []
    loop_ranges: list[tuple[int, int]] = []
    connector_samples = max(3, int(linker_nt / max(stem_count + 1, 1) / 3))
    a_form_diameter = float(secondary.get("a_form_diameter_mm", 0.92))

    for stem_index in range(stem_count):
        anchor = anchors[stem_index + 1]
        connector_start = points[-1]
        for step in range(1, connector_samples + 1):
            points.append(connector_start.lerp(anchor, step / connector_samples))

        tangent = tangents[stem_index + 1].normalized()
        if compact:
            angle = math.tau * (stem_index * 0.38196601125 + rng.uniform(-0.06, 0.06))
            axis = Vec(math.cos(angle), math.sin(angle), 0.16 * math.sin(angle * 1.7)).normalized()
        else:
            normal, binormal = path_frame(tangent)
            angle = (stem_index % 6 - 2.5) * 0.30 + rng.uniform(-0.10, 0.10)
            axis = (normal * math.cos(angle) + binormal * math.sin(angle)).normalized()
        separation_axis = tangent.cross(axis).normalized()
        if separation_axis.length < 1e-6:
            separation_axis = path_frame(axis)[0]

        loop_length = loop_nt[stem_index] * nt_to_mm
        separation = min(a_form_diameter * 0.62, max(0.28, 2.0 * loop_length / math.pi))
        leg_length = stem_bp[stem_index] * nt_to_mm
        base_a = anchor - separation_axis * (separation * 0.5)
        base_b = anchor + separation_axis * (separation * 0.5)
        top_center = anchor + axis * leg_length
        top_a = top_center - separation_axis * (separation * 0.5)
        top_b = top_center + separation_axis * (separation * 0.5)

        points.append(base_a)
        outbound: list[int] = []
        arm_samples = max(5, stem_bp[stem_index] * 2)
        for step in range(arm_samples + 1):
            t = step / arm_samples
            subtle_twist = 0.035 * math.sin(math.tau * stem_bp[stem_index] * t / 11.0)
            point = base_a.lerp(top_a, t) + tangent * subtle_twist
            points.append(point)
            outbound.append(len(points) - 1)

        loop_start = len(points)
        loop_samples = max(8, loop_nt[stem_index] * 3)
        radius = separation * 0.5
        for step in range(1, loop_samples + 1):
            theta = math.pi - math.pi * step / loop_samples
            points.append(top_center + separation_axis * (radius * math.cos(theta)) + axis * (radius * math.sin(theta)))
        loop_ranges.append((loop_start, len(points) - 1))

        inbound: list[int] = []
        for step in range(1, arm_samples + 1):
            t = step / arm_samples
            subtle_twist = -0.035 * math.sin(math.tau * stem_bp[stem_index] * (1.0 - t) / 11.0)
            point = top_b.lerp(base_b, t) + tangent * subtle_twist
            points.append(point)
            inbound.append(len(points) - 1)
        for pair_step in range(stem_bp[stem_index]):
            out_index = outbound[min(len(outbound) - 1, round(pair_step * (len(outbound) - 1) / max(stem_bp[stem_index] - 1, 1)))]
            in_index = inbound[min(len(inbound) - 1, round((stem_bp[stem_index] - 1 - pair_step) * (len(inbound) - 1) / max(stem_bp[stem_index] - 1, 1)))]
            pair_indices.append((out_index, in_index))

    final_start = points[-1]
    for step in range(1, connector_samples + 1):
        points.append(final_start.lerp(anchors[-1], step / connector_samples))

    raw_length = polyline_length(points)
    target_length = total_nt * nt_to_mm
    scale_origin = mean_vec(points) if compact else points[0]
    scale = target_length / raw_length if raw_length else 1.0
    scaled_points = [scale_origin + (point - scale_origin) * scale for point in points]
    base_pairs = [(scaled_points[left].to_tuple(), scaled_points[right].to_tuple()) for left, right in pair_indices]
    measured = polyline_length(scaled_points)
    structure_bounds = bounds(scaled_points)
    paired_nt = 2 * sum(stem_bp)
    report = {
        "model": "deterministic_schematic_stem_loop",
        "sequence_resolved": False,
        "compact": compact,
        "seed": seed,
        "total_nt": total_nt,
        "target_length_mm": target_length,
        "measured_length_mm": measured,
        "length_error_mm": measured - target_length,
        "stem_count": stem_count,
        "hairpin_loop_count": stem_count,
        "junction_count": max(3, stem_count // 8),
        "stem_bp": stem_bp,
        "paired_nt": paired_nt,
        "paired_fraction": paired_nt / total_nt,
        "loop_nt": sum(loop_nt),
        "linker_nt": linker_nt,
        "allocation_check_nt": paired_nt + sum(loop_nt) + linker_nt,
        "base_pair_bridge_count": len(base_pairs),
        "a_form_bp_per_turn": float(secondary.get("a_form_bp_per_turn", 11.0)),
        "a_form_diameter_mm": a_form_diameter,
        "uniform_length_scale": scale,
        "bbox": structure_bounds,
        "sampled_self_clearance_mm": _sampled_self_clearance(scaled_points),
        "biological_intent": (
            "secondary-structure-rich compact RNA schematic" if compact else "moderately structured elongated translating RNA schematic"
        ),
    }
    return {"points": scaled_points, "base_pairs": base_pairs, "loop_ranges": loop_ranges, "report": report}


def split_path_by_lengths(points: list[Vec], segments: list[dict], nt_to_mm: float) -> list[dict]:
    path = SampledPath(points)
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


def conical_spiral_mrna_points(manifest: dict, config: dict) -> tuple[list[Vec], dict]:
    target_length = sum(segment["nt"] for segment in manifest["mrna"]["segments"]) * manifest["units"]["mrna_nt_to_mm"]
    start = Vec(*config.get("start_mm", manifest["mrna"].get("start_mm", (-112.0, 39.0, 0.0))))
    start_radius = float(config.get("start_radius_mm", 28.0))
    end_radius = float(config.get("end_radius_mm", 7.0))
    vertical_rise = float(config.get("vertical_rise_mm", 95.0))
    rise_axis = str(config.get("rise_axis", "y")).lower()
    theta0 = math.radians(float(config.get("start_angle_deg", 0.0)))
    orientation = -1.0 if float(config.get("orientation", -1.0)) < 0 else 1.0
    samples = int(config.get("samples", 1800))
    radial_noise = float(config.get("radial_irregularity_mm", 2.1))
    vertical_noise = float(config.get("vertical_irregularity_mm", 2.4))
    lateral_noise = float(config.get("lateral_irregularity_mm", 3.0))
    z_amp = float(config.get("z_irregularity_mm", 1.2))
    terminal_drift = Vec(*config.get("terminal_drift_mm", (0.0, 0.0, 0.0)))
    center0 = start - Vec(math.cos(theta0), math.sin(theta0), 0.0) * start_radius

    def raw_points(turns: float) -> list[Vec]:
        points: list[Vec] = []
        for index in range(samples + 1):
            t = index / samples
            ease = smoothstep(t)
            radius = start_radius * (1.0 - ease) + end_radius * ease
            envelope = math.sin(math.pi * t) ** 0.70
            angle = theta0 + orientation * math.tau * turns * t
            local_radius = radius + envelope * radial_noise * (
                0.55 * math.sin(math.tau * (4.1 * t + 0.11))
                + 0.30 * math.sin(math.tau * (9.7 * t + 0.39))
                + 0.15 * math.sin(math.tau * (17.0 * t + 0.23))
            )
            axial_wobble = vertical_noise * envelope * math.sin(math.tau * (2.25 * t + 0.33))
            lateral_x = lateral_noise * math.sin(math.tau * (0.70 * t + 0.15)) * envelope
            lateral_y = lateral_noise * 0.72 * math.sin(math.tau * (1.15 * t + 0.42)) * envelope
            z_wobble = z_amp * envelope * (
                0.60 * math.sin(math.tau * (3.5 * t + 0.20))
                + 0.40 * math.sin(math.tau * (8.1 * t + 0.42))
            )
            drift = terminal_drift * ease
            if rise_axis == "z":
                center = center0 + drift + Vec(lateral_x, lateral_y, vertical_rise * ease + axial_wobble)
                points.append(center + Vec(math.cos(angle), math.sin(angle), 0.0) * local_radius + Vec(0.0, 0.0, z_wobble))
            else:
                center = center0 + drift + Vec(lateral_x, vertical_rise * ease + axial_wobble, 0.0)
                z = start.z + z_wobble
                points.append(center + Vec(math.cos(angle), math.sin(angle), 0.0) * local_radius + Vec(0.0, 0.0, z - center.z))
        points[0] = start
        return points

    low = float(config.get("min_turns", 0.5))
    high = float(config.get("max_turns", 5.0))
    while polyline_length(raw_points(high)) < target_length:
        high *= 1.25
    for _ in range(60):
        mid = (low + high) * 0.5
        if polyline_length(raw_points(mid)) < target_length:
            low = mid
        else:
            high = mid
    points = raw_points(high)
    report = {
        "target_length_mm": target_length,
        "turns": high,
        "start_radius_mm": start_radius,
        "end_radius_mm": end_radius,
        "vertical_rise_mm": vertical_rise,
        "start_mm": start.to_tuple(),
        "end_mm": points[-1].to_tuple(),
        "center0_mm": center0.to_tuple(),
        "orientation": orientation,
        "rise_axis": rise_axis,
        "radial_irregularity_mm": radial_noise,
        "vertical_irregularity_mm": vertical_noise,
        "lateral_irregularity_mm": lateral_noise,
        "z_irregularity_mm": z_amp,
        "terminal_drift_mm": terminal_drift.to_tuple(),
    }
    return points, report


def mrna_report_from_segments(manifest: dict, segment_models: list[dict], style: str, path_generator: str) -> dict:
    proxy = manifest.get("procedural_nucleic_acids", {}).get(
        "mrna",
        manifest.get("procedural_nucleic_surfaces", {}).get("mrna", {}),
    )
    segment_reports = []
    all_points: list[Vec] = []
    for segment_model in segment_models:
        segment = segment_model["segment"]
        points = [point if isinstance(point, Vec) else Vec(*point) for point in segment_model["points"]]
        if all_points:
            all_points.extend(points[1:])
        else:
            all_points.extend(points)
        measured = polyline_length(points)
        end_to_end = (points[-1] - points[0]).length if len(points) > 1 else 0.0
        segment_reports.append(
            {
                "name": segment["name"],
                "nt": segment["nt"],
                "target_length_mm": segment["nt"] * manifest["units"]["mrna_nt_to_mm"],
                "measured_length_mm": measured,
                "error_mm": measured - segment["nt"] * manifest["units"]["mrna_nt_to_mm"],
                "end_to_end_mm": end_to_end,
                "tortuosity": measured / end_to_end if end_to_end > 0 else None,
                "bbox_mm": bounds(points)["bbox_mm"],
                "start_mm": list(points[0]),
                "end_mm": list(points[-1]),
            }
        )
    full_bounds = bounds(all_points)
    total_measured = sum(item["measured_length_mm"] for item in segment_reports)
    total_end_to_end = (all_points[-1] - all_points[0]).length if all_points else None
    return {
        "style": style,
        "path_generator": path_generator,
        "segments": segment_reports,
        "total_measured_mm": total_measured,
        "total_end_to_end_mm": total_end_to_end,
        "contour_to_max_extent_ratio": total_measured / max(full_bounds["bbox_mm"]) if full_bounds["bbox_mm"] else None,
        "bbox": full_bounds,
        "marker_interval_nt": 50,
        "surface_vdw_A": proxy.get("surface_vdw_A", 2.0),
        "base_offset_mm": proxy.get("base_offset_mm", 0.16),
    }


def build_mrna_model(manifest: dict) -> dict:
    units = manifest["units"]
    mrna = manifest["mrna"]
    nt_to_mm = units["mrna_nt_to_mm"]
    segments = mrna["segments"]
    proxy = manifest.get("procedural_nucleic_acids", {}).get(
        "mrna",
        manifest.get("procedural_nucleic_surfaces", {}).get("mrna", {}),
    )
    # The polymerase is attached to path point zero.  A nascent transcript's
    # growing end at polymerase is its 3' end, so allocate the mature-mRNA
    # color blocks away from Pol II in 3' UTR -> CDS -> 5' UTR order.  Keep the
    # manifest itself in conventional 5' -> 3' biological order.
    path_segments = (
        list(reversed(segments))
        if proxy.get("elongated_path_origin") == "3_prime_at_polymerase_ii"
        else segments
    )
    if proxy.get("path_mode") == "custom_mrna_path":
        raw_points = [Vec(*point) for point in proxy.get("custom_points_mm", [])]
        if len(raw_points) < 2:
            raise ValueError("custom_mrna_path requires at least two custom_points_mm")
        target_length = sum(segment["nt"] for segment in segments) * nt_to_mm
        source_path = SampledPath(raw_points)
        if source_path.length <= 0:
            raise ValueError("custom_mrna_path has zero length")
        start = raw_points[0]
        scale = target_length / source_path.length
        scaled_points = [start + (point - start) * scale for point in raw_points]
        segment_models = split_path_by_lengths(scaled_points, path_segments, nt_to_mm)
        marker_centers = []
        scaled_path = SampledPath(scaled_points)
        for nt in range(50, sum(segment["nt"] for segment in segments) + 1, 50):
            marker_centers.append(scaled_path.point_at_length(nt * nt_to_mm))
        report = mrna_report_from_segments(
            manifest,
            segment_models,
            "custom_curve_rebaked_mrna",
            "edited_source_curve_scaled_to_exact_contour_length",
        )
        report["custom_path"] = {
            "source_length_mm": source_path.length,
            "target_length_mm": target_length,
            "uniform_scale_about_start": scale,
            "source_points": len(raw_points),
        }
        return {
            "path": SampledPath([point.to_tuple() for point in scaled_points]),
            "segments": segment_models,
            "marker_centers": marker_centers,
            "report": report,
        }
    if proxy.get("path_mode") == "canonical_polymerase_spiral":
        spiral_config = proxy.get("canonical_spiral", {})
        points, spiral_report = conical_spiral_mrna_points(manifest, spiral_config)
        segment_models = split_path_by_lengths(points, path_segments, nt_to_mm)
        marker_centers = []
        for nt in range(50, sum(segment["nt"] for segment in segments) + 1, 50):
            marker_centers.append(SampledPath(points).point_at_length(nt * nt_to_mm))
        report = mrna_report_from_segments(
            manifest,
            segment_models,
            "canonical_v5_polymerase_origin_conical_spiral_mrna",
            "deterministic_conical_spiral_tightening_to_top",
        )
        report["spiral"] = spiral_report
        report["path_origin"] = proxy.get("elongated_path_origin", "5_prime")
        report["path_segment_order_from_origin"] = [segment["name"] for segment in path_segments]
        return {
            "path": SampledPath([point.to_tuple() for point in points]),
            "segments": segment_models,
            "marker_centers": marker_centers,
            "report": report,
        }
    start = tuple(mrna["start_mm"])
    current = (start[0], start[1] + 1.3, start[2])
    segment_specs = [
        {"end_y": 40.6, "amplitude": 1.9, "seed": 201},
        {"end_y": 36.8, "amplitude": 7.4, "seed": 202},
        {"end_y": 31.2, "amplitude": 6.2, "seed": 203},
    ]
    all_points: list[tuple[float, float, float]] = []
    segment_models = []
    marker_centers = []
    nt_cursor = 0
    for segment, spec in zip(segments, segment_specs):
        target = segment["nt"] * nt_to_mm
        count = max(80, int(segment["nt"] / 2))
        points = solve_irregular_segment(current, target, spec["end_y"], spec["amplitude"], count, spec["seed"])
        if all_points:
            all_points.extend([point.to_tuple() for point in points[1:]])
        else:
            all_points.extend([point.to_tuple() for point in points])
        for nt in range(nt_cursor + 50, nt_cursor + segment["nt"] + 1, 50):
            fraction = (nt - nt_cursor) / segment["nt"]
            marker_centers.append(points[min(len(points) - 1, int(fraction * (len(points) - 1)))])
        segment_models.append({"segment": segment, "points": [point.to_tuple() for point in points]})
        nt_cursor += segment["nt"]
        current = points[-1].to_tuple()
    report = mrna_report_from_segments(manifest, segment_models, "optimization_v3_irregular_elongated_mrna", "deterministic_correlated_random_control_path")
    return {
        "path": SampledPath(all_points),
        "segments": segment_models,
        "marker_centers": marker_centers,
        "report": report,
    }


def build_compact_mrna_model(manifest: dict) -> dict:
    mrna = manifest["mrna"]
    nt_to_mm = manifest["units"]["mrna_nt_to_mm"]
    proxy = manifest.get("procedural_nucleic_acids", {}).get("mrna", {})
    total_nt = sum(segment["nt"] for segment in mrna["segments"])
    total_target = total_nt * nt_to_mm
    center = tuple(proxy.get("compact_center_mm", [66.0, 43.0, 0.0]))
    full_points = compact_mrnp_path(
        center=center,
        target_length=total_target,
        width=float(proxy.get("compact_width_mm", 10.5)),
        height=float(proxy.get("compact_height_mm", 8.0)),
        depth=float(proxy.get("compact_depth_mm", 3.4)),
        points_count=max(2200, int(total_nt * 1.35)),
        seed=int(proxy.get("compact_seed", 304)),
    )
    secondary = proxy.get("secondary_structure", {})
    if secondary.get("model") == "deterministic_schematic_stem_loop":
        structure = structured_stem_loop_path(full_points, total_nt, nt_to_mm, proxy, compact=True)
        full_points = structure["points"]
    segment_models = split_path_by_lengths(full_points, mrna["segments"], nt_to_mm)
    all_points = []
    marker_centers = []
    nt_cursor = 0
    for segment_model in segment_models:
        points = [Vec(*point) for point in segment_model["points"]]
        if all_points:
            all_points.extend([point.to_tuple() for point in points[1:]])
        else:
            all_points.extend([point.to_tuple() for point in points])
        segment = segment_model["segment"]
        for nt in range(nt_cursor + 50, nt_cursor + segment["nt"] + 1, 50):
            fraction = (nt - nt_cursor) / segment["nt"]
            marker_centers.append(points[min(len(points) - 1, int(fraction * (len(points) - 1)))])
        nt_cursor += segment["nt"]
    report = mrna_report_from_segments(
        manifest,
        segment_models,
        "optimization_v3_compact_mrnp_rosette_hairpin_mrna",
        "continuous_mRNP_rosette_hairpin_globule",
    )
    report["model_basis"] = {
        "biological_intent": "schematic compact mRNP-like fold with many local hairpin/rosette turns",
        "not_sequence_specific": True,
        "translation_state": "compact non-translating or storage/stress-associated mRNP, for visual contrast with elongated mRNA",
    }
    if secondary.get("model") == "deterministic_schematic_stem_loop":
        report.update(structure["report"])
        report["variant"] = str(proxy.get("compact_variant", "compact_rosette"))
        report["rna_only"] = True
        return {
            "path": SampledPath(all_points),
            "segments": segment_models,
            "marker_centers": marker_centers,
            "secondary_structure": {
                "base_pairs": structure["base_pairs"],
                "loop_ranges": structure["loop_ranges"],
            },
            "report": report,
        }
    return {
        "path": SampledPath(all_points),
        "segments": segment_models,
        "marker_centers": marker_centers,
        "report": report,
    }
