#!/usr/bin/env python3
"""Shared geometry for the detailed DNA/RNA route comparison experiment."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path


ANGSTROM_TO_MM = 0.04
DNA_BP_TO_MM = 0.136
RNA_NT_TO_MM = 0.12
DNA_BP_PER_TURN = 10.5
DNA_VISUAL_BP = 120
RNA_VISUAL_SEGMENTS = [
    {"name": "5' UTR", "component": "utr5", "chain": "A", "nt": 40},
    {"name": "coding", "component": "coding", "chain": "B", "nt": 120},
    {"name": "3' UTR", "component": "utr3", "chain": "C", "nt": 80},
]
FULL_ACTIN_SEGMENTS = [
    {"name": "5' UTR", "nt": 84},
    {"name": "coding sequence", "nt": 1125},
    {"name": "3' UTR", "nt": 643},
]


@dataclass(frozen=True)
class Vec:
    x: float
    y: float
    z: float

    def __add__(self, other: "Vec") -> "Vec":
        return Vec(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other: "Vec") -> "Vec":
        return Vec(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, scalar: float) -> "Vec":
        return Vec(self.x * scalar, self.y * scalar, self.z * scalar)

    __rmul__ = __mul__

    @property
    def length(self) -> float:
        return math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def normalized(self) -> "Vec":
        length = self.length
        if length < 1e-12:
            return Vec(1.0, 0.0, 0.0)
        return self * (1.0 / length)

    def dot(self, other: "Vec") -> float:
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other: "Vec") -> "Vec":
        return Vec(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )

    def lerp(self, other: "Vec", t: float) -> "Vec":
        return self + (other - self) * t

    def to_tuple(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)


class SampledPath:
    def __init__(self, points: list[Vec]):
        self.points = points
        self.cumulative = [0.0]
        for i in range(1, len(points)):
            self.cumulative.append(self.cumulative[-1] + (points[i] - points[i - 1]).length)
        self.length = self.cumulative[-1] if self.cumulative else 0.0

    def point_at(self, distance: float) -> Vec:
        distance = min(max(distance, 0.0), self.length)
        for i in range(1, len(self.cumulative)):
            if self.cumulative[i] >= distance:
                span = self.cumulative[i] - self.cumulative[i - 1]
                t = 0.0 if span == 0.0 else (distance - self.cumulative[i - 1]) / span
                return self.points[i - 1].lerp(self.points[i], t)
        return self.points[-1]

    def tangent_at(self, distance: float) -> Vec:
        distance = min(max(distance, 0.0), self.length)
        for i in range(1, len(self.cumulative)):
            if self.cumulative[i] >= distance:
                return (self.points[i] - self.points[i - 1]).normalized()
        return (self.points[-1] - self.points[-2]).normalized()


def path_frame(tangent: Vec) -> tuple[Vec, Vec]:
    tangent = tangent.normalized()
    normal = Vec(-tangent.y, tangent.x, 0.0)
    if normal.length < 1e-6:
        normal = Vec(0.0, 1.0, 0.0)
    binormal = tangent.cross(normal).normalized()
    return normal.normalized(), binormal


def polyline_length(points: list[Vec]) -> float:
    return sum((points[i] - points[i - 1]).length for i in range(1, len(points)))


def make_rna_path(nt: int) -> SampledPath:
    target_length = nt * RNA_NT_TO_MM

    def points_for_span(span: float) -> list[Vec]:
        points = []
        for i in range(nt + 1):
            t = i / nt
            x = span * t
            y = 0.62 * math.sin(2.0 * math.pi * 2.8 * t) + 0.18 * math.sin(2.0 * math.pi * 9.0 * t)
            z = 0.28 * math.sin(2.0 * math.pi * 1.3 * t + 0.5)
            points.append(Vec(x, y, z))
        return points

    low = 0.0
    high = target_length
    while polyline_length(points_for_span(high)) < target_length:
        high *= 1.5
    for _ in range(50):
        mid = (low + high) * 0.5
        if polyline_length(points_for_span(mid)) < target_length:
            low = mid
        else:
            high = mid
    return SampledPath(points_for_span(high))


def make_dna_path(bp: int) -> SampledPath:
    target_length = bp * DNA_BP_TO_MM
    samples = bp * 4

    def points_for_span(span: float) -> list[Vec]:
        points = []
        for i in range(samples + 1):
            t = i / samples
            x = span * t
            y = 0.28 * math.sin(2.0 * math.pi * t)
            z = 0.10 * math.sin(2.0 * math.pi * 1.7 * t)
            points.append(Vec(x, y, z))
        return points

    low = 0.0
    high = target_length
    while polyline_length(points_for_span(high)) < target_length:
        high *= 1.5
    for _ in range(50):
        mid = (low + high) * 0.5
        if polyline_length(points_for_span(mid)) < target_length:
            low = mid
        else:
            high = mid
    return SampledPath(points_for_span(high))


def mm_to_angstrom(point: Vec) -> tuple[float, float, float]:
    return (point.x / ANGSTROM_TO_MM, point.y / ANGSTROM_TO_MM, point.z / ANGSTROM_TO_MM)


def atom_record(
    name: str,
    resname: str,
    chain: str,
    resseq: int,
    point_mm: Vec,
    element: str,
    component: str,
    radius_mm: float,
) -> dict:
    return {
        "name": name,
        "resname": resname,
        "chain": chain,
        "resseq": resseq,
        "xyz_mm": point_mm.to_tuple(),
        "xyz_A": mm_to_angstrom(point_mm),
        "element": element,
        "component": component,
        "radius_mm": radius_mm,
    }


def build_dna_detail_atoms(bp: int = DNA_VISUAL_BP) -> dict:
    path = make_dna_path(bp)
    strand_vdw_A = 3.2
    base_vdw_A = 2.45
    strand_radius_mm = strand_vdw_A * ANGSTROM_TO_MM
    base_radius_mm = base_vdw_A * ANGSTROM_TO_MM
    strand_center_radius_mm = 0.4 - strand_radius_mm
    strand_subsamples_per_bp = 3
    base_atoms_per_pair = 9
    strand_a = []
    strand_b = []
    base_atoms = []
    strand_spacing = []

    for i in range(bp):
        for sub in range(strand_subsamples_per_bp):
            f = i + sub / strand_subsamples_per_bp
            distance = min(path.length, f * DNA_BP_TO_MM)
            center = path.point_at(distance)
            tangent = path.tangent_at(distance)
            normal, binormal = path_frame(tangent)
            theta = 2.0 * math.pi * f / DNA_BP_PER_TURN
            radial = normal * math.sin(theta) + binormal * math.cos(theta)
            a = center + radial * strand_center_radius_mm
            b = center - radial * strand_center_radius_mm
            resseq = i + 1
            strand_a.append(atom_record("S", "DNA", "A", resseq, a, "P", "strand_A", strand_radius_mm))
            strand_b.append(atom_record("S", "DNA", "B", resseq, b, "P", "strand_B", strand_radius_mm))

    for i in range(1, len(strand_a)):
        a0 = Vec(*strand_a[i - 1]["xyz_mm"])
        a1 = Vec(*strand_a[i]["xyz_mm"])
        strand_spacing.append((a1 - a0).length / ANGSTROM_TO_MM)

    for i in range(bp):
        f = i + 0.5
        distance = min(path.length, f * DNA_BP_TO_MM)
        center = path.point_at(distance)
        tangent = path.tangent_at(distance)
        normal, binormal = path_frame(tangent)
        theta = 2.0 * math.pi * f / DNA_BP_PER_TURN
        radial = normal * math.sin(theta) + binormal * math.cos(theta)
        a = center + radial * (strand_center_radius_mm * 0.92)
        b = center - radial * (strand_center_radius_mm * 0.92)
        for j in range(base_atoms_per_pair):
            fraction = (j + 1) / (base_atoms_per_pair + 1)
            offset = tangent * (0.018 * math.sin((j + 1) * math.pi / 2.0))
            point = a.lerp(b, fraction) + offset
            base_atoms.append(atom_record(f"B{j + 1}", "DBS", "C", i + 1, point, "C", "base_pairs", base_radius_mm))

    atoms = strand_a + strand_b + base_atoms
    return {
        "asset_id": "DNA_DETAIL_PROXY",
        "atoms": atoms,
        "components": {"strand_A": "chain A", "strand_B": "chain B", "base_pairs": "chain C"},
        "surface_vdw_A": {"strand_A": strand_vdw_A, "strand_B": strand_vdw_A, "base_pairs": base_vdw_A},
        "geometry": {
            "bp": bp,
            "axis_length_mm": path.length,
            "bp_spacing_mm": DNA_BP_TO_MM,
            "bp_per_turn": DNA_BP_PER_TURN,
            "strand_center_radius_mm": strand_center_radius_mm,
            "strand_vdw_A": strand_vdw_A,
            "base_vdw_A": base_vdw_A,
            "target_envelope_diameter_mm": 0.8,
            "estimated_envelope_diameter_mm": 2.0 * (strand_center_radius_mm + strand_radius_mm),
            "strand_subsamples_per_bp": strand_subsamples_per_bp,
            "base_atoms_per_pair": base_atoms_per_pair,
            "max_strand_center_spacing_A": max(strand_spacing) if strand_spacing else None,
            "strand_overlap_margin_A": (2.0 * strand_vdw_A - max(strand_spacing)) if strand_spacing else None,
        },
    }


def build_rna_detail_atoms() -> dict:
    segments = []
    atoms = []
    cursor = 0
    vdw_A = 2.5
    radius_mm = vdw_A * ANGSTROM_TO_MM
    total_nt = sum(segment["nt"] for segment in RNA_VISUAL_SEGMENTS)
    path = make_rna_path(total_nt)
    spacing = []

    for segment in RNA_VISUAL_SEGMENTS:
        segment_atoms = []
        for local_nt in range(segment["nt"]):
            global_nt = cursor + local_nt
            distance = min(path.length, global_nt * RNA_NT_TO_MM)
            center = path.point_at(distance)
            tangent = path.tangent_at(distance)
            normal, binormal = path_frame(tangent)
            theta = 2.0 * math.pi * global_nt / 6.0
            radial = normal * math.sin(theta) + binormal * math.cos(theta)
            side = binormal * (0.035 * math.sin(theta * 1.7))
            resseq = local_nt + 1
            component = segment["component"]
            chain = segment["chain"]
            resname = {"utr5": "U5R", "coding": "CDS", "utr3": "U3R"}[component]
            centers = [
                ("P", center, "P"),
                ("S", center + radial * 0.095 + side, "C"),
                ("B1", center + radial * 0.205 + side, "N"),
                ("B2", center + radial * 0.270 - side, "N"),
                ("B3", center + radial * 0.205 - binormal * 0.055, "C"),
            ]
            for name, point, element in centers:
                record = atom_record(name, resname, chain, resseq, point, element, component, radius_mm)
                atoms.append(record)
                segment_atoms.append(record)
        segments.append(
            {
                "name": segment["name"],
                "component": segment["component"],
                "chain": segment["chain"],
                "nt": segment["nt"],
                "target_length_mm": segment["nt"] * RNA_NT_TO_MM,
            }
        )
        cursor += segment["nt"]

    p_centers = [Vec(*atom["xyz_mm"]) for atom in atoms if atom["name"] == "P"]
    for i in range(1, len(p_centers)):
        spacing.append((p_centers[i] - p_centers[i - 1]).length / ANGSTROM_TO_MM)

    return {
        "asset_id": "RNA_DETAIL_PROXY",
        "atoms": atoms,
        "components": {"utr5": "chain A", "coding": "chain B", "utr3": "chain C"},
        "surface_vdw_A": {"utr5": vdw_A, "coding": vdw_A, "utr3": vdw_A},
        "geometry": {
            "nt": total_nt,
            "segments": segments,
            "target_contour_mm": total_nt * RNA_NT_TO_MM,
            "measured_path_mm": path.length,
            "nt_spacing_mm": RNA_NT_TO_MM,
            "pseudoatoms_per_nt": 5,
            "surface_vdw_A": vdw_A,
            "max_backbone_center_spacing_A": max(spacing) if spacing else None,
            "backbone_overlap_margin_A": (2.0 * vdw_A - max(spacing)) if spacing else None,
        },
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


def full_actin_mrna_scale_report() -> dict:
    segments = []
    total = 0.0
    for segment in FULL_ACTIN_SEGMENTS:
        length = segment["nt"] * RNA_NT_TO_MM
        segments.append({"name": segment["name"], "nt": segment["nt"], "target_length_mm": length})
        total += length
    return {"segments": segments, "total_nt": sum(item["nt"] for item in FULL_ACTIN_SEGMENTS), "total_length_mm": total}


def write_generation_report(path: Path, dna: dict, rna: dict) -> dict:
    report = {
        "title": "Detailed procedural nucleic-acid proxy generation",
        "units": {
            "coordinate_units": "angstrom in CIF, millimeter in Blender reports",
            "angstrom_to_mm": ANGSTROM_TO_MM,
            "dna_bp_to_mm": DNA_BP_TO_MM,
            "rna_nt_to_mm": RNA_NT_TO_MM,
        },
        "assets": [
            {
                "asset_id": dna["asset_id"],
                "atom_count": len(dna["atoms"]),
                "components": dna["components"],
                "surface_vdw_A": dna["surface_vdw_A"],
                "geometry": dna["geometry"],
            },
            {
                "asset_id": rna["asset_id"],
                "atom_count": len(rna["atoms"]),
                "components": rna["components"],
                "surface_vdw_A": rna["surface_vdw_A"],
                "geometry": rna["geometry"],
            },
        ],
        "full_actin_mrna_scale_validation": full_actin_mrna_scale_report(),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
