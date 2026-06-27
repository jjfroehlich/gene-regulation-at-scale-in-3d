#!/usr/bin/env python3
"""Old-proxy-derived DNA/RNA pseudoatom geometry with gap-filling additions."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
import procedural_nucleic_geometry as geom  # noqa: E402


ANGSTROM_TO_MM = 0.04
DNA_GAP_ID = "DNA_PROXY_GAP"
MRNA_GAP_ID = "MRNA_PROXY_GAP"
DNA_SURFACE_VDW_A = 2.4
RNA_SURFACE_VDW_A = 2.2
DNA_BASE_ATOMS_PER_PAIR = 9


def mm_to_angstrom(point: geom.Vec) -> tuple[float, float, float]:
    return (point.x / ANGSTROM_TO_MM, point.y / ANGSTROM_TO_MM, point.z / ANGSTROM_TO_MM)


def atom_record(
    name: str,
    resname: str,
    chain: str,
    resseq: int,
    point: geom.Vec,
    element: str,
    component: str,
) -> dict:
    return {
        "name": name,
        "resname": resname,
        "chain": chain,
        "resseq": resseq,
        "xyz_A": mm_to_angstrom(point),
        "component": component,
        "element": element,
    }


def max_spacing_A(points: list[geom.Vec]) -> float:
    if len(points) < 2:
        return 0.0
    return max((points[i] - points[i - 1]).length / ANGSTROM_TO_MM for i in range(1, len(points)))


def add_midpoint_filled_strand(points: list[geom.Vec], chain: str, component: str) -> list[dict]:
    atoms: list[dict] = []
    for i, point in enumerate(points):
        resseq = i + 1
        atoms.append(atom_record("S", "DNA", chain, resseq, point, "P", component))
        if i < len(points) - 1:
            atoms.append(atom_record("M", "DNA", chain, resseq, point.lerp(points[i + 1], 0.5), "P", component))
    return atoms


def build_dna_gap_asset(manifest: dict) -> dict:
    model = geom.build_dna_model(manifest)
    strand_a: list[geom.Vec] = model["strand_a"]
    strand_b: list[geom.Vec] = model["strand_b"]
    atoms = []
    atoms.extend(add_midpoint_filled_strand(strand_a, "A", "strand_A"))
    atoms.extend(add_midpoint_filled_strand(strand_b, "B", "strand_B"))
    base_spacing = []
    for i, (a, b) in enumerate(zip(strand_a, strand_b), start=1):
        previous = None
        for j in range(DNA_BASE_ATOMS_PER_PAIR):
            fraction = (j + 1) / (DNA_BASE_ATOMS_PER_PAIR + 1)
            point = a.lerp(b, fraction)
            atoms.append(atom_record(f"B{j + 1}", "DBS", "C", i, point, "C", "base_pairs"))
            if previous is not None:
                base_spacing.append((point - previous).length / ANGSTROM_TO_MM)
            previous = point
    gap_strand_spacing = max_spacing_A([geom.Vec(*atom["xyz_A"]) * ANGSTROM_TO_MM for atom in atoms if atom["component"] == "strand_A"])
    old_base_atoms_per_pair = int(model["report"].get("base_atoms_per_pair", 5))
    old_base_spacing_A = (
        (2.0 * model["report"]["strand_center_radius_mm"] / (old_base_atoms_per_pair + 1)) / ANGSTROM_TO_MM
    )
    return {
        "asset_id": DNA_GAP_ID,
        "atoms": atoms,
        "components": {"strand_A": "chain A", "strand_B": "chain B", "base_pairs": "chain C"},
        "surface_vdw_A": {"strand_A": DNA_SURFACE_VDW_A, "strand_B": DNA_SURFACE_VDW_A, "base_pairs": DNA_SURFACE_VDW_A},
        "geometry": {
            **model["report"],
            "style": "legacy_pymol_gap_filled_proxy",
            "gap_fill_changes": [
                "one midpoint pseudoatom inserted between adjacent old strand pseudoatoms",
                f"base-pair pseudoatoms increased from {old_base_atoms_per_pair} to {DNA_BASE_ATOMS_PER_PAIR} per bp",
                "old strand center radius and VDW radius preserved to keep the older visual character",
            ],
            "atom_count": len(atoms),
            "old_base_atoms_per_pair": old_base_atoms_per_pair,
            "new_base_atoms_per_pair": DNA_BASE_ATOMS_PER_PAIR,
            "gap_filled_max_strand_center_spacing_A": gap_strand_spacing,
            "gap_filled_strand_overlap_margin_A": 2.0 * DNA_SURFACE_VDW_A - gap_strand_spacing,
            "old_estimated_base_center_spacing_A": old_base_spacing_A,
            "gap_filled_max_base_center_spacing_A": max(base_spacing) if base_spacing else None,
            "gap_filled_base_overlap_margin_A": 2.0 * DNA_SURFACE_VDW_A - max(base_spacing) if base_spacing else None,
        },
    }


def sample_nt_points(points: list[tuple[float, float, float]], nt_count: int, nt_to_mm: float) -> list[tuple[geom.Vec, geom.Vec]]:
    path = geom.SampledPath(points)
    return [(path.point_at_length(min(path.length, i * nt_to_mm)), path.tangent_at_length(min(path.length, i * nt_to_mm))) for i in range(nt_count)]


def build_mrna_gap_asset(manifest: dict) -> dict:
    model = geom.build_mrna_model(manifest)
    nt_to_mm = manifest["units"]["mrna_nt_to_mm"]
    base_offset_mm = 0.12
    sugar_offset_mm = 0.06
    chains = ["A", "B", "C"]
    resnames = ["U5R", "CDS", "U3R"]
    atoms = []
    p_points_by_component: dict[str, list[geom.Vec]] = {}
    for segment_index, segment_model in enumerate(model["segments"]):
        segment = segment_model["segment"]
        chain = chains[segment_index]
        resname = resnames[segment_index]
        component = {"A": "utr5", "B": "coding", "C": "utr3"}[chain]
        p_points: list[geom.Vec] = []
        sampled = sample_nt_points(segment_model["points"], segment["nt"], nt_to_mm)
        for i, (center, tangent) in enumerate(sampled, start=1):
            normal, binormal = geom.path_frame(tangent)
            theta = 2.0 * math.pi * (i - 1) / 6.0
            radial = normal * math.sin(theta) + binormal * math.cos(theta)
            sugar = center + radial * sugar_offset_mm
            base = center + radial * base_offset_mm
            atoms.append(atom_record("P", resname, chain, i, center, "P", component))
            atoms.append(atom_record("S", resname, chain, i, sugar, "C", component))
            atoms.append(atom_record("B", resname, chain, i, base, "N", component))
            p_points.append(center)
            if i > 1:
                previous_center = sampled[i - 2][0]
                atoms.append(atom_record("M", resname, chain, i - 1, previous_center.lerp(center, 0.5), "P", component))
        p_points_by_component[component] = p_points
    max_old_spacing = max((max_spacing_A(points) for points in p_points_by_component.values() if points), default=0.0)
    midpoint_points = []
    for points in p_points_by_component.values():
        for i in range(1, len(points)):
            midpoint_points.append(points[i - 1])
            midpoint_points.append(points[i - 1].lerp(points[i], 0.5))
        if points:
            midpoint_points.append(points[-1])
    max_gap_spacing = max_spacing_A(midpoint_points)
    return {
        "asset_id": MRNA_GAP_ID,
        "atoms": atoms,
        "components": {"utr5": "chain A", "coding": "chain B", "utr3": "chain C"},
        "surface_vdw_A": {"utr5": RNA_SURFACE_VDW_A, "coding": RNA_SURFACE_VDW_A, "utr3": RNA_SURFACE_VDW_A},
        "geometry": {
            **model["report"],
            "style": "legacy_pymol_gap_filled_proxy",
            "gap_fill_changes": [
                "old P/S/base pseudoatoms retained",
                "one phosphate-like connector pseudoatom inserted between adjacent nucleotides",
                "RNA VDW radius increased modestly from 2.1 A to 2.2 A",
            ],
            "atom_count": len(atoms),
            "old_max_backbone_center_spacing_A": max_old_spacing,
            "gap_filled_max_backbone_center_spacing_A": max_gap_spacing,
            "old_backbone_overlap_margin_A": 2.0 * 2.1 - max_old_spacing,
            "gap_filled_backbone_overlap_margin_A": 2.0 * RNA_SURFACE_VDW_A - max_gap_spacing,
            "surface_vdw_A": RNA_SURFACE_VDW_A,
        },
    }


def baseline_metrics(manifest: dict) -> dict:
    dna = geom.build_dna_model(manifest)
    mrna = geom.build_mrna_model(manifest)
    dna_old_base_atoms = int(dna["report"].get("base_atoms_per_pair", 5))
    dna_old_base_spacing_A = (2.0 * dna["report"]["strand_center_radius_mm"] / (dna_old_base_atoms + 1)) / ANGSTROM_TO_MM
    rna_component_spacings = []
    nt_to_mm = manifest["units"]["mrna_nt_to_mm"]
    for segment_model in mrna["segments"]:
        points = [center for center, _ in sample_nt_points(segment_model["points"], segment_model["segment"]["nt"], nt_to_mm)]
        rna_component_spacings.append(max_spacing_A(points))
    dna_max_strand_spacing_A = max_spacing_A(dna["strand_a"])
    rna_max_spacing_A = max(rna_component_spacings) if rna_component_spacings else 0.0
    return {
        "legacy_original": {
            "dna": {
                "axis_length_mm": dna["report"]["axis_length_mm"],
                "represented_bp": dna["report"]["represented_bp"],
                "estimated_envelope_diameter_mm": dna["report"]["estimated_envelope_diameter_mm"],
                "max_strand_center_spacing_A": dna_max_strand_spacing_A,
                "strand_overlap_margin_A": 2.0 * 2.4 - dna_max_strand_spacing_A,
                "base_atoms_per_pair": dna_old_base_atoms,
                "estimated_base_center_spacing_A": dna_old_base_spacing_A,
                "base_overlap_margin_A": 2.0 * 2.4 - dna_old_base_spacing_A,
            },
            "mrna": {
                "segments": mrna["report"]["segments"],
                "total_measured_mm": mrna["report"]["total_measured_mm"],
                "max_backbone_center_spacing_A": rna_max_spacing_A,
                "backbone_overlap_margin_A": 2.0 * 2.1 - rna_max_spacing_A,
            },
        }
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


def full_actin_mrna_scale_report(manifest: dict) -> dict:
    segments = []
    total = 0.0
    nt_to_mm = manifest["units"]["mrna_nt_to_mm"]
    for segment in manifest["mrna"]["segments"]:
        length = segment["nt"] * nt_to_mm
        segments.append({"name": segment["name"], "nt": segment["nt"], "target_length_mm": length})
        total += length
    return {"segments": segments, "total_nt": manifest["mrna"]["total_nt"], "total_length_mm": total}


def write_generation_report(path: Path, manifest: dict, dna: dict, mrna: dict) -> dict:
    report = {
        "title": "Legacy PyMOL proxy gap-filled pseudoatom assets",
        "units": {
            "coordinate_units": "angstrom in CIF",
            "angstrom_to_mm": ANGSTROM_TO_MM,
            "dna_bp_to_mm": manifest["units"]["dna_bp_to_mm"],
            "rna_nt_to_mm": manifest["units"]["mrna_nt_to_mm"],
        },
        "baseline_metrics": baseline_metrics(manifest),
        "assets": [
            {
                "asset_id": dna["asset_id"],
                "atom_count": len(dna["atoms"]),
                "components": dna["components"],
                "surface_vdw_A": dna["surface_vdw_A"],
                "geometry": dna["geometry"],
            },
            {
                "asset_id": mrna["asset_id"],
                "atom_count": len(mrna["atoms"]),
                "components": mrna["components"],
                "surface_vdw_A": mrna["surface_vdw_A"],
                "geometry": mrna["geometry"],
            },
        ],
        "full_actin_mrna_scale_validation": full_actin_mrna_scale_report(manifest),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
