#!/usr/bin/env python3
"""Generate scale-correct pseudoatom mmCIF files for procedural DNA and mRNA."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import procedural_nucleic_geometry as geom


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = Path(os.environ.get("GENE_SCENE_MANIFEST", ROOT / "config" / "scene_manifest.json"))
OUTPUT_DIR = ROOT / "assets" / "procedural_nucleic_acids"
CIF_DIR = OUTPUT_DIR / "cif"
REPORT_PATH = OUTPUT_DIR / os.environ.get("PROCEDURAL_ASSET_REPORT", "procedural_nucleic_assets_report.json")

VARIANT = os.environ.get("PROCEDURAL_ASSET_VARIANT", "").strip().upper()


def asset_id(base: str) -> str:
    return f"{base}_{VARIANT}" if VARIANT else base


DNA_ID = asset_id("DNA_PROXY")
MRNA_ID = asset_id("MRNA_PROXY")
MRNA_COMPACT_ID = asset_id("MRNA_COMPACT_PROXY")
DNA_STRAND_CONNECTORS_PER_STEP = 2
DNA_BASE_STACK_CONNECTORS_PER_STEP = 2
RNA_BACKBONE_CONNECTORS_PER_STEP = 2
RNA_SUGAR_CONNECTORS_PER_STEP = 1
DNA_SURFACE_EXPORT_VDW_A = 2.4
RNA_SURFACE_EXPORT_VDW_A = 2.2


def mm_to_angstrom(point: geom.Vec, angstrom_to_mm: float) -> tuple[float, float, float]:
    return (point.x / angstrom_to_mm, point.y / angstrom_to_mm, point.z / angstrom_to_mm)


def max_spacing_A(points: list[geom.Vec], angstrom_to_mm: float) -> float:
    if len(points) < 2:
        return 0.0
    return max((points[i] - points[i - 1]).length / angstrom_to_mm for i in range(1, len(points)))


def add_atom(
    atoms: list[dict],
    name: str,
    resname: str,
    chain: str,
    resseq: int,
    point: geom.Vec,
    element: str,
    scale: float,
) -> None:
    atoms.append(
        {
            "name": name,
            "resname": resname,
            "chain": chain,
            "resseq": resseq,
            "xyz_A": mm_to_angstrom(point, scale),
            "element": element,
        }
    )


def with_connectors(points: list[geom.Vec], connectors_per_step: int) -> list[geom.Vec]:
    if connectors_per_step <= 0 or len(points) < 2:
        return list(points)
    filled: list[geom.Vec] = []
    for index, point in enumerate(points):
        filled.append(point)
        if index < len(points) - 1:
            next_point = points[index + 1]
            for connector in range(1, connectors_per_step + 1):
                filled.append(point.lerp(next_point, connector / (connectors_per_step + 1)))
    return filled


def pdb_atom_line(serial: int, name: str, resname: str, chain: str, resseq: int, xyz_A: tuple[float, float, float], element: str) -> str:
    x, y, z = xyz_A
    return (
        f"ATOM  {serial:5d} {name:<4} {resname:>3} {chain:1}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}\n"
    )


def write_pdb(path: Path, atoms: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("REMARK Procedural pseudoatom nucleic acid surface proxy\n")
        for serial, atom in enumerate(atoms, start=1):
            handle.write(
                pdb_atom_line(
                    serial,
                    atom["name"],
                    atom["resname"],
                    atom["chain"],
                    atom["resseq"],
                    atom["xyz_A"],
                    atom["element"],
                )
            )
        handle.write("END\n")


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


def dna_atoms(manifest: dict, model: dict) -> list[dict]:
    scale = manifest["units"]["angstrom_to_mm"]
    atoms = []
    sugar_a = model.get("sugar_a", [])
    sugar_b = model.get("sugar_b", [])
    for i, point in enumerate(model["strand_a"], start=1):
        add_atom(atoms, "P", "DNA", "A", i, point, "P", scale)
        if sugar_a:
            add_atom(atoms, "S", "DNA", "A", i, sugar_a[i - 1], "C", scale)
        if i < len(model["strand_a"]):
            next_point = model["strand_a"][i]
            for connector in range(1, DNA_STRAND_CONNECTORS_PER_STEP + 1):
                add_atom(atoms, f"M{connector}", "DNA", "A", i, point.lerp(next_point, connector / (DNA_STRAND_CONNECTORS_PER_STEP + 1)), "P", scale)
            if sugar_a:
                next_sugar = sugar_a[i]
                add_atom(atoms, "SM1", "DNA", "A", i, sugar_a[i - 1].lerp(next_sugar, 0.5), "C", scale)
    for i, point in enumerate(model["strand_b"], start=1):
        add_atom(atoms, "P", "DNA", "B", i, point, "P", scale)
        if sugar_b:
            add_atom(atoms, "S", "DNA", "B", i, sugar_b[i - 1], "C", scale)
        if i < len(model["strand_b"]):
            next_point = model["strand_b"][i]
            for connector in range(1, DNA_STRAND_CONNECTORS_PER_STEP + 1):
                add_atom(atoms, f"M{connector}", "DNA", "B", i, point.lerp(next_point, connector / (DNA_STRAND_CONNECTORS_PER_STEP + 1)), "P", scale)
            if sugar_b:
                next_sugar = sugar_b[i]
                add_atom(atoms, "SM1", "DNA", "B", i, sugar_b[i - 1].lerp(next_sugar, 0.5), "C", scale)
    per_pair = model["report"]["base_atoms_per_pair"]
    rows = model.get("base_pair_rows") or [model["base_pair_atoms"][i:i + per_pair] for i in range(0, len(model["base_pair_atoms"]), per_pair)]
    for row_index, row in enumerate(rows, start=1):
        for base_index, point in enumerate(row, start=1):
            add_atom(atoms, f"B{base_index}", "DBS", "C", row_index, point, "N" if base_index % 2 else "C", scale)
            if row_index < len(rows):
                next_point = rows[row_index][base_index - 1]
                for connector in range(1, DNA_BASE_STACK_CONNECTORS_PER_STEP + 1):
                    add_atom(
                        atoms,
                        f"L{base_index}",
                        "DBS",
                        "C",
                        row_index,
                        point.lerp(next_point, connector / (DNA_BASE_STACK_CONNECTORS_PER_STEP + 1)),
                        "C",
                        scale,
                    )
    return atoms


def sample_nt_points(points: list[tuple[float, float, float]], nt_count: int, nt_to_mm: float) -> list[tuple[geom.Vec, geom.Vec]]:
    path = geom.SampledPath(points)
    sampled = []
    for i in range(nt_count):
        distance = min(path.length, i * nt_to_mm)
        sampled.append((path.point_at_length(distance), path.tangent_at_length(distance)))
    return sampled


def mrna_atoms(manifest: dict, model: dict) -> list[dict]:
    scale = manifest["units"]["angstrom_to_mm"]
    nt_to_mm = manifest["units"]["mrna_nt_to_mm"]
    proxy = manifest.get("procedural_nucleic_surfaces", {}).get("mrna", {})
    base_offset_mm = proxy.get("base_offset_mm", 0.16)
    base_mid_offset_mm = proxy.get("base_mid_offset_mm", 0.115)
    sugar_offset_mm = proxy.get("sugar_offset_mm", 0.065)
    chains = ["A", "B", "C"]
    resnames = ["U5R", "CDS", "U3R"]
    atoms = []
    for segment_index, segment_model in enumerate(model["segments"]):
        segment = segment_model["segment"]
        chain = chains[segment_index]
        resname = resnames[segment_index]
        sampled = sample_nt_points(segment_model["points"], segment["nt"], nt_to_mm)
        sugar_points: list[geom.Vec] = []
        for i, (center, tangent) in enumerate(sampled, start=1):
            normal, binormal = geom.path_frame(tangent)
            theta = 2.0 * math.pi * (i - 1) / 5.5
            radial = normal * math.sin(theta) + binormal * math.cos(theta)
            sugar = center + radial * sugar_offset_mm
            base = center + radial * base_offset_mm + tangent * (0.018 * math.sin(theta * 0.7))
            base_mid = center + radial * base_mid_offset_mm
            add_atom(atoms, "P", resname, chain, i, center, "P", scale)
            add_atom(atoms, "S", resname, chain, i, sugar, "C", scale)
            add_atom(atoms, "B", resname, chain, i, base, "N", scale)
            add_atom(atoms, "BM", resname, chain, i, base_mid, "C", scale)
            sugar_points.append(sugar)
            if i < len(sampled):
                next_center = sampled[i][0]
                for connector in range(1, RNA_BACKBONE_CONNECTORS_PER_STEP + 1):
                    add_atom(atoms, f"M{connector}", resname, chain, i, center.lerp(next_center, connector / (RNA_BACKBONE_CONNECTORS_PER_STEP + 1)), "P", scale)
        for i, sugar in enumerate(sugar_points[:-1], start=1):
            next_sugar = sugar_points[i]
            for connector in range(1, RNA_SUGAR_CONNECTORS_PER_STEP + 1):
                add_atom(atoms, f"SM{connector}", resname, chain, i, sugar.lerp(next_sugar, connector / (RNA_SUGAR_CONNECTORS_PER_STEP + 1)), "C", scale)
    return atoms


def dna_repair_report(manifest: dict, model: dict) -> dict:
    scale = manifest["units"]["angstrom_to_mm"]
    per_pair = model["report"]["base_atoms_per_pair"]
    rows = [model["base_pair_atoms"][i:i + per_pair] for i in range(0, len(model["base_pair_atoms"]), per_pair)]
    base_stack_spacings = []
    for base_index in range(per_pair):
        base_stack_spacings.append(max_spacing_A([row[base_index] for row in rows], scale))
    return {
        "style": "procedural_pymol_surface_proxy_repaired",
        "strand_connectors_per_bp_step": DNA_STRAND_CONNECTORS_PER_STEP,
        "base_stack_connectors_per_bp_step": DNA_BASE_STACK_CONNECTORS_PER_STEP,
        "original_max_strand_center_spacing_A": max_spacing_A(model["strand_a"], scale),
        "repaired_max_strand_center_spacing_A": max_spacing_A(with_connectors(model["strand_a"], DNA_STRAND_CONNECTORS_PER_STEP), scale),
        "original_max_base_stack_spacing_A": max(base_stack_spacings) if base_stack_spacings else None,
        "repaired_max_base_stack_spacing_A": (max(base_stack_spacings) / (DNA_BASE_STACK_CONNECTORS_PER_STEP + 1)) if base_stack_spacings else None,
        "surface_export_vdw_A": DNA_SURFACE_EXPORT_VDW_A,
    }


def mrna_repair_report(manifest: dict, model: dict) -> dict:
    scale = manifest["units"]["angstrom_to_mm"]
    nt_to_mm = manifest["units"]["mrna_nt_to_mm"]
    segment_reports = []
    for segment_model in model["segments"]:
        segment = segment_model["segment"]
        centers = [center for center, _ in sample_nt_points(segment_model["points"], segment["nt"], nt_to_mm)]
        segment_reports.append(
            {
                "name": segment["name"],
                "original_max_backbone_spacing_A": max_spacing_A(centers, scale),
                "repaired_max_backbone_spacing_A": max_spacing_A(with_connectors(centers, RNA_BACKBONE_CONNECTORS_PER_STEP), scale),
            }
        )
    return {
        "style": "procedural_pymol_surface_proxy_repaired",
        "backbone_connectors_per_nt_step": RNA_BACKBONE_CONNECTORS_PER_STEP,
        "sugar_connectors_per_nt_step": RNA_SUGAR_CONNECTORS_PER_STEP,
        "surface_export_vdw_A": RNA_SURFACE_EXPORT_VDW_A,
        "segments": segment_reports,
    }


def main() -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    CIF_DIR.mkdir(parents=True, exist_ok=True)

    dna_model = geom.build_dna_model(manifest)
    mrna_model = geom.build_mrna_model(manifest)
    compact_mrna_model = geom.build_compact_mrna_model(manifest)
    dna_atom_records = dna_atoms(manifest, dna_model)
    mrna_atom_records = mrna_atoms(manifest, mrna_model)
    compact_mrna_atom_records = mrna_atoms(manifest, compact_mrna_model)

    dna_cif_path = CIF_DIR / f"{DNA_ID}.cif"
    mrna_cif_path = CIF_DIR / f"{MRNA_ID}.cif"
    compact_mrna_cif_path = CIF_DIR / f"{MRNA_COMPACT_ID}.cif"
    write_cif(dna_cif_path, DNA_ID, dna_atom_records)
    write_cif(mrna_cif_path, MRNA_ID, mrna_atom_records)
    write_cif(compact_mrna_cif_path, MRNA_COMPACT_ID, compact_mrna_atom_records)

    report = {
        "title": "Procedural nucleic-acid pseudoatom assets",
        "coordinate_units": "angstrom",
        "note": "mmCIF is used for PyMOL export because full-scene Angstrom coordinates exceed fixed-width PDB coordinate columns.",
        "angstrom_to_mm": manifest["units"]["angstrom_to_mm"],
        "assets": [
            {
                "id": DNA_ID,
                "path": str(dna_cif_path.relative_to(ROOT)).replace("\\", "/"),
                "atom_count": len(dna_atom_records),
                "components": {"strand_A": "chain A", "strand_B": "chain B", "base_pairs": "chain C"},
                "geometry": {**dna_model["report"], **dna_repair_report(manifest, dna_model)},
            },
            {
                "id": MRNA_ID,
                "path": str(mrna_cif_path.relative_to(ROOT)).replace("\\", "/"),
                "atom_count": len(mrna_atom_records),
                "components": {"utr5": "chain A", "coding": "chain B", "utr3": "chain C"},
                "geometry": {**mrna_model["report"], **mrna_repair_report(manifest, mrna_model)},
            },
            {
                "id": MRNA_COMPACT_ID,
                "path": str(compact_mrna_cif_path.relative_to(ROOT)).replace("\\", "/"),
                "atom_count": len(compact_mrna_atom_records),
                "components": {"utr5": "chain A", "coding": "chain B", "utr3": "chain C"},
                "geometry": {**compact_mrna_model["report"], **mrna_repair_report(manifest, compact_mrna_model)},
            },
        ],
    }
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
