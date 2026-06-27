#!/usr/bin/env python3
"""Analyze DNA/RNA calibrator mmCIF files for procedural scale checks."""

from __future__ import annotations

import json
import math
import re
import shlex
from collections import defaultdict
from pathlib import Path

from procedural_nucleic_geometry import Vec


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "config" / "scene_manifest.json"
RCSB_DIR = ROOT / "assets" / "rcsb"
OUTPUT_DIR = ROOT / "experiments" / "procedural_nucleic_acids" / "outputs"
REPORT_PATH = OUTPUT_DIR / "nucleic_acid_calibrator_analysis.json"

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


def load_manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


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
        if comp not in NUCLEIC_COMPS:
            continue
        atom_by_name = {item["atom"]: item for item in items}
        chosen = None
        for atom_name in ("P", "C4'", "C3'", "C1'"):
            if atom_name in atom_by_name:
                chosen = atom_by_name[atom_name]
                break
        if chosen is None:
            chosen = {
                "x": sum(item["x"] for item in items) / len(items),
                "y": sum(item["y"] for item in items) / len(items),
                "z": sum(item["z"] for item in items) / len(items),
            }
        points.append({"chain": chain, "seq": seq, "comp": comp, "pos_A": Vec(chosen["x"], chosen["y"], chosen["z"])})
    return points


def natural_seq_key(seq: str) -> tuple[int, str]:
    match = re.search(r"-?\d+", seq)
    if match:
        return (int(match.group(0)), seq)
    return (10**9, seq)


def mean_vector(vectors: list[Vec]) -> Vec:
    total = Vec(0.0, 0.0, 0.0)
    for vector in vectors:
        total = total + vector
    return total * (1.0 / max(1, len(vectors)))


def dot(a: Vec, b: Vec) -> float:
    return a.x * b.x + a.y * b.y + a.z * b.z


def principal_axis(vectors: list[Vec]) -> Vec:
    center = mean_vector(vectors)
    cov = [[0.0, 0.0, 0.0] for _ in range(3)]
    for vector in vectors:
        d = vector - center
        values = (d.x, d.y, d.z)
        for i in range(3):
            for j in range(3):
                cov[i][j] += values[i] * values[j]
    axis = Vec(1.0, 0.37, 0.19).normalized()
    for _ in range(40):
        next_axis = Vec(
            cov[0][0] * axis.x + cov[0][1] * axis.y + cov[0][2] * axis.z,
            cov[1][0] * axis.x + cov[1][1] * axis.y + cov[1][2] * axis.z,
            cov[2][0] * axis.x + cov[2][1] * axis.y + cov[2][2] * axis.z,
        )
        if next_axis.length < 1e-9:
            return Vec(1.0, 0.0, 0.0)
        axis = next_axis.normalized()
    return axis


def span_and_diameter(points: list[Vec]) -> tuple[float, float]:
    if len(points) < 2:
        return 0.0, 0.0
    axis = principal_axis(points)
    center = mean_vector(points)
    projections = [dot(point - center, axis) for point in points]
    span = max(projections) - min(projections)
    max_radius = 0.0
    for point, projection in zip(points, projections):
        closest = center + axis * projection
        max_radius = max(max_radius, (point - closest).length)
    return span, 2.0 * max_radius


def chain_stats(points: list[dict]) -> list[dict]:
    chains = defaultdict(list)
    for point in points:
        chains[point["chain"]].append(point)
    stats = []
    for chain, chain_points in sorted(chains.items()):
        ordered = sorted(chain_points, key=lambda item: natural_seq_key(item["seq"]))
        positions = [item["pos_A"] for item in ordered]
        contour = sum((positions[i] - positions[i - 1]).length for i in range(1, len(positions)))
        span, diameter = span_and_diameter(positions)
        stats.append(
            {
                "chain": chain,
                "residues": len(ordered),
                "first_seq": ordered[0]["seq"],
                "last_seq": ordered[-1]["seq"],
                "contour_A": contour,
                "mean_contour_step_A": contour / max(1, len(ordered) - 1),
                "axis_span_A": span,
                "axis_step_A": span / max(1, len(ordered) - 1),
                "envelope_diameter_A": diameter,
                "bases": sorted({item["comp"] for item in ordered}),
            }
        )
    return stats


def analyze_one(pdb_id: str, expected_type: str, units: dict) -> dict:
    path = RCSB_DIR / f"{pdb_id}.cif"
    if not path.exists():
        return {"pdb_id": pdb_id, "expected_type": expected_type, "status": "missing", "path": str(path)}
    atoms = parse_atom_site(path)
    nucleic_atoms = [atom for atom in atoms if atom["comp"] in NUCLEIC_COMPS]
    protein_atoms = [atom for atom in atoms if atom["comp"] not in NUCLEIC_COMPS and atom["element"] != "H"]
    points = residue_points(atoms)
    positions = [point["pos_A"] for point in points]
    span_A, diameter_A = span_and_diameter(positions) if positions else (0.0, 0.0)
    chains = chain_stats(points)
    longest_chain = max(chains, key=lambda item: item["residues"], default=None)
    return {
        "pdb_id": pdb_id,
        "expected_type": expected_type,
        "status": "ok" if points else "no_nucleic_residue_points",
        "path": str(path.relative_to(ROOT)).replace("\\", "/"),
        "atom_count": len(atoms),
        "nucleic_atom_count": len(nucleic_atoms),
        "protein_or_other_atom_count": len(protein_atoms),
        "nucleic_residue_count": len(points),
        "nucleic_chain_count": len(chains),
        "axis_span_A": span_A,
        "axis_span_mm": span_A * units["angstrom_to_mm"],
        "envelope_diameter_A": diameter_A,
        "envelope_diameter_mm": diameter_A * units["angstrom_to_mm"],
        "chains": chains,
        "longest_chain": longest_chain,
    }


def main() -> int:
    manifest = load_manifest()
    units = manifest["units"]
    entries = []
    for expected_type, assets in manifest.get("nucleic_acid_calibrators", {}).items():
        for asset in assets:
            entry = analyze_one(asset["pdb_id"].upper(), expected_type, units)
            entry["source"] = asset.get("source")
            entry["note"] = asset.get("note")
            entries.append(entry)
    report = {
        "title": "Nucleic-acid PDB calibrator analysis",
        "units": units,
        "scale_targets": {
            "dna_bp_spacing_mm": units["dna_bp_to_mm"],
            "mrna_nt_contour_mm": units["mrna_nt_to_mm"],
            "dna_envelope_diameter_mm": manifest["procedural_nucleic_acids"]["dna"]["envelope_diameter_nm"] * units["nm_to_mm"],
        },
        "entries": entries,
    }
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if all(entry["status"] == "ok" for entry in entries) else 1


if __name__ == "__main__":
    raise SystemExit(main())
