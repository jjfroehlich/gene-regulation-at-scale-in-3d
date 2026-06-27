#!/usr/bin/env python3
"""Export low-quality PyMOL molecular surfaces for scene PDB assets.

Run inside PyMOL:
    pymol -cq -d "run C:/path/to/scripts/export_pymol_surface_assets.py"

Use env var PYMOL_SURFACE_PDBS=1J5E,1JJ2 to export a subset.
For legacy 4V5D exports, only one 70S copy is exported by default (chains A*/B*).
"""

from __future__ import annotations

import json
import os
from pathlib import Path

from pymol import cmd


ROOT = Path(os.environ.get("GENE_SCENE_ROOT", Path.cwd())).resolve()
MANIFEST_PATH = ROOT / "config" / "scene_manifest.json"
RCSB_DIR = ROOT / "assets" / "rcsb"
OUTPUT_DIR = ROOT / "assets" / "pymol_exports" / "surface_assets"

GENERIC_COMPONENTS = {
    "protein": "polymer.protein",
    "nucleic": "polymer.nucleic",
}


def safe_path(path: Path) -> str:
    return path.as_posix()


def selected_pdbs() -> set[str] | None:
    raw = os.environ.get("PYMOL_SURFACE_PDBS", "").strip()
    if not raw:
        return None
    return {item.strip().upper() for item in raw.split(",") if item.strip()}


def configure() -> None:
    cmd.set("surface_quality", 0)
    cmd.set("surface_smooth_edges", 1)
    cmd.set("solvent_radius", 1.4)
    cmd.set("specular", 0)
    cmd.set("ambient", 0.72)
    cmd.set("direct", 0.32)
    cmd.set("shininess", 0)
    cmd.set("two_sided_lighting", 1)


def chain_expression(chains: list[str]) -> str:
    if not chains:
        return "none"
    return "(" + " or ".join(f"chain {chain}" for chain in chains) + ")"


def residue_count(selection: str) -> int:
    seen = set()
    cmd.iterate(selection, "seen.add((chain, resi, resn))", space={"seen": seen})
    return len(seen)


def export_selection(object_name: str, component: str, selection_body: str, out_path: Path) -> dict:
    selection = f"{object_name}_{component}_sel"
    if out_path.exists():
        out_path.unlink()
    cmd.select(selection, f"{object_name} and ({selection_body})")
    atom_count = cmd.count_atoms(selection)
    residue_total = residue_count(selection) if atom_count else 0
    if atom_count:
        cmd.hide("everything", object_name)
        cmd.show("surface", selection)
        cmd.save(safe_path(out_path), selection, state=1)
    cmd.delete(selection)
    return {
        "component": component,
        "selection": selection_body,
        "atom_count": atom_count,
        "residue_count": residue_total,
        "exists": out_path.exists(),
        "bytes": out_path.stat().st_size if out_path.exists() else 0,
    }


def export_generic_components(object_name: str, pdb_id: str, out_dir: Path) -> list[dict]:
    exports = []
    for component, selection_body in GENERIC_COMPONENTS.items():
        obj_path = out_dir / f"{pdb_id}_surface_{component}.obj"
        report = export_selection(object_name, component, selection_body, obj_path)
        report["obj"] = str(obj_path.relative_to(ROOT)).replace("\\", "/")
        exports.append(report)
    return exports


def first_ribosome_copy_chains(chains: list[str]) -> list[str]:
    copy = os.environ.get("PYMOL_RIBOSOME_COPY", "first").strip().lower()
    if copy == "all":
        return chains
    if copy == "second":
        prefixes = {"C", "D"}
    else:
        prefixes = {"A", "B"}
    return [chain for chain in chains if chain and chain[0].upper() in prefixes]


def classify_ribosome_nucleic_chains(object_name: str, copy_selection: str) -> dict[str, list[str]]:
    chains = cmd.get_chains(f"{object_name} and polymer.nucleic and {copy_selection}")
    groups = {"rRNA": [], "mRNA": [], "tRNA": []}
    for chain in chains:
        selection = f"{object_name} and polymer.nucleic and chain {chain} and {copy_selection}"
        count = residue_count(selection)
        if count <= 30:
            groups["mRNA"].append(chain)
        elif 50 <= count <= 90:
            groups["tRNA"].append(chain)
        else:
            groups["rRNA"].append(chain)
    return groups


def export_ribosome_components(object_name: str, pdb_id: str, out_dir: Path) -> list[dict]:
    chains = cmd.get_chains(object_name)
    copy_chains = first_ribosome_copy_chains(chains)
    copy_selection = chain_expression(copy_chains)
    exports = []

    protein_path = out_dir / f"{pdb_id}_surface_ribosomal_protein.obj"
    protein_report = export_selection(
        object_name,
        "ribosomal_protein",
        f"polymer.protein and {copy_selection}",
        protein_path,
    )
    protein_report["obj"] = str(protein_path.relative_to(ROOT)).replace("\\", "/")
    protein_report["copy_chains"] = copy_chains
    exports.append(protein_report)

    nucleic_groups = classify_ribosome_nucleic_chains(object_name, copy_selection)
    rna_path = out_dir / f"{pdb_id}_surface_rRNA.obj"
    rna_report = export_selection(
        object_name,
        "rRNA",
        f"polymer.nucleic and {chain_expression(nucleic_groups['rRNA'])}",
        rna_path,
    )
    rna_report["obj"] = str(rna_path.relative_to(ROOT)).replace("\\", "/")
    rna_report["chains"] = nucleic_groups["rRNA"]
    exports.append(rna_report)

    mrna_path = out_dir / f"{pdb_id}_surface_mRNA.obj"
    mrna_report = export_selection(
        object_name,
        "mRNA",
        f"polymer.nucleic and {chain_expression(nucleic_groups['mRNA'])}",
        mrna_path,
    )
    mrna_report["obj"] = str(mrna_path.relative_to(ROOT)).replace("\\", "/")
    mrna_report["chains"] = nucleic_groups["mRNA"]
    exports.append(mrna_report)

    for chain in nucleic_groups["tRNA"]:
        component = f"tRNA_{chain}"
        trna_path = out_dir / f"{pdb_id}_surface_{component}.obj"
        trna_report = export_selection(
            object_name,
            component,
            f"polymer.nucleic and chain {chain}",
            trna_path,
        )
        trna_report["obj"] = str(trna_path.relative_to(ROOT)).replace("\\", "/")
        trna_report["chains"] = [chain]
        exports.append(trna_report)
    return exports


def export_one(pdb_id: str) -> list[dict]:
    cif_path = RCSB_DIR / f"{pdb_id}.cif"
    if not cif_path.exists():
        raise FileNotFoundError(cif_path)

    object_name = f"{pdb_id}_surface"
    cmd.reinitialize()
    configure()
    cmd.load(safe_path(cif_path), object_name)
    cmd.remove(f"{object_name} and solvent")
    cmd.remove(f"{object_name} and hydro")

    out_dir = OUTPUT_DIR / pdb_id
    out_dir.mkdir(parents=True, exist_ok=True)
    if pdb_id == "4V5D":
        exports = export_ribosome_components(object_name, pdb_id, out_dir)
    else:
        exports = export_generic_components(object_name, pdb_id, out_dir)
    for entry in exports:
        entry["pdb_id"] = pdb_id
    return exports


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    requested = selected_pdbs()
    seen = set()
    exports = []

    for asset in manifest["pdb_assets"]:
        pdb_id = asset["pdb_id"].upper()
        if pdb_id in seen:
            continue
        seen.add(pdb_id)
        if requested is not None and pdb_id not in requested:
            continue
        print(f"Exporting {pdb_id}")
        exports.extend(export_one(pdb_id))

    report = {
        "title": "PyMOL surface asset exports",
        "surface_quality": 0,
        "ribosome_copy": os.environ.get("PYMOL_RIBOSOME_COPY", "first"),
        "exports": exports,
        "skipped": [],
    }
    if requested is None:
        report_path = OUTPUT_DIR / "surface_assets_manifest.json"
    else:
        subset = "_".join(sorted(requested))
        report_path = OUTPUT_DIR / f"surface_assets_manifest_subset_{subset}.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Wrote {report_path}")
    cmd.quit()


main()
