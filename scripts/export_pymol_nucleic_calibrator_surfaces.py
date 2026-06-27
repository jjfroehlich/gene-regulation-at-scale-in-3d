#!/usr/bin/env python3
"""Export protein-stripped PyMOL nucleic-acid surfaces for calibrator structures.

Run inside PyMOL:
    pymol -cq -d "run C:/path/to/scripts/export_pymol_nucleic_calibrator_surfaces.py"
"""

from __future__ import annotations

import json
import os
from pathlib import Path

from pymol import cmd


ROOT = Path(os.environ.get("GENE_SCENE_ROOT", Path.cwd())).resolve()
MANIFEST_PATH = ROOT / "config" / "scene_manifest.json"
RCSB_DIR = ROOT / "assets" / "rcsb"
OUTPUT_DIR = ROOT / "experiments" / "procedural_nucleic_acids" / "assets" / "pymol_calibrator_surfaces"
REPRESENTATIVE_CHAINS = {
    "9IOB": "B",
}


def safe_path(path: Path) -> str:
    return path.as_posix()


def selected_pdbs() -> set[str] | None:
    raw = os.environ.get("PYMOL_CALIBRATOR_PDBS", "").strip()
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


def residue_count(selection: str) -> int:
    seen = set()
    cmd.iterate(selection, "seen.add((chain, resi, resn))", space={"seen": seen})
    return len(seen)


def calibrator_assets(manifest: dict) -> list[dict]:
    assets = []
    for expected_type, group in manifest.get("nucleic_acid_calibrators", {}).items():
        for asset in group:
            item = dict(asset)
            item["expected_type"] = expected_type
            item["pdb_id"] = item["pdb_id"].upper()
            assets.append(item)
    return assets


def export_one(asset: dict) -> dict:
    pdb_id = asset["pdb_id"]
    cif_path = RCSB_DIR / f"{pdb_id}.cif"
    out_dir = OUTPUT_DIR / pdb_id
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{pdb_id}_nucleic_surface.obj"
    object_name = f"{pdb_id}_calibrator"

    if not cif_path.exists():
        return {"pdb_id": pdb_id, "expected_type": asset["expected_type"], "status": "missing_cif", "path": str(cif_path)}

    cmd.reinitialize()
    configure()
    cmd.load(safe_path(cif_path), object_name)
    cmd.remove(f"{object_name} and solvent")
    cmd.remove(f"{object_name} and hydro")
    representative_chain = None
    selection_body = "polymer.nucleic"
    if os.environ.get("PYMOL_CALIBRATOR_REPRESENTATIVE_CHAINS", "") == "1":
        representative_chain = REPRESENTATIVE_CHAINS.get(pdb_id)
        if representative_chain:
            selection_body = f"polymer.nucleic and chain {representative_chain}"
    cmd.select("calibrator_nucleic", f"{object_name} and {selection_body}")
    atom_count = cmd.count_atoms("calibrator_nucleic")
    residues = residue_count("calibrator_nucleic") if atom_count else 0
    if out_path.exists():
        out_path.unlink()
    if atom_count:
        cmd.hide("everything", object_name)
        cmd.show("surface", "calibrator_nucleic")
        cmd.save(safe_path(out_path), "calibrator_nucleic", state=1)
    cmd.delete("calibrator_nucleic")
    return {
        "pdb_id": pdb_id,
        "expected_type": asset["expected_type"],
        "status": "exported" if out_path.exists() and out_path.stat().st_size > 0 else "no_nucleic_surface",
        "atom_count": atom_count,
        "residue_count": residues,
        "selection": selection_body,
        "representative_chain": representative_chain,
        "obj": str(out_path.relative_to(ROOT)).replace("\\", "/"),
        "bytes": out_path.stat().st_size if out_path.exists() else 0,
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    requested = selected_pdbs()
    exports = []
    for asset in calibrator_assets(manifest):
        if requested is not None and asset["pdb_id"] not in requested:
            continue
        print(f"Exporting nucleic calibrator {asset['pdb_id']}")
        try:
            exports.append(export_one(asset))
        except Exception as exc:  # noqa: BLE001 - keep independent calibrators from blocking the rest.
            exports.append({"pdb_id": asset["pdb_id"], "expected_type": asset["expected_type"], "status": "error", "error": repr(exc)})
    report = {"title": "PyMOL nucleic-acid calibrator surfaces", "exports": exports}
    report_path = OUTPUT_DIR / "pymol_calibrator_surfaces_report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Wrote {report_path}")
    cmd.quit()


main()
