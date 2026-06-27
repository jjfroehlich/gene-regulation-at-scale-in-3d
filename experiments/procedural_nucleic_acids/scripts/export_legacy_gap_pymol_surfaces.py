#!/usr/bin/env python3
"""Export V2 gap-filled old-proxy DNA/RNA surfaces through PyMOL."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

from pymol import cmd


ROOT = Path(os.environ.get("GENE_SCENE_ROOT", Path.cwd())).resolve()
SCRIPT_DIR = ROOT / "experiments" / "procedural_nucleic_acids" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))
import legacy_gap_proxy_geometry as gap  # noqa: E402


ASSET_DIR = ROOT / "experiments" / "procedural_nucleic_acids" / "assets" / "legacy_gap_pymol_proxy"
CIF_DIR = ASSET_DIR / "cif"
RAW_DIR = ASSET_DIR / "raw_surfaces"
REPORT_PATH = RAW_DIR / "legacy_gap_pymol_surface_exports_report.json"


def safe_path(path: Path) -> str:
    return path.as_posix()


def configure() -> None:
    cmd.set("surface_quality", 0)
    cmd.set("surface_smooth_edges", 1)
    cmd.set("solvent_radius", 0.0)
    cmd.set("specular", 0)
    cmd.set("ambient", 0.72)
    cmd.set("direct", 0.32)
    cmd.set("shininess", 0)
    cmd.set("two_sided_lighting", 1)


def export_component(object_name: str, asset_id: str, component: str, selection_body: str, vdw_A: float) -> dict:
    out_dir = RAW_DIR / asset_id
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{asset_id}_surface_{component}.obj"
    if out_path.exists():
        out_path.unlink()
    selection = f"{asset_id}_{component}_sel"
    cmd.select(selection, f"{object_name} and ({selection_body})")
    atom_count = cmd.count_atoms(selection)
    if atom_count:
        cmd.alter(selection, f"vdw={float(vdw_A)}")
        cmd.rebuild()
        cmd.hide("everything", object_name)
        cmd.show("surface", selection)
        cmd.save(safe_path(out_path), selection, state=1)
    cmd.delete(selection)
    return {
        "asset_id": asset_id,
        "component": component,
        "selection": selection_body,
        "surface_vdw_A": vdw_A,
        "atom_count": atom_count,
        "obj": str(out_path.relative_to(ROOT)).replace("\\", "/"),
        "exists": out_path.exists(),
        "bytes": out_path.stat().st_size if out_path.exists() else 0,
    }


def export_asset(asset: dict) -> list[dict]:
    asset_id = asset["asset_id"]
    cif_path = CIF_DIR / f"{asset_id}.cif"
    if not cif_path.exists():
        raise FileNotFoundError(cif_path)
    object_name = f"{asset_id}_surface"
    cmd.reinitialize()
    configure()
    cmd.load(safe_path(cif_path), object_name)
    exports = []
    for component, selection in asset["components"].items():
        exports.append(export_component(object_name, asset_id, component, selection, float(asset["surface_vdw_A"][component])))
    return exports


def main() -> None:
    manifest = json.loads((ROOT / "config" / "scene_manifest.json").read_text(encoding="utf-8"))
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    exports = []
    for asset in (gap.build_dna_gap_asset(manifest), gap.build_mrna_gap_asset(manifest)):
        print(f"Exporting V2 legacy gap-filled PyMOL surface {asset['asset_id']}")
        exports.extend(export_asset(asset))
    report = {
        "title": "V2 legacy PyMOL gap-filled surface exports",
        "coordinate_units": "angstrom",
        "angstrom_to_mm": gap.ANGSTROM_TO_MM,
        "surface_quality": 0,
        "solvent_radius": 0.0,
        "exports": exports,
    }
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Wrote {REPORT_PATH}")
    cmd.quit()


main()
