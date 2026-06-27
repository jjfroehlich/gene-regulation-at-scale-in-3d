#!/usr/bin/env python3
"""Export PyMOL surfaces for generated procedural DNA/RNA pseudoatom PDBs.

Run inside PyMOL:
    pymol -cq -d "run C:/path/to/scripts/export_pymol_procedural_nucleic_surfaces.py"
"""

from __future__ import annotations

import json
import os
from pathlib import Path

from pymol import cmd


ROOT = Path(os.environ.get("GENE_SCENE_ROOT", Path.cwd())).resolve()
MANIFEST_PATH = Path(os.environ.get("GENE_SCENE_MANIFEST", ROOT / "config" / "scene_manifest.json"))
CIF_DIR = ROOT / "assets" / "procedural_nucleic_acids" / "cif"
OUTPUT_DIR = ROOT / "assets" / "pymol_exports" / "surface_assets"
VARIANT = os.environ.get("PROCEDURAL_ASSET_VARIANT", "").strip().upper()


def asset_id(base: str) -> str:
    return f"{base}_{VARIANT}" if VARIANT else base


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


def export_component(object_name: str, asset_id: str, component: str, selection_body: str, out_dir: Path) -> dict:
    out_path = out_dir / f"{asset_id}_surface_{component}.obj"
    if out_path.exists():
        out_path.unlink()
    selection = f"{asset_id}_{component}_sel"
    cmd.select(selection, f"{object_name} and ({selection_body})")
    atom_count = cmd.count_atoms(selection)
    if atom_count:
        cmd.hide("everything", object_name)
        cmd.show("surface", selection)
        cmd.save(safe_path(out_path), selection, state=1)
    cmd.delete(selection)
    return {
        "asset_id": asset_id,
        "component": component,
        "selection": selection_body,
        "atom_count": atom_count,
        "obj": str(out_path.relative_to(ROOT)).replace("\\", "/"),
        "exists": out_path.exists(),
        "bytes": out_path.stat().st_size if out_path.exists() else 0,
    }


def export_asset(asset_id: str, cif_path: Path, components: dict[str, str], surface_vdw_A: float) -> list[dict]:
    if not cif_path.exists():
        raise FileNotFoundError(cif_path)
    object_name = f"{asset_id}_surface"
    cmd.reinitialize()
    configure()
    cmd.load(safe_path(cif_path), object_name)
    cmd.alter(object_name, f"vdw={float(surface_vdw_A)}")
    cmd.rebuild()
    out_dir = OUTPUT_DIR / asset_id
    out_dir.mkdir(parents=True, exist_ok=True)
    exports = []
    for component, selection in components.items():
        exports.append(export_component(object_name, asset_id, component, selection, out_dir))
    return exports


def main() -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    proxy = manifest.get("procedural_nucleic_surfaces", manifest.get("procedural_nucleic_acids", {}))
    assets = [
        {
            "asset_id": asset_id("DNA_PROXY"),
            "cif": CIF_DIR / f"{asset_id('DNA_PROXY')}.cif",
            "components": {"strand_A": "chain A", "strand_B": "chain B", "base_pairs": "chain C"},
            "surface_vdw_A": proxy.get("dna", {}).get("surface_vdw_A", 2.4),
        },
        {
            "asset_id": asset_id("MRNA_PROXY"),
            "cif": CIF_DIR / f"{asset_id('MRNA_PROXY')}.cif",
            "components": {"utr5": "chain A", "coding": "chain B", "utr3": "chain C"},
            "surface_vdw_A": proxy.get("mrna", {}).get("surface_vdw_A", 2.2),
        },
        {
            "asset_id": asset_id("MRNA_COMPACT_PROXY"),
            "cif": CIF_DIR / f"{asset_id('MRNA_COMPACT_PROXY')}.cif",
            "components": {"utr5": "chain A", "coding": "chain B", "utr3": "chain C"},
            "surface_vdw_A": proxy.get("mrna", {}).get("surface_vdw_A", 2.2),
        },
    ]
    requested = {
        value.strip().upper()
        for value in os.environ.get("PROCEDURAL_SURFACE_IDS", "").split(",")
        if value.strip()
    }
    if requested:
        assets = [asset for asset in assets if asset["asset_id"].upper() in requested]
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    exports = []
    for asset in assets:
        print(f"Exporting procedural surface {asset['asset_id']}")
        exports.extend(export_asset(asset["asset_id"], asset["cif"], asset["components"], asset["surface_vdw_A"]))
    report = {
        "title": "Procedural PyMOL surface proxy exports",
        "exports": exports,
        "coordinate_units": "angstrom",
        "surface_quality": 0,
        "solvent_radius": 0.0,
    }
    report_path = OUTPUT_DIR / "procedural_nucleic_surfaces_manifest.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Wrote {report_path}")
    cmd.quit()


main()
