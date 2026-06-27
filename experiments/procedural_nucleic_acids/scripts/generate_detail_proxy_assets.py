#!/usr/bin/env python3
"""Generate experiment-local detailed pseudoatom CIFs for PyMOL surface export."""

from __future__ import annotations

import json
from pathlib import Path

import detail_route_geometry as geom


ROOT = Path(__file__).resolve().parents[3]
EXPERIMENT_DIR = ROOT / "experiments" / "procedural_nucleic_acids"
ASSET_DIR = EXPERIMENT_DIR / "assets" / "detail_route_pymol_proxy"
CIF_DIR = ASSET_DIR / "cif"
REPORT_PATH = ASSET_DIR / "detail_proxy_generation_report.json"


def main() -> None:
    dna = geom.build_dna_detail_atoms()
    rna = geom.build_rna_detail_atoms()
    geom.write_cif(CIF_DIR / f"{dna['asset_id']}.cif", dna["asset_id"], dna["atoms"])
    geom.write_cif(CIF_DIR / f"{rna['asset_id']}.cif", rna["asset_id"], rna["atoms"])
    report = geom.write_generation_report(REPORT_PATH, dna, rna)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
