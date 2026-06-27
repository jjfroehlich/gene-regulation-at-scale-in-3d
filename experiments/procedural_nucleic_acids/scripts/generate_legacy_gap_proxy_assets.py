#!/usr/bin/env python3
"""Generate V2 gap-filled old-proxy DNA/RNA pseudoatom CIFs."""

from __future__ import annotations

import json
from pathlib import Path

import legacy_gap_proxy_geometry as gap


ROOT = Path(__file__).resolve().parents[3]
MANIFEST_PATH = ROOT / "config" / "scene_manifest.json"
ASSET_DIR = ROOT / "experiments" / "procedural_nucleic_acids" / "assets" / "legacy_gap_pymol_proxy"
CIF_DIR = ASSET_DIR / "cif"
REPORT_PATH = ASSET_DIR / "legacy_gap_proxy_generation_report.json"


def main() -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    dna = gap.build_dna_gap_asset(manifest)
    mrna = gap.build_mrna_gap_asset(manifest)
    gap.write_cif(CIF_DIR / f"{dna['asset_id']}.cif", dna["asset_id"], dna["atoms"])
    gap.write_cif(CIF_DIR / f"{mrna['asset_id']}.cif", mrna["asset_id"], mrna["atoms"])
    report = gap.write_generation_report(REPORT_PATH, manifest, dna, mrna)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
