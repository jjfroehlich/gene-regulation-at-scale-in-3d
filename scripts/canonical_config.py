#!/usr/bin/env python3
"""Configuration for the canonical scene."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs" / "canonical"
MANIFEST_PATH = ROOT / "config" / "scene_manifest.json"
SCENE_BASENAME = "gene_expression_surface_style"
REPORT_KIND = "canonical_gene_end_transcription_layout"
PROTEIN_AA_CONTOUR_NM = 0.36

BLEND_PATH = OUTPUT_DIR / f"{SCENE_BASENAME}.blend"
PREVIEW_PATH = OUTPUT_DIR / f"preview_{SCENE_BASENAME}.png"
REPORT_PATH = OUTPUT_DIR / "gene_expression_surface_scene_report.json"

DNA_STRICT_CONTACT_NAMES = {
    "RNA polymerase II elongation complex",
    "Nucleosome",
    "Cas9",
    "Transcription factor 1",
    "Transcription factor 3",
    "Transcription factor 4",
    "p53 tetramer bound to DNA",
}
RNA_STRICT_CONTACT_NAMES = {
    "Pumilio RBP",
    "MS2 coat protein MCP",
    "Argonaute",
    "HuR-like RBP",
    "Poly(A)-binding RBP",
}
RNA_BRACKET_CONTACT_NAMES = {"Ribosome large subunit", "Ribosome small subunit"}
DISPLAY_EXEMPT_NAMES = {"mCherry/RFP tag", "Actin protein", "Standalone tRNA"}


def asset_by_name(manifest: dict, name: str) -> dict:
    for asset in manifest.get("pdb_assets", []):
        if asset.get("name") == name:
            return asset
    raise KeyError(f"Unknown PDB asset: {name}")


def write_manifest() -> dict:
    """Load and validate the single tracked canonical manifest.

    The function name is retained because the shared renderer calls it before a
    build. The manifest is now authoritative and is no longer derived from an
    older version.
    """
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(f"Missing canonical manifest: {MANIFEST_PATH}")
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    expected_outputs = {
        "blend": "outputs/canonical/gene_expression_surface_style.blend",
        "preview": "outputs/canonical/preview_gene_expression_surface_style.png",
        "report": "outputs/canonical/gene_expression_surface_scene_report.json",
    }
    if manifest.get("outputs") != expected_outputs:
        raise ValueError(f"Canonical manifest outputs do not match {expected_outputs}")
    return manifest


def main() -> None:
    manifest = write_manifest()
    print(json.dumps({"manifest": str(MANIFEST_PATH), "title": manifest["title"]}, indent=2))


if __name__ == "__main__":
    main()
