#!/usr/bin/env python3
"""Configuration for canonical V6 with transcription completed at the gene end."""

from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs" / "canonical"
V5_MANIFEST_PATH = ROOT / "config" / "scene_manifest_v5.json"
MANIFEST_PATH = ROOT / "config" / "scene_manifest_v6.json"
SCENE_BASENAME = "gene_expression_surface_style_v6"
REPORT_KIND = "canonical_v6_gene_end_transcription_layout"
PROTEIN_AA_CONTOUR_NM = 0.36

BLEND_PATH = OUTPUT_DIR / f"{SCENE_BASENAME}.blend"
PREVIEW_PATH = OUTPUT_DIR / f"preview_{SCENE_BASENAME}.png"
REPORT_PATH = OUTPUT_DIR / "gene_expression_surface_scene_v6_report.json"

DNA_END_MM = (40.73481801985251, -38.46183439342853, -0.000001123608390522432)
V5_MRNA_ORIGIN_MM = (28.117722241405147, 18.13054368897545, 0.1653815095123357)
BRANCH_TRANSLATION_MM = tuple(end - start for end, start in zip(DNA_END_MM, V5_MRNA_ORIGIN_MM))
DNA_LENGTH_MM = 537.744

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


def _translated(point: list[float] | tuple[float, ...]) -> list[float]:
    return [float(value) + delta for value, delta in zip(point, BRANCH_TRANSLATION_MM)]


def apply_canonical_v6(v5_manifest: dict) -> dict:
    manifest = copy.deepcopy(v5_manifest)
    manifest["canonical_version"] = "v6"
    manifest["title"] = "Scale-accurate gene expression scene - canonical V6 gene-end transcription layout"
    manifest["outputs"] = {
        "blend": "outputs/canonical/gene_expression_surface_style_v6.blend",
        "preview": "outputs/canonical/preview_gene_expression_surface_style_v6.png",
        "report": "outputs/canonical/gene_expression_surface_scene_v6_report.json",
    }

    mrna = manifest["procedural_nucleic_acids"]["mrna"]
    mrna["canonical_spiral"]["start_mm"] = list(DNA_END_MM)
    mrna["compact_center_mm"] = _translated(mrna["compact_center_mm"])
    mrna["elongated_path_origin"] = "3_prime_at_polymerase_ii"
    manifest["procedural_nucleic_acids"]["dna"]["path_mode"] = (
        "v6_reader_order_serpentine_with_nucleosome_loop"
    )

    polymerase = asset_by_name(manifest, "RNA polymerase II elongation complex")
    polymerase["location_mm"] = list(DNA_END_MM)
    polymerase["path_anchor"]["fraction"] = 1.0
    polymerase["path_anchor"]["distance_mm"] = DNA_LENGTH_MM
    polymerase["role"] = "RNA Pol II elongation complex at completed gene end"

    for asset in manifest.get("pdb_assets", []):
        if "v5_contact_role" in asset:
            asset["canonical_contact_role"] = asset.pop("v5_contact_role")
        if "v5_color_role" in asset:
            asset["canonical_color_role"] = asset.pop("v5_color_role")
        anchor = asset.get("path_anchor") or {}
        if anchor.get("path") == "mrna":
            asset["location_mm"] = _translated(asset["location_mm"])

    v5_mrna_end = (19.048343264999826, -30.34062402116573, 80.1653815095008)
    v6_mrna_end = _translated(v5_mrna_end)
    actin = asset_by_name(manifest, "Actin protein")
    actin["location_mm"] = [v6_mrna_end[0] - 6.4, v6_mrna_end[1] + 0.6, v6_mrna_end[2] + 27.0]
    actin.pop("arrangement_v5_positioning", None)
    actin.pop("canonical_positioning", None)

    layout = manifest.setdefault("layout_intent", {})
    layout.update(
        {
            "source_layout": "canonical_v5_rigid_branch_translation",
            "dna": "full ACTB promoter-plus-gene span in reader order; RNA Pol II is attached at the gene endpoint",
            "mrna": "full-length actin mRNA has its nascent 3-prime end at Pol II on the completed gene endpoint",
            "branch_transform": "rigid_translation_of_mrna_and_all_downstream_rna_to_protein_objects",
            "branch_translation_mm": list(BRANCH_TRANSLATION_MM),
            "v5_mrna_origin_mm": list(V5_MRNA_ORIGIN_MM),
            "v6_mrna_origin_mm": list(DNA_END_MM),
        }
    )
    return manifest


def write_manifest() -> dict:
    if not V5_MANIFEST_PATH.exists():
        raise FileNotFoundError(f"Missing canonical V5 manifest: {V5_MANIFEST_PATH}")
    manifest = apply_canonical_v6(json.loads(V5_MANIFEST_PATH.read_text(encoding="utf-8")))
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
