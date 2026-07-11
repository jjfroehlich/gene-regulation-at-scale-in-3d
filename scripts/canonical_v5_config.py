#!/usr/bin/env python3
"""Standalone configuration access for the sole canonical V5 scene.

The checked-in V5 manifest is now the canonical, resolved source of truth.  Older
arrangement and canonical configuration modules are intentionally not imported.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs" / "canonical"
MANIFEST_PATH = ROOT / "config" / "scene_manifest_v5.json"
SCENE_BASENAME = "gene_expression_surface_style_v5"
REPORT_KIND = "canonical_v5_reader_order_layout"
PROTEIN_AA_CONTOUR_NM = 0.36

BLEND_PATH = OUTPUT_DIR / f"{SCENE_BASENAME}.blend"
PREVIEW_PATH = OUTPUT_DIR / f"preview_{SCENE_BASENAME}.png"
REPORT_PATH = OUTPUT_DIR / "gene_expression_surface_scene_v5_report.json"

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


def apply_canonical_v5(manifest: dict) -> dict:
    """Apply current-only defaults to an already resolved V5 manifest."""
    manifest["canonical_version"] = "v5"
    manifest["title"] = "Scale-accurate gene expression scene - canonical V5 reader-order layout"
    manifest["units"]["protein_aa_contour_nm"] = PROTEIN_AA_CONTOUR_NM
    manifest["units"]["protein_aa_to_mm"] = PROTEIN_AA_CONTOUR_NM * float(manifest["units"]["nm_to_mm"])
    manifest["outputs"] = {
        "blend": "outputs/canonical/gene_expression_surface_style_v5.blend",
        "preview": "outputs/canonical/preview_gene_expression_surface_style_v5.png",
        "report": "outputs/canonical/gene_expression_surface_scene_v5_report.json",
    }
    manifest.pop("legacy_outputs", None)
    mrna = manifest["procedural_nucleic_acids"]["mrna"]
    mrna.update(
        {
            "tube_radius_mm": 0.150,
            "tube_sides": 16,
            "surface_mode": "twisted_groove",
            "surface_bump_amplitude": 0.0,
            "base_lobe_every_nt": 1,
            "nucleotide_detail_every_nt": 1,
            "base_lobe_radius_mm": 0.095,
            "phosphate_radius_mm": 0.09,
            "sugar_radius_mm": 0.085,
            "rna_base_ellipsoid_radii_mm": [0.135, 0.078, 0.055],
            "base_connector_radius_mm": 0.048,
            "show_radial_nucleotide_detail": False,
            "elongated_path_origin": "3_prime_at_polymerase_ii",
            "direct_smooth_factor": 0.028,
            "direct_smooth_iterations": 1,
            "compact_variant": "compact_rosette",
            "stem_bridge_radius_mm": 0.032,
            "secondary_structure": {
                "model": "deterministic_schematic_stem_loop",
                "sequence_resolved": False,
                "elongated_stem_count": 0,
                "compact_stem_count": 38,
                "elongated_paired_fraction_target": 0.0,
                "compact_paired_fraction_target": 0.58,
                "stem_bp_min": 6,
                "stem_bp_max": 18,
                "hairpin_loop_nt_min": 4,
                "hairpin_loop_nt_max": 12,
                "a_form_bp_per_turn": 11.0,
                "a_form_diameter_mm": 0.92,
                "elongated_seed": 7103,
                "compact_seed": 3838
            },
        }
    )
    manifest.setdefault("layout_intent", {})["source_layout"] = "standalone_canonical_v5"
    manifest["layout_intent"]["contact_policy"] = "strict co-crystal/contact validation for all current canonical DNA/RNA binders"
    manifest["layout_intent"]["compact_reference"] = "RNA-only compact rosette with deterministic schematic paired stems"
    trna = asset_by_name(manifest, "Standalone tRNA")
    trna["path_anchor"]["fraction"] = 0.90
    trna["path_anchor"]["offset_local_mm"] = [0.0, 10.5, 0.0]
    trna["path_anchor"]["roll_deg"] = -52.0
    trna["path_anchor"]["binding_mode"] = "display_tRNA_side_view"
    trna["rotation_deg"] = [0.0, 0.0, -52.0]
    trna["display_plane_normal"] = [-13.0, -18.0, 15.0]
    trna["display_roll_deg"] = 18.0
    trna["label_offset_mm"] = [3.2, 2.8, 0.45]
    trna["strict_contact_required"] = False
    trna["contact_role"] = "display_exempt"
    trna["role"] = "incoming tRNA displayed beside translating ribosome"
    overview_labels = {
        "RNA polymerase II elongation complex": "RNA Pol II (2E2I)",
        "Ribosome small subunit": "30S subunit (1J5E)",
        "Ribosome large subunit": "large subunit (1JJ2)",
        "Standalone tRNA": "tRNA-Phe (4TNA)",
        "Nucleosome": "nucleosome (1AOI)",
        "Cas9": "Cas9-DNA (4UN3)",
        "Actin protein": "actin (1J6Z)",
        "Transcription factor 1": "ZBTB24 (6ML2)",
        "Transcription factor 4": "R2R3 MYB-DNA (6KKS)",
        "p53 tetramer bound to DNA": "p53-DNA (3TS8)",
        "Transcription factor 3": "FOXM1-DBD (3G73)",
        "Argonaute": "Argonaute (1U04)",
        "Poly(A)-binding RBP": "PABP-RNA (1CVJ)",
        "HuR-like RBP": "HuR-RNA (4ED5)",
        "Pumilio RBP": "PUM2-RNA (3Q0Q)",
        "MS2 coat protein MCP": "MS2-RNA (1ZDH)",
        "mCherry/RFP tag": "mCherry (2H5Q)",
    }
    for asset in manifest.get("pdb_assets", []):
        if asset.get("name") in overview_labels:
            asset["label_text"] = overview_labels[asset["name"]]
    ms2 = asset_by_name(manifest, "MS2 coat protein MCP")
    ms2["path_anchor"]["fraction"] = 0.29
    return manifest


def write_manifest(base_manifest_path: Path | None = None) -> dict:
    """Normalize the checked-in resolved V5 manifest in place.

    ``base_manifest_path`` is accepted for compatibility with the former writer;
    the versioned resolved manifest remains authoritative.
    """
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(f"Missing canonical manifest: {MANIFEST_PATH}")
    manifest = apply_canonical_v5(json.loads(MANIFEST_PATH.read_text(encoding="utf-8")))
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
