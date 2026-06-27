#!/usr/bin/env python3
"""Canonical V5 manifest derivation with reader-order DNA layout."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import arrangement_v1_config as v1
import canonical_v4_config as v4
import procedural_nucleic_geometry as geom


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs" / "canonical"
MANIFEST_PATH = ROOT / "config" / "scene_manifest_v5.json"
SCENE_BASENAME = "gene_expression_surface_style_v5"
REPORT_KIND = "canonical_v5_reader_order_layout"
PROTEIN_AA_CONTOUR_NM = 0.36

BLEND_PATH = OUTPUT_DIR / f"{SCENE_BASENAME}.blend"
PREVIEW_PATH = OUTPUT_DIR / f"preview_{SCENE_BASENAME}.png"
REPORT_PATH = OUTPUT_DIR / "gene_expression_surface_scene_v5_report.json"

DNA_STRICT_CONTACT_NAMES = set(v4.DNA_STRICT_CONTACT_NAMES)
RNA_STRICT_CONTACT_NAMES = set(v4.RNA_STRICT_CONTACT_NAMES)
RNA_BRACKET_CONTACT_NAMES = set(v4.RNA_BRACKET_CONTACT_NAMES)
DISPLAY_EXEMPT_NAMES = set(v4.DISPLAY_EXEMPT_NAMES)
LABEL_TEXT_BY_ASSET = {
    "RNA polymerase II elongation complex": "RNA Pol II elongation\ncomplex (2E2I)",
    "Ribosome small subunit": "30S ribosomal\nsubunit (1J5E)",
    "Ribosome large subunit": "large ribosomal\nsubunit (1JJ2)",
    "Standalone tRNA": "yeast tRNA-Phe\n(4TNA)",
    "Nucleosome": "nucleosome core +\n146 bp DNA (1AOI)",
    "Cas9": "Cas9-DNA\ncomplex (4UN3)",
    "Actin protein": "actin protein\nproduct (1J6Z)",
    "Transcription factor 1": "ZBTB24 zinc\nfingers + DNA (6ML2)",
    "Transcription factor 4": "R2R3 MYB TF +\nDNA (6KKS)",
    "p53 tetramer bound to DNA": "p53 tetramer +\nDNA (3TS8)",
    "Transcription factor 3": "FOXM1 DNA-binding\ndomain (3G73)",
    "Argonaute": "Argonaute\n(1U04)",
    "Poly(A)-binding RBP": "poly(A)-binding\nprotein + RNA (1CVJ)",
    "HuR-like RBP": "HuR RRM domains\n+ RNA (4ED5)",
    "Pumilio RBP": "Pumilio2 RNA-\nbinding domain (3Q0Q)",
    "MS2 coat protein MCP": "MS2 coat protein\n+ RNA (1ZDH)",
    "mCherry/RFP tag": "mCherry/RFP tag\n(2H5Q)",
}
NON_PRODUCT_PROTEIN_COLOR = "protein_uniform_slate"
PRODUCT_PROTEIN_COLOR = "protein_product_coral"
DNA_GUIDE_NUCLEIC_COLOR = "dna_guide_amber"
RNA_GUIDE_NUCLEIC_COLOR = "rna_guide_orange"


def asset_by_name(manifest: dict, name: str) -> dict:
    return v1.asset_by_name(manifest, name)


def set_v5_dna_contact_anchors(manifest: dict) -> None:
    anchors = {
        "Transcription factor 4": (0.030, 8.0, "co_crystal_dna_axis", None),
        "Cas9": (0.050, 24.0, "co_crystal_dna_axis", None),
        "Transcription factor 1": (0.070, -20.0, "co_crystal_dna_axis", None),
        "RNA polymerase II elongation complex": (v1.TSS_FRACTION, -18.0, "co_crystal_dna_axis", v1.TSS_DISTANCE_MM),
        "p53 tetramer bound to DNA": (0.095, 12.0, "co_crystal_dna_axis", None),
        "Nucleosome": (0.165, 0.0, "wrapped_nucleosome_loop", None),
        "Transcription factor 3": (0.205, -5.0, "co_crystal_dna_axis", None),
    }
    for name, (fraction, roll, binding_mode, distance_mm) in anchors.items():
        asset = asset_by_name(manifest, name)
        asset["anchor_kind"] = "nucleic"
        v4.set_contact_anchor(asset, "dna", fraction, roll, binding_mode, distance_mm=distance_mm)
        asset["rotation_deg"] = [0.0, 0.0, roll]
        v4.mark_contact_role(asset, "dna_strict_contact", True)
        asset["v5_contact_role"] = "dna_strict_contact"
        if name == "Nucleosome":
            asset["protein_color"] = "tf_purple"


def tighten_labels(manifest: dict) -> None:
    label_offsets = {
        "Transcription factor 4": [0.0, 3.0, 0.35],
        "Cas9": [0.0, 3.2, 0.35],
        "Transcription factor 1": [0.0, 3.0, 0.35],
        "RNA polymerase II elongation complex": [1.4, -3.2, 0.4],
        "p53 tetramer bound to DNA": [0.0, 3.2, 0.35],
        "Nucleosome": [-3.2, 3.0, 0.4],
        "Transcription factor 3": [0.0, 3.2, 0.35],
        "Pumilio RBP": [2.4, 2.0, 0.35],
        "MS2 coat protein MCP": [2.4, 2.0, 0.35],
        "Argonaute": [2.4, 2.0, 0.35],
        "HuR-like RBP": [2.4, 2.0, 0.35],
        "Poly(A)-binding RBP": [2.4, 2.0, 0.35],
        "Standalone tRNA": [2.4, 2.0, 0.35],
        "Ribosome large subunit": [3.0, 2.2, 0.45],
        "Ribosome small subunit": [3.0, 2.2, 0.45],
        "mCherry/RFP tag": [2.0, 1.8, 0.4],
        "Actin protein": [0.0, 3.6, 0.45],
    }
    for asset in manifest.get("pdb_assets", []):
        if asset["name"] in label_offsets:
            asset["label_offset_mm"] = label_offsets[asset["name"]]
            asset["label_size_mm"] = 1.55
        if asset["name"] in LABEL_TEXT_BY_ASSET:
            asset["label_text"] = LABEL_TEXT_BY_ASSET[asset["name"]]


def place_actin_above_translation_end(manifest: dict) -> None:
    mrna_model = geom.build_mrna_model(manifest)
    mrna_end = mrna_model["path"].points[-1]
    actin = asset_by_name(manifest, "Actin protein")
    actin.pop("path_anchor", None)
    offset = [-6.4, 0.6, 20.0]
    actin["location_mm"] = [mrna_end.x + offset[0], mrna_end.y + offset[1], mrna_end.z + offset[2]]
    actin["rotation_deg"] = [72.0, 0.0, 18.0]
    actin["arrangement_v5_positioning"] = {
        "mode": "above_translation_end",
        "mrna_end_mm": [mrna_end.x, mrna_end.y, mrna_end.z],
        "offset_mm": offset,
    }


def apply_v5_color_roles(manifest: dict) -> None:
    for asset in manifest.get("pdb_assets", []):
        asset["protein_color"] = PRODUCT_PROTEIN_COLOR if asset["name"] == "Actin protein" else NON_PRODUCT_PROTEIN_COLOR
        path = (asset.get("path_anchor") or {}).get("path")
        if path == "dna":
            asset["nucleic_color"] = DNA_GUIDE_NUCLEIC_COLOR
        elif path == "mrna":
            asset["nucleic_color"] = RNA_GUIDE_NUCLEIC_COLOR
        asset["v5_color_role"] = "highlighted_protein_product" if asset["name"] == "Actin protein" else "uniform_context_protein"


def apply_canonical_v5(base_manifest: dict) -> dict:
    manifest = copy.deepcopy(v4.apply_canonical_v4(base_manifest))
    manifest["title"] = "Scale-accurate gene expression scene - canonical V5 reader-order layout"
    manifest["canonical_version"] = "v5"
    manifest["units"]["protein_aa_contour_nm"] = PROTEIN_AA_CONTOUR_NM
    manifest["units"]["protein_aa_to_mm"] = PROTEIN_AA_CONTOUR_NM * float(manifest["units"]["nm_to_mm"])
    manifest["outputs"] = {
        "blend": "outputs/canonical/gene_expression_surface_style_v5.blend",
        "preview": "outputs/canonical/preview_gene_expression_surface_style_v5.png",
        "report": "outputs/canonical/gene_expression_surface_scene_v5_report.json",
    }
    manifest["layout_intent"] = {
        "source_layout": "canonical_v4_contact_layout",
        "dna": "full ACTB promoter-plus-gene span in reader order from top-left toward bottom-right",
        "mrna": "full-length actin mRNA starts at the V5 Pol II/TSS point and keeps the elongated spiral",
        "contact_policy": "V4 strict co-crystal/contact validation retained; V5 reclusters DNA binders near the gene start",
        "compact_reference": "compact full-length mRNP-like mRNA is dynamically placed between the elongated mRNA end and actin",
        "wire_policy": "editable source curves are kept for rebaking but hidden in viewport and render",
    }

    dna_proxy = manifest["procedural_nucleic_acids"]["dna"]
    dna_proxy["path_mode"] = "v5_reader_order_serpentine_with_nucleosome_loop"
    dna_proxy["full_gene_serpentine"].update(
        {
            "target_bp": v1.ACTB_WITH_PROMOTER_BP,
            "width_mm": 112.0,
            "height_mm": 74.0,
            "lanes": 5,
            "samples_per_lane": 16,
            "nucleosome_loop_fraction": 0.165,
            "nucleosome_loop_roll_offset_deg": 0.0,
            "reader_order": "promoter_top_left_to_gene_end_bottom_right",
        }
    )

    tss_point = v1.path_point(manifest, v1.TSS_DISTANCE_MM)
    mrna_start = [tss_point[0], tss_point[1], tss_point[2]]
    manifest["mrna"]["start_mm"] = mrna_start
    mrna_proxy = manifest["procedural_nucleic_acids"]["mrna"]
    mrna_proxy["arrangement_v1"]["start_mm"] = mrna_start
    mrna_proxy["arrangement_v1"]["terminal_drift_mm"] = [6.55, -51.35, 0.0]
    mrna_proxy["compact_position_mode"] = "between_mrna_end_and_actin"
    mrna_proxy["compact_center_mm"] = [15.5, -30.0, 91.0]

    place_actin_above_translation_end(manifest)
    set_v5_dna_contact_anchors(manifest)
    v4.set_v4_mrna_contact_anchors(manifest)
    v4.mark_display_assets(manifest)
    apply_v5_color_roles(manifest)
    tighten_labels(manifest)
    return manifest


def write_manifest(base_manifest_path: Path = ROOT / "config" / "scene_manifest.json") -> dict:
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    manifest = apply_canonical_v5(json.loads(base_manifest_path.read_text(encoding="utf-8")))
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
