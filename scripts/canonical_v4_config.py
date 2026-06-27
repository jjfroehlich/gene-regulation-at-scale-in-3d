#!/usr/bin/env python3
"""Canonical V4 manifest derivation from the arrangement V2 layout."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import arrangement_v1_config as v1
import arrangement_v2_config as v2


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "outputs" / "canonical"
MANIFEST_PATH = ROOT / "config" / "scene_manifest_v4.json"
SCENE_BASENAME = "gene_expression_surface_style_v4"
REPORT_KIND = "canonical_v4_from_arrangement_v2"

BLEND_PATH = OUTPUT_DIR / f"{SCENE_BASENAME}.blend"
PREVIEW_PATH = OUTPUT_DIR / f"preview_{SCENE_BASENAME}.png"
REPORT_PATH = OUTPUT_DIR / "gene_expression_surface_scene_v4_report.json"

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
    "Standalone tRNA",
}

RNA_BRACKET_CONTACT_NAMES = {
    "Ribosome large subunit",
    "Ribosome small subunit",
}

DISPLAY_EXEMPT_NAMES = {
    "mCherry/RFP tag",
    "Actin protein",
}


def asset_by_name(manifest: dict, name: str) -> dict:
    return v1.asset_by_name(manifest, name)


def mark_contact_role(asset: dict, role: str, strict: bool) -> None:
    asset["v4_contact_role"] = role
    asset["strict_contact_required"] = strict


def set_contact_anchor(
    asset: dict,
    path: str,
    fraction: float,
    roll_deg: float,
    binding_mode: str,
    *,
    distance_mm: float | None = None,
    local_offset_mm: list[float] | None = None,
) -> None:
    v1.set_path_anchor(
        asset,
        path,
        fraction,
        [0.0, 0.0, 0.0],
        local_offset_mm or [0.0, 0.0, 0.0],
        roll_deg,
        binding_mode,
        distance_mm,
    )


def set_v4_dna_contact_anchors(manifest: dict) -> None:
    anchors = {
        "Transcription factor 4": (0.045, 8.0, "co_crystal_dna_axis", None),
        "Cas9": (0.075, 24.0, "co_crystal_dna_axis", None),
        "Transcription factor 1": (0.105, -20.0, "co_crystal_dna_axis", None),
        "RNA polymerase II elongation complex": (v1.TSS_FRACTION, -18.0, "co_crystal_dna_axis", v1.TSS_DISTANCE_MM),
        "p53 tetramer bound to DNA": (0.205, 12.0, "co_crystal_dna_axis", None),
        "Transcription factor 3": (0.295, -5.0, "co_crystal_dna_axis", None),
        "Nucleosome": (0.42, 0.0, "wrapped_nucleosome_loop", None),
    }
    for name, (fraction, roll, binding_mode, distance_mm) in anchors.items():
        asset = asset_by_name(manifest, name)
        asset["anchor_kind"] = "nucleic"
        set_contact_anchor(asset, "dna", fraction, roll, binding_mode, distance_mm=distance_mm)
        asset["rotation_deg"] = [0.0, 0.0, roll]
        mark_contact_role(asset, "dna_strict_contact", True)
        if name == "Nucleosome":
            asset["protein_color"] = "tf_purple"


def set_v4_mrna_contact_anchors(manifest: dict) -> None:
    strict_anchors = {
        "Pumilio RBP": (0.20, -12.0, "co_crystal_rna_axis"),
        "MS2 coat protein MCP": (0.34, 16.0, "co_crystal_rna_axis"),
        "Argonaute": (0.54, 8.0, "protein_surface_rna_contact"),
        "HuR-like RBP": (0.64, 10.0, "co_crystal_rna_axis"),
        "Poly(A)-binding RBP": (0.76, -8.0, "co_crystal_rna_axis"),
        "Standalone tRNA": (0.915, -20.0, "co_crystal_rna_axis"),
    }
    for name, (fraction, roll, binding_mode) in strict_anchors.items():
        asset = asset_by_name(manifest, name)
        asset["anchor_kind"] = "nucleic"
        if binding_mode == "protein_surface_rna_contact":
            asset["anchor_kind"] = "protein_surface"
            asset["binding_side_mm"] = [0.0, 1.0, 0.0]
        elif name != "Standalone tRNA":
            asset["hide_alignment_nucleic"] = True
        set_contact_anchor(asset, "mrna", fraction, roll, binding_mode)
        asset["rotation_deg"] = [0.0, 0.0, roll]
        mark_contact_role(asset, "mrna_strict_contact", True)

    ribosome_anchors = {
        "Ribosome large subunit": (0.890, [0.0, 0.0, 0.0], 20.0, [0.0, -1.0, 0.0]),
        "Ribosome small subunit": (0.895, [0.0, 0.0, 0.0], 20.0, [0.0, 1.0, 0.0]),
    }
    for name, (fraction, local_offset, roll, binding_side) in ribosome_anchors.items():
        asset = asset_by_name(manifest, name)
        asset["anchor_kind"] = "nucleic"
        asset["binding_side_mm"] = binding_side
        set_contact_anchor(asset, "mrna", fraction, roll, "co_crystal_rna_bracket", local_offset_mm=local_offset)
        asset["rotation_deg"] = [74.0, 0.0, roll]
        mark_contact_role(asset, "mrna_bracket_contact", True)

    tag = asset_by_name(manifest, "mCherry/RFP tag")
    set_contact_anchor(tag, "mrna", 0.355, -18.0, "display_tag_near_mrna", local_offset_mm=[2.2, 1.2, 2.0])
    tag["rotation_deg"] = [0.0, 0.0, -18.0]
    mark_contact_role(tag, "display_exempt", False)


def mark_display_assets(manifest: dict) -> None:
    actin = asset_by_name(manifest, "Actin protein")
    mark_contact_role(actin, "product_exempt", False)
    actin["label_offset_mm"] = [0.0, 8.0, 0.4]


def apply_canonical_v4(base_manifest: dict) -> dict:
    manifest = copy.deepcopy(v2.apply_arrangement_v2(base_manifest))
    manifest["title"] = "Scale-accurate gene expression scene - canonical V4 ACTB full-gene contact layout"
    manifest["canonical_version"] = "v4"
    manifest["outputs"] = {
        "blend": "outputs/canonical/gene_expression_surface_style_v4.blend",
        "preview": "outputs/canonical/preview_gene_expression_surface_style_v4.png",
        "report": "outputs/canonical/gene_expression_surface_scene_v4_report.json",
    }
    manifest["layout_intent"] = {
        "source_layout": "experiments/arrangement_v2",
        "dna": "full ACTB canonical gene span plus 500 bp upstream promoter, folded into the arrangement V2 serpentine base",
        "mrna": "full-length actin mRNA starts at RNA polymerase II/TSS and rises as the arrangement V2 spiral",
        "contact_policy": "co-crystal nucleic guides are aligned to procedural DNA/RNA paths with near-zero guide offsets",
        "compact_reference": "separate compact full-length mRNP-like mRNA remains visible as a secondary reference",
    }

    mrna_proxy = manifest["procedural_nucleic_acids"]["mrna"]
    tss_point = v1.path_point(manifest, v1.TSS_DISTANCE_MM)
    mrna_start = [tss_point[0], tss_point[1], tss_point[2]]
    manifest["mrna"]["start_mm"] = mrna_start
    mrna_proxy["arrangement_v1"]["start_mm"] = mrna_start
    mrna_proxy["compact_center_mm"] = [78.0, 42.0, 18.0]

    set_v4_dna_contact_anchors(manifest)
    set_v4_mrna_contact_anchors(manifest)
    mark_display_assets(manifest)
    return manifest


def write_manifest(base_manifest_path: Path = ROOT / "config" / "scene_manifest.json") -> dict:
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    manifest = apply_canonical_v4(json.loads(base_manifest_path.read_text(encoding="utf-8")))
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
