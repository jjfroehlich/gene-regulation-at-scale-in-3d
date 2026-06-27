#!/usr/bin/env python3
"""Configuration for the arrangement V2 protein-placement experiment."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import arrangement_v1_config as v1
import procedural_nucleic_geometry as geom


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT_DIR = ROOT / "experiments" / "arrangement_v2"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
MANIFEST_PATH = EXPERIMENT_DIR / "arrangement_v2_manifest.json"
EXPERIMENT_KEY = "arrangement_v2"
SCENE_BASENAME = "gene_expression_arrangement_v2"
REPORT_KIND = "arrangement_v2_protein_arrangement_experiment"
LAYOUT_INTENT_OVERRIDES = {
    "translation": "RNA-bound proteins follow the RNA spiral, ribosome is near the upper spiral, and actin is centered above the RNA",
    "protein_arrangement": "DNA-binding proteins bind the DNA path; RNA-binding proteins bind the RNA path; actin is deliberately placed as the central product above translation",
}


def asset_by_name(manifest: dict, name: str) -> dict:
    return v1.asset_by_name(manifest, name)


def set_path_anchor(*args, **kwargs) -> None:
    v1.set_path_anchor(*args, **kwargs)


def set_v2_dna_anchors(manifest: dict) -> None:
    anchors = {
        "Transcription factor 4": (0.045, [0.0, 0.0, 0.0], [0.0, -4.5, 1.4], 8.0, "co_crystal_dna_axis", None),
        "Cas9": (0.075, [0.0, 0.0, 0.0], [0.0, 5.0, 1.2], 24.0, "co_crystal_dna_axis", None),
        "Transcription factor 1": (0.105, [0.0, 0.0, 0.0], [0.0, -4.2, 1.1], -20.0, "co_crystal_dna_axis", None),
        "RNA polymerase II elongation complex": (v1.TSS_FRACTION, [0.0, 0.0, 0.0], [0.0, 2.2, 0.8], -18.0, "co_crystal_dna_axis", v1.TSS_DISTANCE_MM),
        "p53 tetramer bound to DNA": (0.205, [0.0, 0.0, 0.0], [0.0, 5.4, 1.3], 12.0, "co_crystal_dna_axis", None),
        "Transcription factor 3": (0.295, [0.0, 0.0, 0.0], [0.0, -5.0, 1.0], -5.0, "co_crystal_dna_axis", None),
        "Nucleosome": (0.42, [0.0, 0.0, 0.0], [0.0, 4.0, 0.0], 0.0, "wrapped_nucleosome_loop", None),
    }
    for name, (fraction, offset, local_offset, roll, binding_mode, distance_mm) in anchors.items():
        asset = asset_by_name(manifest, name)
        set_path_anchor(asset, "dna", fraction, offset, local_offset, roll, binding_mode, distance_mm)
        asset["rotation_deg"] = [0.0, 0.0, roll]
        if name == "Nucleosome":
            asset["protein_color"] = "tf_purple"


def set_v2_mrna_anchors(manifest: dict) -> None:
    anchors = {
        "Pumilio RBP": (0.20, [0.0, 0.0, 0.0], [3.0, 1.8, 1.2], -12.0, "rna_tangent_frame"),
        "MS2 coat protein MCP": (0.34, [0.0, 0.0, 0.0], [-3.4, 1.9, 2.8], 16.0, "rna_tangent_frame"),
        "mCherry/RFP tag": (0.355, [0.0, 0.0, 0.0], [6.0, 3.8, 5.6], -18.0, "manual_tangent_frame"),
        "Argonaute": (0.54, [0.0, 0.0, 0.0], [-5.2, 1.7, 2.0], 8.0, "rna_tangent_frame"),
        "HuR-like RBP": (0.64, [0.0, 0.0, 0.0], [5.2, 1.8, 2.0], 10.0, "rna_tangent_frame"),
        "Poly(A)-binding RBP": (0.76, [0.0, 0.0, 0.0], [-3.2, 1.6, 4.0], -8.0, "rna_tangent_frame"),
        "Ribosome large subunit": (0.89, [0.0, 0.0, 0.0], [-1.4, 3.8, -4.0], 20.0, "rna_tangent_frame"),
        "Ribosome small subunit": (0.895, [0.0, 0.0, 0.0], [2.0, 6.6, 4.8], 20.0, "rna_tangent_frame"),
        "Standalone tRNA": (0.915, [0.0, 0.0, 0.0], [8.0, 5.5, 0.8], -20.0, "rna_tangent_frame"),
    }
    for name, (fraction, offset, local_offset, roll, binding_mode) in anchors.items():
        asset = asset_by_name(manifest, name)
        set_path_anchor(asset, "mrna", fraction, offset, local_offset, roll, binding_mode)
        if name.startswith("Ribosome"):
            asset["rotation_deg"] = [74.0, 0.0, 20.0]
        else:
            asset["rotation_deg"] = [0.0, 0.0, roll]


def center_actin_above_rna(manifest: dict) -> None:
    mrna_model = geom.build_mrna_model(manifest)
    mrna_path = mrna_model["path"]
    top = max(mrna_path.points, key=lambda point: point.z)
    min_x = min(point.x for point in mrna_path.points)
    max_x = max(point.x for point in mrna_path.points)
    min_y = min(point.y for point in mrna_path.points)
    max_y = max(point.y for point in mrna_path.points)
    center_x = (min_x + max_x) * 0.5
    center_y = (min_y + max_y) * 0.5
    actin = asset_by_name(manifest, "Actin protein")
    actin.pop("path_anchor", None)
    actin["location_mm"] = [center_x, center_y, top.z + 18.0]
    actin["rotation_deg"] = [72.0, 0.0, 18.0]
    actin["label_offset_mm"] = [0.0, 8.0, 0.4]
    actin["arrangement_v2_positioning"] = {
        "mode": "centered_above_rna_spiral",
        "rna_top_mm": [top.x, top.y, top.z],
        "center_xy_mm": [center_x, center_y],
        "vertical_offset_mm": 18.0,
    }


def apply_arrangement_v2(base_manifest: dict) -> dict:
    manifest = copy.deepcopy(v1.apply_arrangement_v1(base_manifest))
    manifest["title"] = manifest["title"].replace("arrangement experiment V1", "arrangement experiment V2 protein placement")
    manifest["outputs"] = {
        "blend": "experiments/arrangement_v2/outputs/gene_expression_arrangement_v2.blend",
        "preview": "experiments/arrangement_v2/outputs/preview_gene_expression_arrangement_v2.png",
        "report": "experiments/arrangement_v2/outputs/gene_expression_arrangement_v2_report.json",
    }
    set_v2_dna_anchors(manifest)
    set_v2_mrna_anchors(manifest)
    center_actin_above_rna(manifest)
    return manifest


def write_manifest(base_manifest_path: Path = ROOT / "config" / "scene_manifest.json") -> dict:
    EXPERIMENT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = apply_arrangement_v2(json.loads(base_manifest_path.read_text(encoding="utf-8")))
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
