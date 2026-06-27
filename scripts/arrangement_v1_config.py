#!/usr/bin/env python3
"""Configuration helpers for the polymerase-origin spiral RNA arrangement experiment."""

from __future__ import annotations

import copy
import json
import os
from pathlib import Path

import procedural_nucleic_geometry as geom


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT_DIR = ROOT / "experiments" / "arrangement_v1"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
MANIFEST_PATH = EXPERIMENT_DIR / "arrangement_v1_manifest.json"
EDITED_PATHS_PATH = EXPERIMENT_DIR / "edited_source_paths.json"
EXPERIMENT_KEY = "arrangement_v1"
SCENE_BASENAME = "gene_expression_arrangement_v1"
REPORT_KIND = "arrangement_v1_experiment"
ASSET_VARIANT = "ARRANGEMENT_V1"
DNA_ASSET_ID = f"DNA_PROXY_{ASSET_VARIANT}"
MRNA_ASSET_ID = f"MRNA_PROXY_{ASSET_VARIANT}"
ACTB_CANONICAL_GENE_BP = 3454
UPSTREAM_PROMOTER_BP = 500
ACTB_WITH_PROMOTER_BP = ACTB_CANONICAL_GENE_BP + UPSTREAM_PROMOTER_BP
TSS_DISTANCE_MM = UPSTREAM_PROMOTER_BP * 0.136
TSS_FRACTION = UPSTREAM_PROMOTER_BP / ACTB_WITH_PROMOTER_BP


def actb_gene_annotation() -> dict:
    segments = [{"name": "ACTB upstream promoter", "kind": "promoter", "length_bp": 500, "material": "dna_promoter_blue"}]
    for index, length in enumerate([78, 129, 240, 439, 182, 744], start=1):
        segments.append({"name": f"ACTB exon {index}", "kind": "exon", "length_bp": length, "material": "dna_exon_orange"})
        introns = [860, 134, 441, 95, 112]
        if index <= len(introns):
            segments.append({"name": f"ACTB intron {index}", "kind": "intron", "length_bp": introns[index - 1], "material": "dna_intron_olive"})
    return {
        "basis": "human ACTB canonical 3454 bp genomic transcript span plus 500 bp upstream promoter",
        "actb_canonical_gene_bp": ACTB_CANONICAL_GENE_BP,
        "upstream_promoter_bp": UPSTREAM_PROMOTER_BP,
        "total_bp": ACTB_WITH_PROMOTER_BP,
        "tss_bp": UPSTREAM_PROMOTER_BP,
        "segments": segments,
    }


def asset_by_name(manifest: dict, name: str) -> dict:
    for asset in manifest["pdb_assets"]:
        if asset["name"] == name:
            return asset
    raise KeyError(name)


def set_path_anchor(
    asset: dict,
    path: str,
    fraction: float,
    offset_mm: list[float] | None = None,
    offset_local_mm: list[float] | None = None,
    roll_deg: float = 0.0,
    binding_mode: str = "manual_tangent_frame",
    distance_mm: float | None = None,
) -> None:
    asset["path_anchor"] = {
        "path": path,
        "fraction": fraction,
        "offset_mm": offset_mm or [0.0, 0.0, 0.0],
        "offset_local_mm": offset_local_mm or [0.0, 0.0, 0.0],
        "roll_deg": roll_deg,
        "binding_mode": binding_mode,
    }
    if distance_mm is not None:
        asset["path_anchor"]["distance_mm"] = distance_mm


def path_point(manifest: dict, distance_mm: float) -> tuple[float, float, float]:
    path = geom.SampledPath(geom.catmull_rom(geom.dna_controls(manifest), 32))
    point = path.point_at_length(distance_mm)
    return point.to_tuple()


def apply_arrangement_v1(base_manifest: dict) -> dict:
    manifest = copy.deepcopy(base_manifest)
    manifest["title"] = f"{manifest['title']} - arrangement experiment V1"
    manifest["outputs"] = {
        "blend": f"experiments/{EXPERIMENT_KEY}/outputs/{SCENE_BASENAME}.blend",
        "preview": f"experiments/{EXPERIMENT_KEY}/outputs/preview_{SCENE_BASENAME}.png",
        "report": f"experiments/{EXPERIMENT_KEY}/outputs/{SCENE_BASENAME}_report.json",
    }

    dna_proxy = manifest["procedural_nucleic_acids"]["dna"]
    dna_proxy.update(
        {
            "style": "direct_blender_surface_proxy_polished",
            "strand_sides": 12,
            "strand_radius_mm": 0.12,
            "base_pair_every_bp": 1,
            "base_pair_radius_mm": 0.082,
            "base_pair_sides": 10,
            "base_stack_radius_mm": 0.172,
            "base_stack_sides": 12,
            "surface_bump_amplitude": 0.05,
            "direct_voxel_size_mm": 0.052,
            "direct_smooth_factor": 0.035,
            "direct_smooth_iterations": 1,
            "direct_unified_mesh": True,
        }
    )
    dna_proxy["path_mode"] = "full_gene_serpentine_with_nucleosome_loop"
    dna_proxy["nucleosome_loop"]["enabled"] = True
    dna_proxy["nucleosome_loop"]["source_pdb_id"] = "1AOI"
    dna_proxy["gene_annotation"] = actb_gene_annotation()
    dna_proxy["full_gene_serpentine"] = {
        "center_mm": [0.0, -9.0, 0.0],
        "target_bp": ACTB_WITH_PROMOTER_BP,
        "width_mm": 112.0,
        "height_mm": 74.0,
        "lanes": 6,
        "samples_per_lane": 16,
        "irregularity_mm": 3.2,
        "z_irregularity_mm": 0.42,
        "phase": 0.53,
        "length_variation": 0.23,
        "nucleosome_loop_fraction": 0.42,
        "nucleosome_loop_roll_offset_deg": 0.0,
        "nucleosome_entry_bridge_mm": 8.0,
        "nucleosome_exit_bridge_mm": 8.0,
    }

    tss_point = path_point(manifest, TSS_DISTANCE_MM)
    mrna_start = [tss_point[0], tss_point[1], tss_point[2] + 2.0]
    manifest["mrna"]["start_mm"] = mrna_start
    mrna_proxy = manifest["procedural_nucleic_acids"]["mrna"]
    mrna_proxy.update(
        {
            "style": "direct_blender_surface_proxy_polished",
            "tube_sides": 12,
            "tube_radius_mm": 0.165,
            "surface_bump_amplitude": 0.1,
            "base_lobe_every_nt": 1,
            "base_lobe_radius_mm": 0.095,
            "base_lobe_sides": 8,
            "base_offset_mm": 0.17,
            "nucleotide_detail_every_nt": 1,
            "phosphate_radius_mm": 0.09,
            "sugar_radius_mm": 0.085,
            "rna_base_ellipsoid_radii_mm": [0.135, 0.078, 0.055],
            "base_connector_radius_mm": 0.048,
            "nucleotide_detail_rings": 5,
            "nucleotide_detail_sides": 8,
            "direct_voxel_size_mm": 0.05,
            "direct_smooth_factor": 0.028,
            "direct_smooth_iterations": 1,
            "direct_unified_mesh": True,
        }
    )
    mrna_proxy["path_mode"] = "arrangement_v1_polymerase_spiral"
    mrna_proxy["arrangement_v1"] = {
        "start_mm": mrna_start,
        "start_radius_mm": 22.0,
        "end_radius_mm": 7.0,
        "vertical_rise_mm": 80.0,
        "rise_axis": "z",
        "start_angle_deg": 0.0,
        "orientation": 1.0,
        "samples": 2400,
        "min_turns": 0.8,
        "max_turns": 4.5,
        "radial_irregularity_mm": 2.2,
        "vertical_irregularity_mm": 2.0,
        "lateral_irregularity_mm": 3.6,
        "z_irregularity_mm": 1.0,
    }

    if os.environ.get("ARRANGEMENT_V1_USE_EDITED_PATHS") == "1" and EDITED_PATHS_PATH.exists():
        edited = json.loads(EDITED_PATHS_PATH.read_text(encoding="utf-8"))
        dna_points = edited.get("DNA_source_path", {}).get("points_mm", [])
        mrna_points = edited.get("MRNA_source_path", {}).get("points_mm", [])
        if len(dna_points) >= 2:
            dna_proxy["path_mode"] = "custom_dna_controls"
            dna_proxy["custom_controls_mm"] = dna_points
            if edited.get("nucleosome_loop"):
                dna_proxy["nucleosome_loop"].update(edited["nucleosome_loop"])
        if len(mrna_points) >= 2:
            mrna_proxy["path_mode"] = "custom_mrna_path"
            mrna_proxy["custom_points_mm"] = mrna_points

    dna_anchor_fractions = {
        "Cas9": (0.07, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 24.0, "co_crystal_dna_axis", None),
        "Transcription factor 4": (0.04, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 8.0, "co_crystal_dna_axis", None),
        "Transcription factor 1": (0.11, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], -20.0, "co_crystal_dna_axis", None),
        "p53 tetramer bound to DNA": (0.18, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 12.0, "co_crystal_dna_axis", None),
        "Transcription factor 3": (0.26, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], -5.0, "co_crystal_dna_axis", None),
        "Nucleosome": (0.42, [0.0, 0.0, 0.0], [0.0, 4.0, 0.0], 0.0, "wrapped_nucleosome_loop", None),
        "RNA polymerase II elongation complex": (TSS_FRACTION, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], -18.0, "co_crystal_dna_axis", TSS_DISTANCE_MM),
    }
    for name, (fraction, offset, local_offset, roll, binding_mode, distance_mm) in dna_anchor_fractions.items():
        asset = asset_by_name(manifest, name)
        set_path_anchor(asset, "dna", fraction, offset, local_offset, roll, binding_mode, distance_mm)
        asset["rotation_deg"] = [0.0, 0.0, roll]
        if name == "RNA polymerase II elongation complex":
            asset["label_offset_mm"] = [1.5, -13.0, 0.4]
        elif name == "Nucleosome":
            asset["protein_color"] = "tf_purple"
            asset["label_offset_mm"] = [-12.0, 8.0, 0.4]

    mrna_anchor_fractions = {
        "Pumilio RBP": (0.22, [0.0, 0.0, 0.0], [2.0, 1.2, 0.0], -10.0, "rna_tangent_frame"),
        "MS2 coat protein MCP": (0.36, [0.0, 0.0, 0.0], [-2.0, 1.4, 3.0], 16.0, "rna_tangent_frame"),
        "mCherry/RFP tag": (0.37, [0.0, 0.0, 0.0], [5.0, 3.0, 6.0], -18.0, "manual_tangent_frame"),
        "Argonaute": (0.58, [0.0, 0.0, 0.0], [-4.0, 1.4, 2.0], 8.0, "rna_tangent_frame"),
        "HuR-like RBP": (0.67, [0.0, 0.0, 0.0], [4.0, 1.3, 1.0], 10.0, "rna_tangent_frame"),
        "Poly(A)-binding RBP": (0.76, [0.0, 0.0, 0.0], [-2.0, 1.2, 4.0], -8.0, "rna_tangent_frame"),
        "Ribosome large subunit": (0.90, [0.0, 0.0, 0.0], [0.0, 2.0, -3.2], 20.0, "rna_tangent_frame"),
        "Ribosome small subunit": (0.90, [0.0, 0.0, 0.0], [0.0, 5.4, 4.2], 20.0, "rna_tangent_frame"),
        "Standalone tRNA": (0.92, [0.0, 0.0, 0.0], [7.0, 4.8, 0.0], -20.0, "rna_tangent_frame"),
        "Actin protein": (1.0, [0.0, 0.0, 0.0], [0.0, 8.5, 15.0], 20.0, "manual_tangent_frame"),
    }
    for name, (fraction, offset, local_offset, roll, binding_mode) in mrna_anchor_fractions.items():
        asset = asset_by_name(manifest, name)
        set_path_anchor(asset, "mrna", fraction, offset, local_offset, roll, binding_mode)
        if name.startswith("Ribosome"):
            asset["rotation_deg"] = [74.0, 0.0, 20.0]
        elif name == "Actin protein":
            asset["rotation_deg"] = [75.0, 0.0, roll]
            asset["label_offset_mm"] = [0.0, 8.0, 0.4]
        else:
            asset["rotation_deg"] = [0.0, 0.0, roll]

    for segment in manifest["mrna"]["segments"]:
        segment["label_location_mm"] = [mrna_start[0] - 10.0, mrna_start[1] + 12.0, 0.3]

    return manifest


def write_manifest(base_manifest_path: Path = ROOT / "config" / "scene_manifest.json") -> dict:
    EXPERIMENT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = apply_arrangement_v1(json.loads(base_manifest_path.read_text(encoding="utf-8")))
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
