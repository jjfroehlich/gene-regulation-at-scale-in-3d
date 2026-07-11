#!/usr/bin/env python3
"""Deterministic, scale-accurate scene-RNA surfaces and compact mRNP variants."""

from __future__ import annotations

import copy
import json
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))

import procedural_nucleic_geometry as geom  # noqa: E402


MANIFEST_PATH = ROOT / "config" / "scene_manifest_v5.json"
TOTAL_NT = 1852


def load_manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def crescent_guide() -> list[geom.Vec]:
    return [
        geom.Vec(
            (3.5 + 0.45 * math.sin(math.tau * 5 * i / 420)) * math.cos(math.radians(30 + 300 * i / 420)),
            0.78 * (3.5 + 0.45 * math.sin(math.tau * 5 * i / 420)) * math.sin(math.radians(30 + 300 * i / 420)),
            1.1 * math.sin(math.tau * 3 * i / 420),
        )
        for i in range(421)
    ]


def bilobal_guide() -> list[geom.Vec]:
    return [
        geom.Vec(
            3.7 * math.sin(math.tau * 4 * i / 600),
            2.9 * math.sin(math.tau * 5 * i / 600 + 0.55),
            1.5 * math.sin(math.tau * 7 * i / 600 + 0.2),
        )
        for i in range(601)
    ]


def variant_specs(manifest: dict) -> list[dict]:
    mrna = manifest["procedural_nucleic_acids"]["mrna"]
    conical, _ = geom.conical_spiral_mrna_points(manifest, mrna["canonical_spiral"])
    rosette = geom.compact_mrnp_path(
        center=(0, 0, 0), target_length=TOTAL_NT * manifest["units"]["mrna_nt_to_mm"],
        width=10.5, height=8.0, depth=3.4, points_count=2500, seed=304,
    )
    surfaces = [
        ("surface_smooth_tube", "Smooth backbone surface", "smooth", 0.145,
         "A clean, continuous rounded envelope with no explicit nucleotide relief."),
        ("surface_soft_molecular", "Soft molecular surface", "soft_molecular", 0.155,
         "Low-frequency, shallow radius variation suggests a hydrated molecular surface without spikes."),
        ("surface_twisted_groove", "Subtle twisted-groove surface", "twisted_groove", 0.150,
         "A shallow rotating oval groove adds directionality while keeping the silhouette smooth."),
    ]
    specs = [dict(key=k, group="elongated", title=t, intent=intent, guide=conical,
                  compact=False, surface_mode=mode, tube_radius_mm=radius, seed=seed)
             for seed, (k, t, mode, radius, intent) in enumerate(surfaces, 7101)]
    specs.extend([
        {
            "key": "compact_rosette", "group": "compact", "title": "RNA-only compact rosette",
            "intent": "A globular RNA-only polymer with many local hairpins; retained as a protein-free baseline.",
            "guide": rosette, "compact": True, "stems": 38, "paired_fraction": 0.58, "seed": 3838,
            "architecture": "rna_only_rosette", "protein_lobes": [],
        },
        {
            "key": "compact_crescent", "group": "compact", "title": "RNA-only compact crescent",
            "intent": "A crescent-like RNA-only compact polymer; retained as a protein-free silhouette baseline.",
            "guide": crescent_guide(), "compact": True, "stems": 42, "paired_fraction": 0.64, "seed": 4242,
            "architecture": "rna_only_crescent", "protein_lobes": [],
        },
        {
            "key": "compact_multi_domain", "group": "compact", "title": "RNA-only multi-domain globule",
            "intent": "A densely paired RNA-only stress test with several interwoven domains.",
            "guide": bilobal_guide(), "compact": True, "stems": 48, "paired_fraction": 0.70, "seed": 4747,
            "architecture": "rna_only_multi_domain", "protein_lobes": [],
        },
        {
            "key": "rnp_ejc_clamped", "group": "compact", "title": "EJC-inspired clamped mRNP",
            "intent": "A compact RNA rosette carrying several localized protein clamps, inspired by EJC-bound mRNA.",
            "guide": rosette, "compact": True, "stems": 38, "paired_fraction": 0.58, "seed": 3838,
            "architecture": "distributed_local_clamps",
            "protein_lobes": [(-3.0, 1.8, 0.3, 1.8, 1.2, 1.0, -18), (0.0, 3.0, 0.5, 1.6, 1.1, 1.0, 12),
                              (3.1, 1.0, -0.2, 1.8, 1.25, 1.0, 28), (2.3, -2.4, 0.5, 1.5, 1.0, 0.9, -12),
                              (-2.4, -2.2, -0.3, 1.55, 1.05, 0.9, 20)],
        },
        {
            "key": "rnp_srp_scaffold", "group": "compact", "title": "SRP-inspired scaffold mRNP",
            "intent": "A crescent RNA scaffold embraced by elongated protein saddles, inspired by SRP RNA-protein recognition surfaces.",
            "guide": crescent_guide(), "compact": True, "stems": 42, "paired_fraction": 0.64, "seed": 4242,
            "architecture": "rna_scaffold_with_protein_saddles",
            "protein_lobes": [(-1.6, 0.7, 0.2, 3.3, 1.45, 1.15, 25), (2.4, -0.7, -0.1, 2.4, 1.25, 1.0, -32),
                              (0.5, 2.7, 0.5, 1.45, 1.0, 0.9, 8)],
        },
        {
            "key": "rnp_telomerase_bilobal", "group": "compact", "title": "Telomerase-inspired bilobal mRNP",
            "intent": "Two asymmetric protein-rich lobes bridged and wrapped by structured RNA, inspired by telomerase RNP architecture.",
            "guide": bilobal_guide(), "compact": True, "stems": 46, "paired_fraction": 0.68, "seed": 4848,
            "architecture": "bilobal_protein_rich_particle",
            "protein_lobes": [(-3.2, 0.3, -0.35, 2.35, 1.9, 1.35, 12), (-1.7, 1.45, -0.15, 1.65, 1.35, 1.05, -20),
                              (3.1, -0.2, -0.3, 2.45, 1.8, 1.4, -15), (1.5, -1.55, -0.45, 1.55, 1.2, 1.0, 22)],
        },
    ])
    return specs


def _bbox(points: list[geom.Vec]) -> dict:
    mins = [min(getattr(p, axis) for p in points) for axis in ("x", "y", "z")]
    maxs = [max(getattr(p, axis) for p in points) for axis in ("x", "y", "z")]
    return {"min_mm": mins, "max_mm": maxs, "bbox_mm": [b - a for a, b in zip(mins, maxs)]}


def build_variant(manifest: dict, spec: dict) -> dict:
    nt_to_mm = manifest["units"]["mrna_nt_to_mm"]
    if spec["group"] == "elongated":
        raw = [geom.Vec(*p) if not isinstance(p, geom.Vec) else p for p in spec["guide"]]
        scale = (TOTAL_NT * nt_to_mm) / geom.polyline_length(raw)
        origin = raw[0]
        points = [origin + (p - origin) * scale for p in raw]
        segments = geom.split_path_by_lengths(points, manifest["mrna"]["segments"], nt_to_mm)
        report = {
            "allocation_check_nt": sum(x["nt"] for x in manifest["mrna"]["segments"]),
            "measured_length_mm": geom.polyline_length(points), "stem_count": 0, "base_pair_bridge_count": 0,
            "paired_fraction": 0.0, "bbox": _bbox(points), "surface_mode": spec["surface_mode"],
            "tube_radius_mm": spec["tube_radius_mm"], "explicit_stem_loops": False,
            "controlled_centerline": "identical canonical conical centerline for all three surface candidates",
        }
        base_pairs = []
    else:
        settings = copy.deepcopy(manifest["procedural_nucleic_acids"]["mrna"])
        settings["secondary_structure"] = {
            "model": "deterministic_schematic_stem_loop", "sequence_resolved": False,
            "elongated_stem_count": spec["stems"], "compact_stem_count": spec["stems"],
            "elongated_paired_fraction_target": spec["paired_fraction"], "compact_paired_fraction_target": spec["paired_fraction"],
            "stem_bp_min": 6, "stem_bp_max": 18, "hairpin_loop_nt_min": 4, "hairpin_loop_nt_max": 12,
            "a_form_bp_per_turn": 11.0, "a_form_diameter_mm": 0.92,
            "elongated_seed": spec["seed"], "compact_seed": spec["seed"],
        }
        structure = geom.structured_stem_loop_path(spec["guide"], TOTAL_NT, nt_to_mm, settings, compact=True)
        points, base_pairs = structure["points"], structure["base_pairs"]
        segments = geom.split_path_by_lengths(points, manifest["mrna"]["segments"], nt_to_mm)
        report = copy.deepcopy(structure["report"])
        report.update({"architecture": spec["architecture"], "protein_lobe_count": len(spec["protein_lobes"]),
                       "schematic_rnp": bool(spec["protein_lobes"]), "sequence_resolved": False})
    report.update({"key": spec["key"], "group": spec["group"], "title": spec["title"], "intent": spec["intent"],
                   "segments": [x["segment"] for x in segments], "end_to_end_mm": (points[-1] - points[0]).length})
    return {"key": spec["key"], "group": spec["group"], "title": spec["title"], "points": points,
            "base_pairs": base_pairs, "segments": segments, "report": report,
            "protein_lobes": spec.get("protein_lobes", [])}


def build_all_variants(manifest: dict | None = None) -> list[dict]:
    manifest = manifest or load_manifest()
    return [build_variant(manifest, spec) for spec in variant_specs(manifest)]


if __name__ == "__main__":
    print(json.dumps([v["report"] for v in build_all_variants()], indent=2))
