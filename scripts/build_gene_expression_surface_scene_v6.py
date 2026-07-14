#!/usr/bin/env python3
"""Build canonical V6 by applying an isolated gene-end layout to the V5 renderer."""

from __future__ import annotations

import copy
import json
import math
import os
import sys
from pathlib import Path

import bpy
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Vector


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import build_gene_expression_surface_scene_v5 as renderer  # noqa: E402
import canonical_v6_config as v6  # noqa: E402
import procedural_nucleic_geometry as geometry  # noqa: E402


ORIGINAL_VALIDATE = renderer.validate_v5_report
ORIGINAL_PLACE_OVERVIEW_LABELS = renderer.place_overview_labels
ORIGINAL_FIT_CAMERA = renderer.fit_camera_to_renderables
ORIGINAL_READER_ORDER_VALIDATION = renderer.reader_order_validation

DETAIL_PREVIEWS = {
    "full_overview": v6.OUTPUT_DIR / "preview_gene_expression_surface_style_v6_full_overview.png",
    "p53_dna": v6.OUTPUT_DIR / "preview_gene_expression_surface_style_v6_p53_dna.png",
    "polymerase_gene_end": v6.OUTPUT_DIR / "preview_gene_expression_surface_style_v6_polymerase_gene_end.png",
    "nucleosome_loop": v6.OUTPUT_DIR / "preview_gene_expression_surface_style_v6_nucleosome_loop.png",
    "ribosome_trna": v6.OUTPUT_DIR / "preview_gene_expression_surface_style_v6_ribosome_trna.png",
    "actin_product": v6.OUTPUT_DIR / "preview_gene_expression_surface_style_v6_actin_product.png",
    "cas9_dna": v6.OUTPUT_DIR / "preview_gene_expression_surface_style_v6_cas9_dna.png",
}

DETAIL_TITLES = {
    "p53_dna": "p53 tetramer + DNA (3TS8)",
    "nucleosome_loop": "nucleosome core + wrapped DNA (1AOI)",
    "polymerase_gene_end": "RNA polymerase II at gene end + nascent RNA 3′ end (2E2I)",
    "ribosome_trna": "ribosome + tRNA",
    "actin_product": "ACTB protein (1J6Z)",
    "cas9_dna": "Cas9 + guide/target DNA (4UN3)",
}

FOCUS_OBJECT_PATTERNS = {
    "p53_dna": ("p53 tetramer bound to DNA (3TS8)",),
    "nucleosome_loop": ("Nucleosome (1AOI)",),
    "polymerase_gene_end": ("RNA polymerase II elongation complex (2E2I)", "detail_polymerase_gene_end_rna"),
    "ribosome_trna": (
        "Ribosome small subunit (1J5E)",
        "Ribosome large subunit (1JJ2)",
        "Standalone tRNA (4TNA)",
        "detail_ribosome_mrna",
    ),
    "actin_product": ("Actin protein (1J6Z)",),
    "cas9_dna": ("Cas9 (4UN3)",),
}

PRIMARY_CALLOUTS = {
    "label_DNA_v6": {"text": "Actb promoter + gene DNA 3954 bp", "view_position": (0.825, 0.29), "span": (0.16, 0.42)},
    "label_mRNA_v6": {"text": "Actb mRNA 1852 nt", "view_position": (0.825, 0.63), "span": (0.44, 0.82)},
    "label_ACTB_primary_v6": {"text": "ACTB protein 375 aa", "view_position": (0.825, 0.92), "span": (0.87, 0.97)},
}
COMPACT_CALLOUT = {"object": "label_compact_mrna_v6", "text": "mRNA compact", "offset": (0.014, 0.014)}


def build_v6_dna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = renderer.direct_nucleic_meshes.build_dna_meshes(manifest, collections, materials)
    report["source_curve"] = renderer.add_source_path_curve("V6_DNA_source_path", path, collections["DNA"], materials["black"])
    start = path.points[0]
    renderer.base.create_text(
        "label_DNA_v6", "Actb promoter + gene DNA 3954 bp",
        (start.x + 4.0, start.y - 3.0, start.z + 0.25), 1.55,
        materials["black"], collections["Labels"],
    )
    return {"path": path, "report": report, "features": {"nucleosome_loop": report.get("nucleosome_loop")}}


def _rigid_translation_report(v6_path) -> dict:
    v5_manifest = json.loads(v6.V5_MANIFEST_PATH.read_text(encoding="utf-8"))
    v5_path = geometry.build_mrna_model(v5_manifest)["path"]
    delta = Vector(v6.BRANCH_TRANSLATION_MM)
    maximum_error = 0.0
    rows = []
    for index in range(21):
        fraction = index / 20.0
        old = v5_path.point_at_length(v5_path.length * fraction)
        new = v6_path.point_at_length(v6_path.length * fraction)
        error = (Vector((new.x, new.y, new.z)) - (Vector((old.x, old.y, old.z)) + delta)).length
        maximum_error = max(maximum_error, error)
        rows.append({"fraction": fraction, "error_mm": error})
    return {
        "basis": "21 centerline samples from tracked canonical V5 manifest",
        "translation_mm": list(v6.BRANCH_TRANSLATION_MM),
        "maximum_sample_error_mm": maximum_error,
        "tolerance_mm": 1e-5,
        "samples": rows,
        "passed": maximum_error <= 1e-5,
    }


def build_v6_mrna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = renderer.direct_nucleic_meshes.build_mrna_meshes(manifest, collections, materials)
    report["source_curve"] = renderer.add_source_path_curve("V6_mRNA_source_path", path, collections["mRNA"], materials["black"])
    report["rigid_branch_translation"] = _rigid_translation_report(path)
    start = path.points[0]
    renderer.base.create_text(
        "label_mRNA_v6", "Actb mRNA 1852 nt",
        (start.x + 3.0, start.y + 3.0, start.z + 0.3), 1.55,
        materials["black"], collections["Labels"],
    )
    return {"path": path, "report": report}


def build_v6_compact_mrna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = renderer.direct_nucleic_meshes.build_compact_mrna_meshes(manifest, collections, materials)
    center = renderer.contact_helpers.path_center(path)
    renderer.base.create_text(
        "label_compact_mrna_v6", "mRNA compact",
        (center.x, center.y + 3.2, center.z + 0.35), 1.45,
        materials["label_grey"], collections["Labels"],
    )
    return {"path": path, "report": report, "center_mm": [center.x, center.y, center.z]}


def source_curve_visibility_validation(report: dict) -> dict:
    rows = []
    failures = []
    for path_name, object_name in (("dna", "V6_DNA_source_path"), ("mrna", "V6_mRNA_source_path")):
        obj = bpy.data.objects.get(object_name)
        row = {
            "path": path_name,
            "object": object_name,
            "exists": obj is not None,
            "hide_viewport": bool(obj.hide_viewport) if obj else None,
            "hide_render": bool(obj.hide_render) if obj else None,
        }
        rows.append(row)
        if obj is None or not obj.hide_viewport or not obj.hide_render:
            failures.append({"object": object_name, "reason": "missing_or_visible_source_curve"})
    return {"policy": "V6 source paths are retained but hidden", "rows": rows, "failures": failures}


def reader_order_validation(dna_path) -> dict:
    result = ORIGINAL_READER_ORDER_VALIDATION(dna_path)
    result["path_mode"] = "v6_reader_order_serpentine_with_gene_end_polymerase"
    return result


def dna_cluster_validation(asset_reports: list[dict]) -> dict:
    early_names = sorted(v6.DNA_STRICT_CONTACT_NAMES - {"RNA polymerase II elongation complex"})
    by_name = {item["name"]: item for item in asset_reports}
    rows = []
    failures = []
    for name in early_names:
        item = by_name.get(name)
        attachment = (item or {}).get("attachment_empty") or {}
        fraction = attachment.get("fraction")
        rows.append({"name": name, "fraction": fraction, "max_allowed_fraction": 0.205})
        if item is None or attachment.get("path") != "dna" or fraction is None or float(fraction) > 0.2055:
            failures.append({"name": name, "reason": "outside_early_gene_cluster", "fraction": fraction})
    return {"policy": "all DNA regulatory binders except gene-end Pol II remain in the first 20.5 percent", "rows": rows, "failures": failures}


def add_v6_cameras(manifest, dna_path, mrna_path, compact_path, asset_reports, collections, materials) -> dict[str, str]:
    renderer.base.add_lighting_and_camera(collections, materials)
    dna_center = renderer.contact_helpers.path_center(dna_path)
    mrna_center = renderer.contact_helpers.path_center(mrna_path)
    polymerase = renderer._asset_report_location(asset_reports, "RNA polymerase II elongation complex")
    p53 = renderer._asset_report_location(asset_reports, "p53 tetramer bound to DNA")
    nucleosome = renderer._asset_report_location(asset_reports, "Nucleosome")
    ribosome = sum((renderer._asset_report_location(asset_reports, name) for name in (
        "Ribosome large subunit", "Ribosome small subunit", "Standalone tRNA"
    )), Vector()) / 3.0
    actin = renderer._asset_report_location(asset_reports, "Actin protein")
    cas9 = renderer._asset_report_location(asset_reports, "Cas9")
    scene_center = (dna_center + mrna_center + actin) / 3.0
    names = {
        "full_overview": renderer.create_camera(
            "Camera_v6_full_overview",
            (scene_center.x - 54.0, scene_center.y - 142.0, scene_center.z + 62.0),
            scene_center, 205.0,
        ),
        "p53_dna": renderer.create_camera("Camera_v6_p53_dna", (p53.x - 8.0, p53.y - 12.0, p53.z + 9.0), p53, 9.0),
        "polymerase_gene_end": renderer.create_camera(
            "Camera_v6_polymerase_gene_end", (polymerase.x - 10.0, polymerase.y - 14.0, polymerase.z + 10.0), polymerase, 11.0
        ),
        "nucleosome_loop": renderer.create_camera(
            "Camera_v6_nucleosome_loop", (nucleosome.x - 9.0, nucleosome.y - 13.0, nucleosome.z + 10.0), nucleosome, 10.0
        ),
        "ribosome_trna": renderer.create_camera(
            "Camera_v6_ribosome_trna", (ribosome.x - 13.0, ribosome.y - 18.0, ribosome.z + 15.0), ribosome, 17.0
        ),
        "actin_product": renderer.create_camera(
            "Camera_v6_actin_product", (actin.x - 7.0, actin.y - 10.0, actin.z + 8.0), actin, 7.0
        ),
        "cas9_dna": renderer.create_camera(
            "Camera_v6_cas9_dna", (cas9.x - 9.0, cas9.y - 13.0, cas9.z + 10.0), cas9, renderer.DETAIL_SHARED_ORTHO_SCALE_MM
        ),
    }
    overview_camera = bpy.data.objects[names["full_overview"]]
    overview_camera["overview_elevation_deg"] = math.degrees(math.atan2(62.0, math.hypot(54.0, 142.0)))
    overview_camera["overview_target_mm"] = list(scene_center)
    bpy.context.scene.camera = bpy.data.objects[names["full_overview"]]
    return names


def add_detail_context_curves(mrna_path, collections: dict, materials: dict) -> dict:
    collection = collections.get("Detail context")
    if collection is None:
        collection = bpy.data.collections.new("Detail context")
        bpy.context.scene.collection.children.link(collection)
        collections["Detail context"] = collection
    return {
        "polymerase_gene_end": renderer._detail_path_curve(
            "detail_polymerase_gene_end_rna", mrna_path, 0.0, 0.055, materials["rna_gold"], collection
        ),
        "ribosome_trna": renderer._detail_path_curve(
            "detail_ribosome_mrna", mrna_path, 0.855, 0.935, materials["rna_gold"], collection
        ),
    }


def place_overview_labels(camera_name: str, collections: dict, materials: dict) -> dict:
    dummy = None
    if bpy.data.objects.get("label_ACTB_primary_v5") is None:
        dummy = bpy.data.objects.new("label_ACTB_primary_v5", None)
        bpy.context.scene.collection.objects.link(dummy)
        dummy.hide_render = True
    if bpy.data.objects.get("label_ACTB_primary_v6") is None:
        obj = renderer.base.create_text(
            "label_ACTB_primary_v6", "ACTB protein 375 aa", (0.0, 0.0, 0.0), 1.55,
            materials["black"], collections["Labels"], align="LEFT",
        )
        obj.rotation_euler = bpy.data.objects[camera_name].rotation_euler
    result = ORIGINAL_PLACE_OVERVIEW_LABELS(camera_name, collections, materials)
    if dummy is not None:
        bpy.data.objects.remove(dummy, do_unlink=True)
    return result


def fit_camera_to_renderables(camera_name: str, margin_fraction: float = 0.075) -> dict:
    result = ORIGINAL_FIT_CAMERA(camera_name, margin_fraction)
    if camera_name == "Camera_v6_full_overview" and result.get("fit_applied"):
        camera = bpy.data.objects[camera_name]
        camera.data.ortho_scale *= 1.08
        camera.location += camera.matrix_world.to_3x3() @ Vector((float(camera.data.ortho_scale) * 0.065, 0.0, 0.0))
        result["new_location_mm"] = list(camera.location)
        result["new_ortho_scale_mm"] = float(camera.data.ortho_scale)
        result["v6_right_gutter"] = True
    return result


def projected_centerline_separation_validation() -> dict:
    camera = bpy.data.objects["Camera_v6_full_overview"]
    scene = bpy.context.scene

    def projected_points(object_name: str):
        obj = bpy.data.objects[object_name]
        return [
            world_to_camera_view(scene, camera, obj.matrix_world @ point.co.xyz)
            for spline in obj.data.splines
            for point in spline.points
        ]

    bpy.context.view_layer.update()
    dna_points = projected_points("V6_DNA_source_path")
    mrna_points = projected_points("V6_mRNA_source_path")
    attachment_fraction = 0.02
    attachment_sample_count = math.ceil(len(mrna_points) * attachment_fraction)
    evaluated_mrna_points = mrna_points[attachment_sample_count:]
    nearest_distances = [
        min(math.hypot(mrna.x - dna.x, mrna.y - dna.y) for dna in dna_points)
        for mrna in evaluated_mrna_points
    ]
    ordered = sorted(nearest_distances)
    midpoint = len(ordered) // 2
    median = (
        ordered[midpoint]
        if len(ordered) % 2
        else (ordered[midpoint - 1] + ordered[midpoint]) * 0.5
    )
    close_fraction = sum(distance < 0.02 for distance in nearest_distances) / len(nearest_distances)
    median_minimum = 0.10
    close_fraction_maximum = 0.31
    passed = median >= median_minimum and close_fraction <= close_fraction_maximum
    return {
        "policy": "sampled DNA and mRNA centerlines remain visually separated except at their biological Pol II attachment",
        "camera": camera.name,
        "camera_elevation_deg": float(camera.get("overview_elevation_deg", 0.0)),
        "dna_sample_count": len(dna_points),
        "mrna_sample_count": len(mrna_points),
        "excluded_attachment_fraction": attachment_fraction,
        "excluded_attachment_sample_count": attachment_sample_count,
        "evaluated_mrna_sample_count": len(evaluated_mrna_points),
        "median_normalized_separation": median,
        "minimum_median_normalized_separation": median_minimum,
        "fraction_within_0_02_normalized_units": close_fraction,
        "maximum_fraction_within_0_02_normalized_units": close_fraction_maximum,
        "previous_v6_fraction_within_0_02_normalized_units": 0.40274885464389837,
        "minimum_normalized_separation": min(nearest_distances),
        "passed": passed,
    }


def validate_v6_report(report: dict) -> None:
    report["canonical_version"] = "v6"
    failures = []
    dna_end = Vector(report["reader_order_validation"]["end_mm"])
    mrna_origin = Vector(report["mrna"]["segments"][0]["start_mm"])
    origin_error = (dna_end - mrna_origin).length
    report["gene_end_transcription_validation"] = {
        "dna_end_mm": list(dna_end),
        "mrna_origin_mm": list(mrna_origin),
        "origin_error_mm": origin_error,
        "origin_tolerance_mm": 1e-5,
    }
    if origin_error > 1e-5:
        failures.append({"reason": "mrna_origin_not_at_gene_end", "error_mm": origin_error})

    assets = {item["name"]: item for item in report.get("pdb_assets", [])}
    pol = assets.get("RNA polymerase II elongation complex") or {}
    attachment = pol.get("attachment_empty") or {}
    fraction = attachment.get("fraction")
    if fraction is None or not math.isclose(float(fraction), 1.0, abs_tol=1e-6):
        failures.append({"reason": "polymerase_not_at_gene_end", "fraction": fraction})
    rigid = report.get("mrna", {}).get("rigid_branch_translation", {})
    if not rigid.get("passed"):
        failures.append({"reason": "mrna_not_rigid_v5_translation", "validation": rigid})
    origin_distance = (report.get("polymerase_rna_origin") or {}).get("distance_mm")
    if origin_distance is None or float(origin_distance) > 0.35:
        failures.append({"reason": "polymerase_not_at_rna_origin", "distance_mm": origin_distance})

    v5_manifest = json.loads(v6.V5_MANIFEST_PATH.read_text(encoding="utf-8"))
    expected_fractions = {
        asset["name"]: float(asset["path_anchor"]["fraction"])
        for asset in v5_manifest.get("pdb_assets", [])
        if (asset.get("path_anchor") or {}).get("path") == "mrna"
    }
    attachment_rows = []
    for name, expected in expected_fractions.items():
        actual = ((assets.get(name) or {}).get("attachment_empty") or {}).get("fraction")
        passed = actual is not None and math.isclose(float(actual), expected, abs_tol=1e-6)
        attachment_rows.append({"name": name, "expected_fraction": expected, "actual_fraction": actual, "passed": passed})
        if not passed:
            failures.append({"reason": "changed_rna_attachment_fraction", "name": name, "expected": expected, "actual": actual})
    report["branch_attachment_validation"] = {"rows": attachment_rows, "passed": all(row["passed"] for row in attachment_rows)}
    separation = projected_centerline_separation_validation()
    report["overview_projected_separation_validation"] = separation
    if not separation["passed"]:
        failures.append({"reason": "insufficient_overview_dna_mrna_separation", "validation": separation})

    compatibility = copy.deepcopy(report)
    compatibility["source_manifest"] = str(v6.V5_MANIFEST_PATH)
    compatibility["mrna"]["path_origin"] = "3_prime_at_polymerase_ii"
    detail = compatibility.get("detail_render_specs", {})
    if "polymerase_gene_end" in detail:
        detail["polymerase_rna_start"] = detail.pop("polymerase_gene_end")
    ORIGINAL_VALIDATE(compatibility)
    if failures:
        raise RuntimeError(f"Canonical V6 validation failed: {failures}")


def configure_renderer() -> None:
    renderer.v5 = v6
    renderer.OUTPUT_DIR = v6.OUTPUT_DIR
    renderer.BLEND_PATH = v6.BLEND_PATH
    renderer.PREVIEW_PATH = v6.PREVIEW_PATH
    renderer.REPORT_PATH = v6.REPORT_PATH
    renderer.DETAIL_PREVIEWS = DETAIL_PREVIEWS
    renderer.DETAIL_TITLES = DETAIL_TITLES
    renderer.FOCUS_OBJECT_PATTERNS = FOCUS_OBJECT_PATTERNS
    renderer.PRIMARY_CALLOUTS = PRIMARY_CALLOUTS
    renderer.COMPACT_CALLOUT = COMPACT_CALLOUT
    renderer.PROTEIN_AA_CONTOUR_NM = v6.PROTEIN_AA_CONTOUR_NM
    renderer.build_v5_dna = build_v6_dna
    renderer.build_v5_mrna = build_v6_mrna
    renderer.build_v5_compact_mrna = build_v6_compact_mrna
    renderer.source_curve_visibility_validation = source_curve_visibility_validation
    renderer.reader_order_validation = reader_order_validation
    renderer.dna_cluster_validation = dna_cluster_validation
    renderer.add_v5_cameras = add_v6_cameras
    renderer.add_detail_context_curves = add_detail_context_curves
    renderer.place_overview_labels = place_overview_labels
    renderer.fit_camera_to_renderables = fit_camera_to_renderables
    renderer.validate_v5_report = validate_v6_report


def finalize_v6_identifiers() -> None:
    """Remove renderer-internal V5 identifiers from the saved V6 scene/report."""
    rename_map = {}
    datablock_groups = (
        bpy.data.objects,
        bpy.data.meshes,
        bpy.data.materials,
        bpy.data.lights,
        bpy.data.cameras,
    )
    for group in datablock_groups:
        for datablock in list(group):
            old_name = datablock.name
            if old_name == "Camera":
                datablock.name = "Camera_v6_base"
                rename_map[old_name] = datablock.name
            elif old_name.startswith("v5_"):
                datablock.name = "v6_" + old_name[3:]
                rename_map[old_name] = datablock.name

    for obj in bpy.data.objects:
        for key in list(obj.keys()):
            if key.startswith("v5_"):
                obj["v6_" + key[3:]] = obj[key]
                del obj[key]
        for modifier in obj.modifiers:
            if modifier.name.startswith("v5_"):
                modifier.name = "v6_" + modifier.name[3:]

    bpy.context.scene["canonical_version"] = "v6"

    def rewrite_names(value):
        if isinstance(value, dict):
            return {key: rewrite_names(item) for key, item in value.items()}
        if isinstance(value, list):
            return [rewrite_names(item) for item in value]
        if isinstance(value, str):
            return rename_map.get(value, value)
        return value

    report = json.loads(v6.REPORT_PATH.read_text(encoding="utf-8"))
    v6.REPORT_PATH.write_text(json.dumps(rewrite_names(report), indent=2), encoding="utf-8")
    bpy.ops.wm.save_as_mainfile(filepath=str(v6.BLEND_PATH))


if __name__ == "__main__":
    configure_renderer()
    os.environ.pop("V5_OVERVIEW_ONLY", None)
    renderer.main()
    finalize_v6_identifiers()
