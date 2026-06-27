#!/usr/bin/env python3
"""Build the arrangement V1 experiment scene without touching canonical outputs."""

from __future__ import annotations

import json
import importlib
import os
import sys
from pathlib import Path

import bpy
from mathutils import Vector


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

arrangement = importlib.import_module(os.environ.get("ARRANGEMENT_CONFIG_MODULE", "arrangement_v1_config"))  # noqa: E402
import blender_nucleic_meshes as direct_nucleic_meshes  # noqa: E402
import build_gene_expression_scene as base  # noqa: E402
import build_gene_expression_surface_scene as scene  # noqa: E402
import procedural_nucleic_geometry as nucleic_geometry  # noqa: E402


EXPERIMENT_KEY = getattr(arrangement, "EXPERIMENT_KEY", "arrangement_v1")
SCENE_BASENAME = getattr(arrangement, "SCENE_BASENAME", f"gene_expression_{EXPERIMENT_KEY}")
REPORT_KIND = getattr(arrangement, "REPORT_KIND", f"{EXPERIMENT_KEY}_experiment")
BLEND_PATH = arrangement.OUTPUT_DIR / f"{SCENE_BASENAME}.blend"
PREVIEW_PATH = arrangement.OUTPUT_DIR / f"preview_{SCENE_BASENAME}.png"
REPORT_PATH = arrangement.OUTPUT_DIR / f"{SCENE_BASENAME}_report.json"
DETAIL_PREVIEWS = {
    "polymerase_rna_start": arrangement.OUTPUT_DIR / f"preview_{EXPERIMENT_KEY}_polymerase_rna_start.png",
    "nucleosome_loop": arrangement.OUTPUT_DIR / f"preview_{EXPERIMENT_KEY}_nucleosome_loop.png",
    "mrna_spiral": arrangement.OUTPUT_DIR / f"preview_{EXPERIMENT_KEY}_mrna_spiral.png",
    "ribosome_top": arrangement.OUTPUT_DIR / f"preview_{EXPERIMENT_KEY}_ribosome_top.png",
}


def path_center(path) -> Vector:
    points = path.points
    min_v = Vector((min(point.x for point in points), min(point.y for point in points), min(point.z for point in points)))
    max_v = Vector((max(point.x for point in points), max(point.y for point in points), max(point.z for point in points)))
    return (min_v + max_v) * 0.5


def add_source_path_curve(name: str, path, collection, material) -> dict:
    points = [(point.x, point.y, point.z + 0.16) for point in path.points]
    obj = base.create_curve(name, points, 0.028, material, collection, resolution=1)
    obj.hide_render = True
    obj["source_path_for_rebake"] = True
    obj["sampled_points"] = len(points)
    obj["curve_role"] = "editable_source_path"
    return {"object": obj.name, "points": len(points), "hide_render": True}


def build_experiment_dna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = direct_nucleic_meshes.build_dna_meshes(manifest, collections, materials)
    source_curve_report = add_source_path_curve("DNA_source_path", path, collections["DNA"], materials["black"])
    report["source_curve"] = source_curve_report
    base.create_text(f"label_DNA_{EXPERIMENT_KEY}", "full-gene ACTB serpentine DNA base", (-58.0, -55.0, 0.2), 3.0, materials["black"], collections["Labels"])
    return {"path": path, "report": report, "features": {"nucleosome_loop": report.get("nucleosome_loop")}}


def build_experiment_mrna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = direct_nucleic_meshes.build_mrna_meshes(manifest, collections, materials)
    source_curve_report = add_source_path_curve("MRNA_source_path", path, collections["mRNA"], materials["black"])
    report["source_curve"] = source_curve_report
    start = path.points[0]
    base.create_text(f"label_mRNA_{EXPERIMENT_KEY}", "mRNA spiral from Pol II", (start.x + 4.0, start.y + 8.0, start.z + 0.3), 2.6, materials["black"], collections["Labels"])
    return {"path": path, "report": report}


def add_experiment_cameras(manifest, mrna_path) -> dict[str, str]:
    mrna_center = path_center(mrna_path)
    ribosome_point = mrna_path.point_at_length(mrna_path.length * 0.88)
    actin_point = mrna_path.point_at_length(mrna_path.length)
    loop_report = nucleic_geometry.dna_nucleosome_loop_report(manifest) or {}
    loop_center = loop_report.get("center_mm") or [56.5, 27.0, 0.0]
    full_target = Vector((0.0, -8.0, 26.0))
    camera_specs = {
        f"Camera_{EXPERIMENT_KEY}_full": ((-92.0, -132.0, 92.0), full_target, 152.0),
        f"Camera_{EXPERIMENT_KEY}_polymerase_rna_start": (
            (mrna_path.points[0].x - 28.0, mrna_path.points[0].y - 42.0, mrna_path.points[0].z + 31.0),
            Vector((mrna_path.points[0].x, mrna_path.points[0].y, mrna_path.points[0].z + 4.0)),
            34.0,
        ),
        f"Camera_{EXPERIMENT_KEY}_nucleosome_loop": (
            (loop_center[0], loop_center[1], loop_center[2] + 45.0 if len(loop_center) > 2 else 45.0),
            Vector((loop_center[0], loop_center[1], loop_center[2] if len(loop_center) > 2 else 0.0)),
            9.5,
        ),
        f"Camera_{EXPERIMENT_KEY}_mrna_spiral": ((mrna_center.x - 46.0, mrna_center.y - 70.0, mrna_center.z + 54.0), mrna_center, 92.0),
        f"Camera_{EXPERIMENT_KEY}_ribosome_top": (
            ((ribosome_point.x + actin_point.x) * 0.5 - 28.0, (ribosome_point.y + actin_point.y) * 0.5 - 42.0, (ribosome_point.z + actin_point.z) * 0.5 + 34.0),
            Vector(((ribosome_point.x + actin_point.x) * 0.5, (ribosome_point.y + actin_point.y) * 0.5, (ribosome_point.z + actin_point.z) * 0.5)),
            54.0,
        ),
    }
    for name, (location, target, ortho_scale) in camera_specs.items():
        camera_data = bpy.data.cameras.new(name)
        camera = bpy.data.objects.new(name, camera_data)
        camera.location = location
        direction = target - camera.location
        camera.rotation_euler = direction.to_track_quat("-Z", "Y").to_euler()
        camera_data.type = "ORTHO"
        camera_data.ortho_scale = ortho_scale
        bpy.context.scene.collection.objects.link(camera)
    bpy.context.scene.camera = bpy.data.objects[f"Camera_{EXPERIMENT_KEY}_full"]
    return {name: camera_name for name, camera_name in {
        "full": f"Camera_{EXPERIMENT_KEY}_full",
        "polymerase_rna_start": f"Camera_{EXPERIMENT_KEY}_polymerase_rna_start",
        "nucleosome_loop": f"Camera_{EXPERIMENT_KEY}_nucleosome_loop",
        "mrna_spiral": f"Camera_{EXPERIMENT_KEY}_mrna_spiral",
        "ribosome_top": f"Camera_{EXPERIMENT_KEY}_ribosome_top",
    }.items()}


def hide_non_scale_labels_for_render() -> int:
    hidden = 0
    for obj in bpy.data.objects:
        if obj.type == "FONT" and not obj.name.startswith("label_scale_"):
            obj.hide_render = True
            hidden += 1
    return hidden


def attachment_table(asset_reports: list[dict]) -> list[dict]:
    rows = []
    for item in asset_reports:
        attachment = item.get("attachment_empty")
        if not attachment:
            continue
        rows.append(
            {
                "name": item["name"],
                "pdb_id": item["pdb_id"],
                "empty": attachment["empty"],
                "path": attachment["path"],
                "distance_mm": attachment["distance_mm"],
                "fraction": attachment["fraction"],
                "offset_local_mm": attachment["offset_local_mm"],
                "roll_deg": attachment["roll_deg"],
                "binding_mode": attachment["binding_mode"],
            }
        )
    return rows


def polymerase_rna_origin_report(asset_reports: list[dict], mrna_path) -> dict | None:
    polymerase = next((item for item in asset_reports if item.get("name") == "RNA polymerase II elongation complex"), None)
    if not polymerase:
        return None
    polymerase_location = Vector(polymerase["location_mm"])
    start = mrna_path.points[0]
    rna_start = Vector((start.x, start.y, start.z))
    distance = (polymerase_location - rna_start).length
    return {
        "polymerase_location_mm": [polymerase_location.x, polymerase_location.y, polymerase_location.z],
        "rna_start_mm": [rna_start.x, rna_start.y, rna_start.z],
        "distance_mm": distance,
    }


def arrangement_binding_checks(asset_reports: list[dict], mrna_path) -> dict:
    dna_names = {
        "RNA polymerase II elongation complex",
        "Nucleosome",
        "Cas9",
        "Transcription factor 1",
        "Transcription factor 3",
        "Transcription factor 4",
        "p53 tetramer bound to DNA",
    }
    rna_names = {
        "Pumilio RBP",
        "MS2 coat protein MCP",
        "mCherry/RFP tag",
        "Argonaute",
        "HuR-like RBP",
        "Poly(A)-binding RBP",
        "Ribosome large subunit",
        "Ribosome small subunit",
        "Standalone tRNA",
    }
    by_name = {item["name"]: item for item in asset_reports}
    dna_paths = {
        name: ((by_name.get(name) or {}).get("attachment_empty") or {}).get("path")
        for name in sorted(dna_names)
    }
    rna_paths = {
        name: ((by_name.get(name) or {}).get("attachment_empty") or {}).get("path")
        for name in sorted(rna_names)
    }
    mrna_z_values = [point.z for point in mrna_path.points]
    mrna_min_z = min(mrna_z_values)
    mrna_max_z = max(mrna_z_values)
    top_quartile_z = mrna_min_z + 0.75 * (mrna_max_z - mrna_min_z)
    ribosome_locations = {
        name: by_name[name]["location_mm"]
        for name in ("Ribosome large subunit", "Ribosome small subunit")
        if name in by_name
    }
    actin = by_name.get("Actin protein")
    actin_location = actin.get("location_mm") if actin else None
    return {
        "dna_binding_assets_on_dna_path": dna_paths,
        "rna_binding_assets_on_mrna_path": rna_paths,
        "dna_binding_all_on_dna": all(path == "dna" for path in dna_paths.values()),
        "rna_binding_all_on_mrna": all(path == "mrna" for path in rna_paths.values()),
        "mrna_z_extent_mm": mrna_max_z - mrna_min_z,
        "mrna_top_quartile_z_mm": top_quartile_z,
        "ribosome_locations_mm": ribosome_locations,
        "ribosome_in_top_quartile": all(location[2] >= top_quartile_z for location in ribosome_locations.values()),
        "actin_location_mm": actin_location,
        "actin_above_rna": bool(actin_location and actin_location[2] > mrna_max_z),
    }


def emphasize_nucleosome_core(materials: dict) -> list[str]:
    emphasized = []
    material = materials["tf_purple"]
    for obj in bpy.data.objects:
        if obj.name.startswith("Nucleosome (1AOI) surface protein"):
            obj.data.materials.clear()
            obj.data.materials.append(material)
            obj["arrangement_v1_emphasized_core_material"] = True
            emphasized.append(obj.name)
    return emphasized


def main() -> None:
    arrangement.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if not arrangement.MANIFEST_PATH.exists():
        arrangement.write_manifest()
    manifest = json.loads(arrangement.MANIFEST_PATH.read_text(encoding="utf-8"))

    scene.clean_scene()
    scene.configure_scene()
    collections = base.make_collections(manifest["collections"])
    materials = base.make_materials()
    scene.soften_materials(materials)

    dna = build_experiment_dna(manifest, collections, materials)
    mrna = build_experiment_mrna(manifest, collections, materials)
    path_context = {
        "dna": dna["path"],
        "dna_features": dna["features"],
        "dna_anchor_mode": "closest_xy",
        "mrna": mrna["path"],
    }

    report = {
        "title": manifest["title"],
        "kind": REPORT_KIND,
        "units": manifest["units"],
        "source_manifest": str(arrangement.MANIFEST_PATH),
        "surface_asset_dir": str(scene.REDUCED_SURFACE_DIR),
        "outputs": {
            "blend": str(BLEND_PATH),
            "preview": str(PREVIEW_PATH),
            "report": str(REPORT_PATH),
            "detail_previews": {name: str(path) for name, path in DETAIL_PREVIEWS.items()},
        },
    }
    report["layout_intent"] = {
        "dna": "full ACTB gene plus 500 bp promoter folded into a rounded serpentine DNA base with segment colors",
        "mrna": "full-length actin mRNA starts at RNA polymerase II and follows a tightening irregular spiral upward in Z",
        "translation": "ribosome is anchored near the top of the mRNA spiral and actin floats above the mRNA terminus",
        "nucleic_acid_pipeline": "direct Blender mesh generation for DNA and mRNA",
    }
    report["layout_intent"].update(getattr(arrangement, "LAYOUT_INTENT_OVERRIDES", {}))
    report["dna"] = dna["report"]
    report["mrna"] = mrna["report"]
    report["scale_bars"] = scene.add_scale_bars(manifest, collections, materials)
    report["pdb_assets"] = scene.build_assets(manifest, collections, materials, path_context)
    report["nucleosome_core_material_emphasis"] = emphasize_nucleosome_core(materials)
    report["attachments"] = attachment_table(report["pdb_assets"])
    report["polymerase_rna_origin"] = polymerase_rna_origin_report(report["pdb_assets"], mrna["path"])
    report["arrangement_binding_checks"] = arrangement_binding_checks(report["pdb_assets"], mrna["path"])
    camera_names = add_experiment_cameras(manifest, mrna["path"])
    report["render_label_policy"] = {
        "hidden_non_scale_label_count": hide_non_scale_labels_for_render(),
        "scale_bar_labels_visible": True,
    }
    scene.validate_scene(report)

    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_PATH))
    bpy.context.scene.render.filepath = str(PREVIEW_PATH)
    bpy.ops.render.render(write_still=True)
    for key, path in DETAIL_PREVIEWS.items():
        camera = bpy.data.objects.get(camera_names[key])
        if camera:
            bpy.context.scene.camera = camera
            bpy.context.scene.render.filepath = str(path)
            bpy.ops.render.render(write_still=True)
    print(f"Wrote {BLEND_PATH}")
    print(f"Wrote {PREVIEW_PATH}")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
