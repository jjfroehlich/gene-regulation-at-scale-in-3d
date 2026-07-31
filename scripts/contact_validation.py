#!/usr/bin/env python3
"""Shared contact, attachment, and camera helpers for the canonical scene."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import bpy
from mathutils import Vector
from mathutils.kdtree import KDTree


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import blender_nucleic_meshes as direct_nucleic_meshes  # noqa: E402
import scene_core as base  # noqa: E402
import surface_assets as scene  # noqa: E402
import canonical_config as canonical  # noqa: E402
import procedural_nucleic_geometry as nucleic_geometry  # noqa: E402


OUTPUT_DIR = canonical.OUTPUT_DIR
BLEND_PATH = canonical.BLEND_PATH
PREVIEW_PATH = canonical.PREVIEW_PATH
REPORT_PATH = canonical.REPORT_PATH
DETAIL_PREVIEWS = {
    "full_overview": OUTPUT_DIR / "preview_gene_expression_surface_style_full_overview.png",
    "polymerase_rna_start": OUTPUT_DIR / "preview_gene_expression_surface_style_polymerase_rna_start.png",
    "nucleosome_loop": OUTPUT_DIR / "preview_gene_expression_surface_style_nucleosome_loop.png",
    "p53_dna": OUTPUT_DIR / "preview_gene_expression_surface_style_p53_dna.png",
    "ribosome_trna": OUTPUT_DIR / "preview_gene_expression_surface_style_ribosome_trna.png",
    "actin_product": OUTPUT_DIR / "preview_gene_expression_surface_style_actin_product.png",
}

STRICT_GUIDE_OFFSET_MAX_MM = 0.25
BRACKET_GUIDE_OFFSET_MAX_MM = 1.50
FINAL_SURFACE_GAP_MAX_MM = 0.22
FINAL_MESH_GAP_MAX_MM = 0.28


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


def build_current_dna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = direct_nucleic_meshes.build_dna_meshes(manifest, collections, materials)
    report["source_curve"] = add_source_path_curve("current_DNA_source_path", path, collections["DNA"], materials["black"])
    base.create_text(
        "label_DNA_current",
        "ACTB DNA + promoter: 3954 bp",
        (-62.0, -55.0, 0.2),
        2.8,
        materials["black"],
        collections["Labels"],
    )
    return {"path": path, "report": report, "features": {"nucleosome_loop": report.get("nucleosome_loop")}}


def build_current_mrna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = direct_nucleic_meshes.build_mrna_meshes(manifest, collections, materials)
    report["source_curve"] = add_source_path_curve("current_mRNA_source_path", path, collections["mRNA"], materials["black"])
    start = path.points[0]
    base.create_text(
        "label_mRNA_current",
        "actin mRNA from Pol II: 1852 nt",
        (start.x + 4.0, start.y + 8.0, start.z + 0.3),
        2.4,
        materials["black"],
        collections["Labels"],
    )
    return {"path": path, "report": report}


def build_current_compact_mrna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = direct_nucleic_meshes.build_compact_mrna_meshes(manifest, collections, materials)
    compact_center = path_center(path)
    base.create_text(
        "label_compact_mrna_current",
        "compact full-length structured RNA reference",
        (compact_center.x, compact_center.y + 8.0, compact_center.z + 0.4),
        2.2,
        materials["label_grey"],
        collections["Labels"],
    )
    return {"path": path, "report": report}


def create_camera(name: str, location: tuple[float, float, float], target: Vector, ortho_scale: float) -> str:
    camera_data = bpy.data.cameras.new(name)
    camera = bpy.data.objects.new(name, camera_data)
    camera.location = location
    direction = target - camera.location
    camera.rotation_euler = direction.to_track_quat("-Z", "Y").to_euler()
    camera_data.type = "ORTHO"
    camera_data.ortho_scale = ortho_scale
    bpy.context.scene.collection.objects.link(camera)
    return camera.name


def add_shared_cameras(manifest: dict, dna_path, mrna_path, compact_path, collections: dict, materials: dict) -> dict[str, str]:
    base.add_lighting_and_camera(collections, materials)
    mrna_center = path_center(mrna_path)
    compact_center = path_center(compact_path)
    ribosome_point = mrna_path.point_at_length(mrna_path.length * 0.892)
    actin_point = mrna_path.point_at_length(mrna_path.length)
    loop_report = nucleic_geometry.dna_nucleosome_loop_report(manifest) or {}
    loop_center = loop_report.get("center_mm") or [0.0, -9.0, 0.0]

    camera_names = {
        "full_overview": create_camera(
            "Camera_shared_full_overview",
            (-96.0, -142.0, 118.0),
            Vector((12.0, -8.0, 34.0)),
            190.0,
        ),
        "polymerase_rna_start": create_camera(
            "Camera_shared_polymerase_rna_start",
            (mrna_path.points[0].x - 28.0, mrna_path.points[0].y - 42.0, mrna_path.points[0].z + 31.0),
            Vector((mrna_path.points[0].x, mrna_path.points[0].y, mrna_path.points[0].z + 4.0)),
            28.0,
        ),
        "nucleosome_loop": create_camera(
            "Camera_shared_nucleosome_loop",
            (loop_center[0] - 8.0, loop_center[1] - 12.0, (loop_center[2] if len(loop_center) > 2 else 0.0) + 9.0),
            Vector((loop_center[0], loop_center[1], loop_center[2] if len(loop_center) > 2 else 0.0)),
            8.0,
        ),
        "mrna_spiral": create_camera(
            "Camera_shared_mrna_spiral",
            (mrna_center.x - 48.0, mrna_center.y - 72.0, mrna_center.z + 58.0),
            mrna_center,
            92.0,
        ),
        "ribosome_top_translation": create_camera(
            "Camera_shared_ribosome_top_translation",
            (8.0, -40.0, 156.0),
            Vector(
                (
                    (ribosome_point.x + actin_point.x) * 0.5,
                    (ribosome_point.y + actin_point.y) * 0.5,
                    (ribosome_point.z + actin_point.z) * 0.5,
                )
            ),
            26.0,
        ),
        "compact_mrna": create_camera(
            "Camera_shared_compact_mrna",
            (compact_center.x - 12.0, compact_center.y - 30.0, compact_center.z + 32.0),
            compact_center,
            24.0,
        ),
    }
    bpy.context.scene.camera = bpy.data.objects[camera_names["full_overview"]]
    return camera_names


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
                "offset_mm": attachment["offset_mm"],
                "offset_local_mm": attachment["offset_local_mm"],
                "roll_deg": attachment["roll_deg"],
                "binding_mode": attachment["binding_mode"],
            }
        )
    return rows


def vector_length(values) -> float:
    return math.sqrt(sum(float(value) * float(value) for value in values or []))


def guide_offset_length_mm(attachment: dict) -> float:
    # World and local offsets are specified in orthonormal frames in millimeters;
    # combine conservatively for the validation metric.
    return math.sqrt(vector_length(attachment.get("offset_mm")) ** 2 + vector_length(attachment.get("offset_local_mm")) ** 2)


def component_center(component: dict) -> Vector:
    min_v = Vector(component["min_mm"])
    max_v = Vector(component["max_mm"])
    return (min_v + max_v) * 0.5


def path_kdtree(path) -> KDTree:
    tree = KDTree(len(path.points))
    for index, point in enumerate(path.points):
        tree.insert(Vector((point.x, point.y, point.z)), index)
    tree.balance()
    return tree


def sampled_component_vertices(component: dict, max_count: int = 8000) -> list[Vector]:
    obj = bpy.data.objects.get(component.get("object", ""))
    if obj is None or obj.type != "MESH":
        return []
    vertices = obj.data.vertices
    if not vertices:
        return []
    stride = max(1, math.ceil(len(vertices) / max_count))
    return [obj.matrix_world @ vertex.co for index, vertex in enumerate(vertices) if index % stride == 0][:max_count]


def mesh_kdtree_from_components(components: list[dict], max_per_component: int = 25000) -> KDTree | None:
    vertices: list[Vector] = []
    for component in components:
        if component.get("object"):
            vertices.extend(sampled_component_vertices(component, max_per_component))
    if not vertices:
        return None
    tree = KDTree(len(vertices))
    for index, vertex in enumerate(vertices):
        tree.insert(vertex, index)
    tree.balance()
    return tree


def visible_contact_components(item: dict) -> list[dict]:
    components = [component for component in item.get("components", []) if not component.get("alignment_guide_only")]
    protein_components = [component for component in components if component.get("component") == "protein"]
    if protein_components:
        return protein_components
    return components


def final_surface_contact_validation(
    asset_reports: list[dict],
    path_context: dict,
    contact_radii: dict[str, float],
    target_reports: dict[str, dict],
) -> dict:
    strict_names = set(canonical.DNA_STRICT_CONTACT_NAMES) | set(canonical.RNA_STRICT_CONTACT_NAMES) | set(canonical.RNA_BRACKET_CONTACT_NAMES)
    trees = {name: path_kdtree(path_context[name]) for name in ("dna", "mrna")}
    mesh_trees = {
        name: mesh_kdtree_from_components(target_reports[name].get("components", []))
        for name in ("dna", "mrna")
    }
    rows = []
    failures = []
    for item in asset_reports:
        if item.get("name") not in strict_names:
            continue
        attachment = item.get("attachment_empty") or {}
        path_name = attachment.get("path")
        if path_name not in trees:
            failures.append({"name": item.get("name"), "reason": "missing_final_contact_path", "path": path_name})
            continue
        components = visible_contact_components(item)
        if not components:
            failures.append({"name": item.get("name"), "reason": "no_visible_contact_components"})
            continue
        tree = trees[path_name]
        min_distance = float("inf")
        nearest_component = None
        nearest_vertex = None
        nearest_path_point = None
        min_mesh_distance = float("inf")
        nearest_mesh_point = None
        sampled_vertices = 0
        for component in components:
            vertices = sampled_component_vertices(component)
            sampled_vertices += len(vertices)
            for vertex in vertices:
                path_point, _index, distance = tree.find(vertex)
                if distance < min_distance:
                    min_distance = distance
                    nearest_component = component.get("object")
                    nearest_vertex = vertex
                    nearest_path_point = Vector(path_point)
                mesh_tree = mesh_trees.get(path_name)
                if mesh_tree is not None:
                    mesh_point, _mesh_index, mesh_distance = mesh_tree.find(vertex)
                    if mesh_distance < min_mesh_distance:
                        min_mesh_distance = mesh_distance
                        nearest_mesh_point = Vector(mesh_point)
        if nearest_vertex is None or nearest_path_point is None:
            failures.append({"name": item.get("name"), "reason": "no_sampled_contact_vertices"})
            continue
        mesh_gap = min_mesh_distance if nearest_mesh_point is not None else None
        radius = contact_radii[path_name]
        surface_gap = min_distance - radius
        row = {
            "name": item["name"],
            "pdb_id": item["pdb_id"],
            "path": path_name,
            "nearest_component": nearest_component,
            "sampled_vertices": sampled_vertices,
            "target_surface_radius_mm": radius,
            "min_centerline_distance_mm": min_distance,
            "surface_gap_mm": surface_gap,
            "max_allowed_surface_gap_mm": FINAL_SURFACE_GAP_MAX_MM,
            "min_target_mesh_distance_mm": mesh_gap,
            "max_allowed_target_mesh_distance_mm": FINAL_MESH_GAP_MAX_MM,
            "nearest_asset_vertex_mm": [nearest_vertex.x, nearest_vertex.y, nearest_vertex.z],
            "nearest_path_point_mm": [nearest_path_point.x, nearest_path_point.y, nearest_path_point.z],
            "nearest_target_mesh_vertex_mm": [nearest_mesh_point.x, nearest_mesh_point.y, nearest_mesh_point.z] if nearest_mesh_point else None,
        }
        rows.append(row)
        if surface_gap > FINAL_SURFACE_GAP_MAX_MM:
            failures.append(
                {
                    "name": item["name"],
                    "reason": "visible_surface_not_binding_procedural_nucleic",
                    "path": path_name,
                    "surface_gap_mm": surface_gap,
                    "max_mm": FINAL_SURFACE_GAP_MAX_MM,
                }
            )
        if mesh_gap is not None and mesh_gap > FINAL_MESH_GAP_MAX_MM:
            failures.append(
                {
                    "name": item["name"],
                    "reason": "visible_surface_not_touching_procedural_nucleic_mesh",
                    "path": path_name,
                    "mesh_gap_mm": mesh_gap,
                    "max_mm": FINAL_MESH_GAP_MAX_MM,
                }
            )
    return {
        "target_surface_radii_mm": contact_radii,
        "max_allowed_surface_gap_mm": FINAL_SURFACE_GAP_MAX_MM,
        "max_allowed_target_mesh_distance_mm": FINAL_MESH_GAP_MAX_MM,
        "rows": rows,
        "failures": failures,
    }


def contact_validation(asset_reports: list[dict]) -> dict:
    by_name = {item["name"]: item for item in asset_reports}
    strict_sets = {
        "dna": sorted(canonical.DNA_STRICT_CONTACT_NAMES),
        "mrna": sorted(canonical.RNA_STRICT_CONTACT_NAMES),
        "mrna_bracket": sorted(canonical.RNA_BRACKET_CONTACT_NAMES),
    }
    expected_path = {"dna": "dna", "mrna": "mrna", "mrna_bracket": "mrna"}
    max_offset = {
        "dna": STRICT_GUIDE_OFFSET_MAX_MM,
        "mrna": STRICT_GUIDE_OFFSET_MAX_MM,
        "mrna_bracket": BRACKET_GUIDE_OFFSET_MAX_MM,
    }
    rows = []
    failures = []
    for role, names in strict_sets.items():
        for name in names:
            item = by_name.get(name)
            if item is None:
                failures.append({"name": name, "reason": "missing_asset_report"})
                continue
            attachment = item.get("attachment_empty")
            if not attachment:
                failures.append({"name": name, "reason": "missing_attachment", "role": role})
                continue
            binding_alignment = item.get("binding_alignment") or {}
            offset = guide_offset_length_mm(attachment)
            expected = expected_path[role]
            row = {
                "name": name,
                "pdb_id": item["pdb_id"],
                "role": role,
                "expected_path": expected,
                "actual_path": attachment.get("path"),
                "guide_anchor_offset_mm": offset,
                "max_allowed_offset_mm": max_offset[role],
                "anchor_kind": item.get("anchor_kind"),
                "path_anchor_mode": item.get("path_anchor_mode"),
                "binding_mode": attachment.get("binding_mode"),
                "binding_alignment_mode": binding_alignment.get("mode"),
                "co_crystal_guide_available": binding_alignment.get("co_crystal_guide_available", True),
                "distance_mm": attachment.get("distance_mm"),
                "fraction": attachment.get("fraction"),
            }
            rows.append(row)
            if attachment.get("path") != expected:
                failures.append({"name": name, "reason": "wrong_attachment_path", "expected": expected, "actual": attachment.get("path")})
            protein_only_rna_contact = (
                role == "mrna"
                and attachment.get("binding_mode") == "protein_surface_rna_contact"
                and item.get("anchor_kind") == "protein_surface"
            )
            if item.get("anchor_kind") != "nucleic" and not protein_only_rna_contact:
                failures.append({"name": name, "reason": "non_nucleic_anchor", "anchor_kind": item.get("anchor_kind")})
            if offset > max_offset[role]:
                failures.append({"name": name, "reason": "guide_anchor_offset_too_large", "offset_mm": offset, "max_mm": max_offset[role]})
            alignment_mode = row["binding_alignment_mode"]
            if protein_only_rna_contact:
                if alignment_mode != "protein_surface_contact_to_procedural_rna_path_without_co_crystal_guide":
                    failures.append({"name": name, "reason": "missing_protein_surface_rna_contact_alignment", "alignment_mode": alignment_mode})
            elif role in {"mrna", "mrna_bracket"} and alignment_mode != "co_crystal_rna_interface_surface_to_procedural_rna_surface":
                failures.append({"name": name, "reason": "missing_rna_binding_side_alignment", "alignment_mode": alignment_mode})
            if role == "dna" and name != "Nucleosome" and alignment_mode != "co_crystal_dna_interface_surface_to_procedural_dna_surface":
                failures.append({"name": name, "reason": "missing_dna_binding_side_alignment", "alignment_mode": alignment_mode})

    display_rows = []
    for name in sorted(canonical.DISPLAY_EXEMPT_NAMES):
        item = by_name.get(name)
        if not item:
            continue
        attachment = item.get("attachment_empty")
        display_rows.append(
            {
                "name": name,
                "pdb_id": item["pdb_id"],
                "role": "display_exempt",
                "actual_path": attachment.get("path") if attachment else None,
                "guide_anchor_offset_mm": guide_offset_length_mm(attachment) if attachment else None,
                "strict_contact_required": False,
            }
        )

    dna_paths = {row["name"]: row["actual_path"] for row in rows if row["role"] == "dna"}
    rna_paths = {row["name"]: row["actual_path"] for row in rows if row["role"] in {"mrna", "mrna_bracket"}}
    rna_display_paths = {row["name"]: row["actual_path"] for row in display_rows if row["actual_path"] == "mrna"}
    nucleosome_check = None
    nucleosome = by_name.get("Nucleosome")
    if nucleosome:
        protein_component = next((component for component in nucleosome.get("components", []) if component.get("component") == "protein"), None)
        target_center = ((nucleosome.get("binding_alignment") or {}).get("target_center_mm"))
        if protein_component and target_center:
            core_center = component_center(protein_component)
            target = Vector(target_center)
            distance = (core_center - target).length
            nucleosome_check = {
                "core_center_mm": [core_center.x, core_center.y, core_center.z],
                "wrapped_dna_center_mm": [target.x, target.y, target.z],
                "core_to_wrap_center_distance_mm": distance,
                "max_allowed_distance_mm": 1.2,
            }
            if distance > 1.2:
                failures.append({"name": "Nucleosome", "reason": "core_not_centered_in_dna_wrap", "distance_mm": distance, "max_mm": 1.2})

    return {
        "policy": "strict procedural nucleic contact; co-crystal guide alignment where available, protein-surface contact for protein-only RBP structures",
        "strict_offset_max_mm": STRICT_GUIDE_OFFSET_MAX_MM,
        "ribosome_bracket_offset_max_mm": BRACKET_GUIDE_OFFSET_MAX_MM,
        "nucleosome_core_in_wrap": nucleosome_check,
        "dna_binding_assets_on_dna_path": dna_paths,
        "rna_binding_assets_on_mrna_path": rna_paths,
        "rna_display_assets_on_mrna_path": rna_display_paths,
        "dna_binding_all_on_dna": all(path == "dna" for path in dna_paths.values()),
        "rna_binding_all_on_mrna": all(path == "mrna" for path in rna_paths.values()),
        "strict_contact_rows": rows,
        "display_exempt_rows": display_rows,
        "strict_contact_failures": failures,
    }


def polymerase_rna_origin_report(asset_reports: list[dict], mrna_path) -> dict | None:
    polymerase = next((item for item in asset_reports if item.get("name") == "RNA polymerase II elongation complex"), None)
    if not polymerase:
        return None
    attachment = polymerase.get("attachment_empty") or {}
    polymerase_location = Vector(polymerase["location_mm"])
    tss_values = attachment.get("path_point_mm") or attachment.get("anchor_point_mm") or polymerase["location_mm"]
    polymerase_tss = Vector(tss_values)
    start = mrna_path.points[0]
    rna_start = Vector((start.x, start.y, start.z))
    return {
        "polymerase_location_mm": [polymerase_location.x, polymerase_location.y, polymerase_location.z],
        "polymerase_tss_anchor_mm": [polymerase_tss.x, polymerase_tss.y, polymerase_tss.z],
        "rna_start_mm": [rna_start.x, rna_start.y, rna_start.z],
        "distance_mm": (polymerase_tss - rna_start).length,
        "protein_contact_to_rna_start_distance_mm": (polymerase_location - rna_start).length,
    }


def hide_non_scale_labels_for_render() -> int:
    hidden = 0
    for obj in bpy.data.objects:
        if obj.type == "FONT" and not obj.name.startswith("label_scale_"):
            obj.hide_render = True
            hidden += 1
    return hidden


def validate_contact_report(report: dict) -> None:
    failures = report.get("contact_validation", {}).get("strict_contact_failures", [])
    origin_distance = (report.get("polymerase_rna_origin") or {}).get("distance_mm")
    if origin_distance is None:
        failures.append({"name": "RNA polymerase II elongation complex", "reason": "missing_polymerase_rna_origin"})
    elif origin_distance > 0.35:
        failures.append(
            {
                "name": "RNA polymerase II elongation complex",
                "reason": "polymerase_not_at_rna_origin",
                "distance_mm": origin_distance,
                "max_mm": 0.35,
            }
        )
    if failures:
        raise RuntimeError(f"Canonical contact validation failed: {failures}")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if not canonical.MANIFEST_PATH.exists():
        canonical.write_manifest()
    manifest = json.loads(canonical.MANIFEST_PATH.read_text(encoding="utf-8"))

    scene.clean_scene()
    scene.configure_scene()
    collections = base.make_collections(manifest["collections"])
    materials = base.make_materials()
    scene.soften_materials(materials)

    dna = build_current_dna(manifest, collections, materials)
    mrna = build_current_mrna(manifest, collections, materials)
    compact_mrna = build_current_compact_mrna(manifest, collections, materials)
    mrna_settings = manifest.get("procedural_nucleic_acids", {}).get("mrna", {})
    mrna_base_radii = mrna_settings.get("rna_base_ellipsoid_radii_mm", [0.135, 0.078, 0.055])
    contact_radii = {
        "dna": float(dna["report"].get("estimated_envelope_diameter_mm", 0.88)) * 0.5,
        "mrna": max(
            float(mrna["report"].get("tube_radius_mm", 0.165)) + 0.08,
            float(mrna_settings.get("base_offset_mm", 0.16)) + float(max(mrna_base_radii)),
        ),
    }
    path_context = {
        "dna": dna["path"],
        "dna_features": dna["features"],
        "dna_anchor_mode": "closest_xy",
        "dna_contact_radius_mm": contact_radii["dna"],
        "mrna": mrna["path"],
        "mrna_contact_radius_mm": contact_radii["mrna"],
    }

    report = {
        "title": manifest["title"],
        "kind": canonical.REPORT_KIND,
        "units": manifest["units"],
        "source_manifest": str(canonical.MANIFEST_PATH),
        "surface_asset_dir": str(scene.REDUCED_SURFACE_DIR),
        "nucleic_acid_pipeline": "direct_blender_meshes_for_procedural_dna_mrna_and_compact_mrna",
        "reduction": scene.load_reduction_summary(),
        "outputs": {
            "blend": str(BLEND_PATH),
            "preview": str(PREVIEW_PATH),
            "report": str(REPORT_PATH),
            "detail_previews": {name: str(path) for name, path in DETAIL_PREVIEWS.items()},
        },
        "layout_intent": manifest.get("layout_intent", {}),
    }
    report["dna"] = dna["report"]
    report["mrna"] = mrna["report"]
    report["compact_mrna"] = compact_mrna["report"]
    report["scale_bars"] = scene.add_scale_bars(manifest, collections, materials)
    report["pdb_assets"] = scene.build_assets(manifest, collections, materials, path_context)
    report["attachments"] = attachment_table(report["pdb_assets"])
    report["polymerase_rna_origin"] = polymerase_rna_origin_report(report["pdb_assets"], mrna["path"])
    report["contact_validation"] = contact_validation(report["pdb_assets"])
    final_contact = final_surface_contact_validation(
        report["pdb_assets"],
        path_context,
        contact_radii,
        {"dna": dna["report"], "mrna": mrna["report"]},
    )
    report["contact_validation"]["final_surface_contact"] = final_contact
    report["contact_validation"]["strict_contact_failures"].extend(final_contact["failures"])
    camera_names = add_shared_cameras(manifest, dna["path"], mrna["path"], compact_mrna["path"], collections, materials)
    report["render_label_policy"] = {
        "hidden_non_scale_label_count": hide_non_scale_labels_for_render(),
        "scale_bar_labels_visible": True,
    }
    scene.validate_scene(report)
    validate_contact_report(report)

    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_PATH))

    primary_camera = bpy.data.objects[camera_names["full_overview"]]
    bpy.context.scene.camera = primary_camera
    bpy.context.scene.render.filepath = str(PREVIEW_PATH)
    bpy.ops.render.render(write_still=True)
    for key, path in DETAIL_PREVIEWS.items():
        camera = bpy.data.objects.get(camera_names[key])
        if camera:
            bpy.context.scene.camera = camera
            bpy.context.scene.render.filepath = str(path)
            bpy.ops.render.render(write_still=True)
    bpy.context.scene.camera = primary_camera
    print(f"Wrote {BLEND_PATH}")
    print(f"Wrote {PREVIEW_PATH}")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    raise RuntimeError("contact_validation.py is a shared module; run run_canonical_workflow.ps1")
