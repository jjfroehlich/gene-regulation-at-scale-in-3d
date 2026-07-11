#!/usr/bin/env python3
"""Build comparison renders for elongated and compact RNA structure variants."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import bpy
from mathutils import Vector


ROOT = Path(__file__).resolve().parents[3]
EXPERIMENT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = EXPERIMENT_ROOT / "outputs"
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import blender_nucleic_meshes as nucleic_meshes  # noqa: E402
import scene_core as base  # noqa: E402
import surface_assets as scene  # noqa: E402
import rna_variant_geometry as variant_geometry  # noqa: E402


BLEND_PATH = OUTPUT_DIR / "rna_structure_variants.blend"
REPORT_PATH = OUTPUT_DIR / "rna_structure_variants_report.json"
ELONGATED_PATH = OUTPUT_DIR / "rna_structure_variants_elongated.png"
COMPACT_PATH = OUTPUT_DIR / "rna_structure_variants_compact.png"
DETAIL_DIR = OUTPUT_DIR / "details"


def material(name: str, color: tuple[float, float, float, float]) -> bpy.types.Material:
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = color
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf:
        bsdf.inputs["Base Color"].default_value = color
        bsdf.inputs["Roughness"].default_value = 0.48
        bsdf.inputs["Specular IOR Level"].default_value = 0.22
    return mat


def configure_render() -> None:
    scene.configure_scene()
    render = bpy.context.scene.render
    render.engine = "BLENDER_EEVEE"
    render.resolution_x = 1800
    render.resolution_y = 1000
    render.resolution_percentage = 100
    render.image_settings.file_format = "PNG"
    render.film_transparent = False
    bpy.context.scene.world.color = (0.93, 0.95, 0.965)
    world = bpy.context.scene.world
    world.use_nodes = True
    world.node_tree.nodes["Background"].inputs["Color"].default_value = (0.93, 0.95, 0.965, 1.0)
    world.node_tree.nodes["Background"].inputs["Strength"].default_value = 0.72


def add_lights() -> None:
    for name, location, energy, size in (
        ("Key", (-80.0, -70.0, 130.0), 1600.0, 90.0),
        ("Fill", (90.0, 40.0, 90.0), 900.0, 70.0),
    ):
        data = bpy.data.lights.new(name, "AREA")
        data.energy = energy
        data.shape = "DISK"
        data.size = size
        obj = bpy.data.objects.new(name, data)
        obj.location = location
        bpy.context.scene.collection.objects.link(obj)


def create_camera(name: str, center: tuple[float, float], ortho_scale: float) -> bpy.types.Object:
    data = bpy.data.cameras.new(name)
    data.type = "ORTHO"
    data.ortho_scale = ortho_scale
    camera = bpy.data.objects.new(name, data)
    camera.location = (center[0], center[1], 180.0)
    camera.rotation_euler = (0.0, 0.0, 0.0)
    bpy.context.scene.collection.objects.link(camera)
    return camera


def translated(points, offset: Vector) -> list[tuple[float, float, float]]:
    return [(point.x + offset.x, point.y + offset.y, point.z + offset.z) for point in points]


def surface_tube(name, points, offset, radius, mode, mat, collection):
    """A smooth connected tube whose centerline is unchanged between surface trials."""
    pts = [Vector(p) + offset for p in points]
    sides, vertices, faces = 16, [], []
    for i, point in enumerate(pts):
        tangent = (pts[min(i + 1, len(pts) - 1)] - pts[max(i - 1, 0)]).normalized()
        normal = tangent.cross(Vector((0, 0, 1)))
        if normal.length < 1e-5:
            normal = tangent.cross(Vector((0, 1, 0)))
        normal.normalize(); binormal = tangent.cross(normal).normalized()
        for side in range(sides):
            angle = math.tau * side / sides
            scale = 1.0
            if mode == "soft_molecular":
                scale += 0.075 * math.sin(math.tau * i / 41.0) + 0.035 * math.sin(math.tau * i / 13.0 + 0.8)
                scale += 0.018 * math.cos(3 * angle + i * 0.09)
            elif mode == "twisted_groove":
                scale += 0.085 * math.cos(2 * (angle - i * 0.055)) + 0.018 * math.sin(math.tau * i / 19.0)
            vertices.append(tuple(point + radius * scale * (math.cos(angle) * normal + math.sin(angle) * binormal)))
    for i in range(len(pts) - 1):
        for side in range(sides):
            a = i * sides + side; b = i * sides + (side + 1) % sides
            c = (i + 1) * sides + (side + 1) % sides; d = (i + 1) * sides + side
            faces.append((a, b, c, d))
    obj, _ = nucleic_meshes.create_mesh_object(name, vertices, faces, mat, collection, "rna_surface_trial", polish=False)
    for poly in obj.data.polygons:
        poly.use_smooth = True
    return obj


def protein_lobe(name, center, scale, rotation_deg, mat, collection):
    bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=4, radius=1.0, location=center)
    obj = bpy.context.object
    obj.name = name
    for owner in list(obj.users_collection):
        owner.objects.unlink(obj)
    collection.objects.link(obj)
    obj.scale = scale
    obj.rotation_euler[2] = math.radians(rotation_deg)
    obj.data.materials.append(mat)
    for poly in obj.data.polygons:
        poly.use_smooth = True
    return obj


def add_variant(variant: dict, offset: Vector, materials: dict, parent: bpy.types.Collection) -> dict:
    collection = bpy.data.collections.new(variant["key"])
    parent.children.link(collection)
    objects = []
    bbox = variant["report"]["bbox"]
    bbox_center = Vector(
        tuple((low + high) * 0.5 for low, high in zip(bbox["min_mm"], bbox["max_mm"]))
    )
    geometry_offset = offset - bbox_center
    for segment_model in variant["segments"]:
        segment = segment_model["segment"]
        color_key = segment.get("color", "orange")
        if variant["group"] == "elongated":
            obj = surface_tube(f"{variant['key']} {segment['name']}", segment_model["points"], geometry_offset,
                               variant["report"]["tube_radius_mm"], variant["report"]["surface_mode"],
                               materials[color_key], collection)
        else:
            obj = base.create_curve(f"{variant['key']} {segment['name']}",
                translated([Vector(point) for point in segment_model["points"]], geometry_offset),
                0.13, materials[color_key], collection, resolution=2)
        obj["variant"] = variant["key"]
        objects.append(obj)

    pair_segments = []
    pair_stride = 1 if variant["group"] == "elongated" else 2
    for left, right in variant["base_pairs"][::pair_stride]:
        pair_segments.append((Vector(left) + geometry_offset, Vector(right) + geometry_offset))
    pair_radius = 0.045 if variant["group"] == "elongated" else 0.032
    vertices, faces = nucleic_meshes.cylinder_segments_mesh(pair_segments, pair_radius, 7)
    if pair_segments:
        pair_obj, _ = nucleic_meshes.create_mesh_object(f"{variant['key']} base pairs", vertices, faces,
            materials["base_pairs"], collection, "rna_variant_base_pair_rungs", polish=False)
        pair_obj["variant"] = variant["key"]
        objects.append(pair_obj)

    for index, lobe in enumerate(variant.get("protein_lobes", [])):
        x, y, z, sx, sy, sz, angle = lobe
        objects.append(protein_lobe(f"{variant['key']} protein lobe {index + 1}",
            offset + Vector((x, y, z + 0.15)), (sx, sy, sz), angle,
            materials["protein_a" if index % 2 == 0 else "protein_b"], collection))

    label_y = offset.y + bbox["bbox_mm"][1] * 0.5 + 7.0
    if variant["group"] == "elongated":
        metrics = f"{variant['title']}\n{variant['report']['surface_mode'].replace('_', ' ')} | no explicit stem loops | 222.24 mm"
    else:
        metrics = (f"{variant['title']}\n{variant['report']['stem_count']} stems | "
                   f"{variant['report']['paired_fraction']:.0%} paired | {variant['report']['protein_lobe_count']} protein lobes")
    label = base.create_text(
        f"label {variant['key']}",
        metrics,
        (offset.x, label_y, 8.0),
        0.88 if variant["group"] == "compact" else 1.45,
        materials["text"],
        collection,
    )
    objects.append(label)
    return {
        "key": variant["key"],
        "collection": collection.name,
        "objects": [obj.name for obj in objects],
        "panel_center_mm": list(offset),
        "geometry_recentering_mm": list(-bbox_center),
        "rendered_base_pair_rungs": len(pair_segments),
        "surface_detail_center_mm": list(Vector(variant["points"][len(variant["points"]) // 3]) + geometry_offset),
    }


def set_group_visibility(group: str, variants: list[dict], object_rows: list[dict]) -> None:
    group_by_key = {variant["key"]: variant["group"] for variant in variants}
    for row in object_rows:
        visible = group_by_key[row["key"]] == group
        collection = bpy.data.collections[row["collection"]]
        collection.hide_render = not visible
        collection.hide_viewport = not visible


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    DETAIL_DIR.mkdir(parents=True, exist_ok=True)
    scene.clean_scene()
    configure_render()
    add_lights()
    manifest = variant_geometry.load_manifest()
    variants = variant_geometry.build_all_variants(manifest)

    materials = {
        "olive": material("RNA 5UTR teal", (0.10, 0.62, 0.58, 1.0)),
        "orange": material("RNA coding coral", (0.95, 0.40, 0.14, 1.0)),
        "yellow_olive": material("RNA 3UTR gold", (0.82, 0.68, 0.20, 1.0)),
        "base_pairs": material("RNA base pairs", (0.28, 0.36, 0.48, 1.0)),
        "protein_a": material("RNP protein slate", (0.34, 0.43, 0.56, 1.0)),
        "protein_b": material("RNP protein blue", (0.48, 0.58, 0.70, 1.0)),
        "text": material("Text", (0.08, 0.09, 0.11, 1.0)),
    }
    parent = bpy.data.collections.new("RNA structure variants")
    bpy.context.scene.collection.children.link(parent)

    elongated_offsets = [Vector((-78.0, 0.0, 0.0)), Vector((0.0, 0.0, 0.0)), Vector((78.0, 0.0, 0.0))]
    compact_offsets = [Vector((-40.0, 13.0, 0.0)), Vector((0.0, 13.0, 0.0)), Vector((40.0, 13.0, 0.0)),
                       Vector((-40.0, -15.0, 0.0)), Vector((0.0, -15.0, 0.0)), Vector((40.0, -15.0, 0.0))]
    rows = []
    for variant, offset in zip([v for v in variants if v["group"] == "elongated"], elongated_offsets):
        rows.append(add_variant(variant, offset, materials, parent))
    for variant, offset in zip([v for v in variants if v["group"] == "compact"], compact_offsets):
        rows.append(add_variant(variant, offset, materials, parent))

    elongated_camera = create_camera("Camera elongated", (0.0, 0.0), 220.0)
    compact_camera = create_camera("Camera compact", (0.0, 0.0), 112.0)

    for label in [obj for obj in bpy.data.objects if obj.type == "FONT"]:
        label.rotation_euler = (0.0, 0.0, 0.0)

    set_group_visibility("elongated", variants, rows)
    bpy.context.scene.camera = elongated_camera
    bpy.context.scene.render.filepath = str(ELONGATED_PATH)
    bpy.ops.render.render(write_still=True)

    set_group_visibility("compact", variants, rows)
    bpy.context.scene.camera = compact_camera
    bpy.context.scene.render.filepath = str(COMPACT_PATH)
    bpy.ops.render.render(write_still=True)

    detail_outputs = {}
    render = bpy.context.scene.render
    render.resolution_x = 1200
    render.resolution_y = 1200
    report_by_key = {variant["key"]: variant["report"] for variant in variants}
    for row in rows:
        for other in rows:
            collection = bpy.data.collections[other["collection"]]
            collection.hide_render = other["key"] != row["key"]
            collection.hide_viewport = other["key"] != row["key"]
        collection = bpy.data.collections[row["collection"]]
        labels = [obj for obj in collection.objects if obj.type == "FONT"]
        for label in labels:
            label.hide_render = True
        details = report_by_key[row["key"]]
        bbox = details["bbox"]["bbox_mm"]
        detail_center = row["panel_center_mm"]
        detail_scale = max(bbox[0], bbox[1]) * 1.35
        if details["group"] == "elongated":
            detail_center = row["surface_detail_center_mm"]
            detail_scale = 5.0
        camera = create_camera(
            f"Camera detail {row['key']}",
            (detail_center[0], detail_center[1]),
            detail_scale,
        )
        bpy.context.scene.camera = camera
        detail_path = DETAIL_DIR / f"{row['key']}.png"
        bpy.context.scene.render.filepath = str(detail_path)
        bpy.ops.render.render(write_still=True)
        detail_outputs[row["key"]] = str(detail_path)
        for label in labels:
            label.hide_render = False

    render.resolution_x = 1800
    render.resolution_y = 1000

    for row in rows:
        collection = bpy.data.collections[row["collection"]]
        collection.hide_render = False
        collection.hide_viewport = False

    report = {
        "title": "RNA structure generation variants",
        "canonical_scene_changed": False,
        "manifest": str(variant_geometry.MANIFEST_PATH),
        "total_nt": variant_geometry.TOTAL_NT,
        "target_contour_length_mm": variant_geometry.TOTAL_NT * manifest["units"]["mrna_nt_to_mm"],
        "visual_encoding": {
            "5_utr": "teal",
            "coding_sequence": "coral",
            "3_utr": "gold",
            "base_pair_rungs": "slate; compact candidates render every second modeled pair",
            "protein_surfaces": "blue-slate; schematic architecture inspired by deposited RNP structures",
            "comparison_scale": "constant within each comparison row; compact and elongated rows use separate cameras",
        },
        "preliminary_recommendations": {
            "elongated_clean_reference": "surface_smooth_tube",
            "elongated_surface_candidate": "surface_soft_molecular",
            "elongated_directional_alternative": "surface_twisted_groove",
            "compact_mrna_like_distribution": "rnp_ejc_clamped",
            "compact_rna_scaffold": "rnp_srp_scaffold",
            "compact_particle_silhouette": "rnp_telomerase_bilobal",
            "interpretation": (
                "The elongated row is a controlled surface study with one unchanged centerline. "
                "The compact area compares three RNA-only baselines with distributed clamps, an RNA scaffold with protein saddles, and a bilobal protein-rich particle. "
                "These are architectural analogies, not literal ACTB mRNP compositions."
            ),
        },
        "variants": [variant["report"] for variant in variants],
        "scene_objects": rows,
        "outputs": {
            "blend": str(BLEND_PATH),
            "elongated_comparison": str(ELONGATED_PATH),
            "compact_comparison": str(COMPACT_PATH),
            "details": detail_outputs,
            "report": str(REPORT_PATH),
        },
    }
    failures = []
    for variant in variants:
        details = variant["report"]
        if details["allocation_check_nt"] != variant_geometry.TOTAL_NT:
            failures.append({"variant": variant["key"], "reason": "nucleotide allocation"})
        if not math.isclose(details["measured_length_mm"], report["target_contour_length_mm"], rel_tol=1e-6, abs_tol=1e-6):
            failures.append({"variant": variant["key"], "reason": "contour length", "actual": details["measured_length_mm"]})
        if variant["group"] == "elongated" and (details["stem_count"] != 0 or details["base_pair_bridge_count"] != 0):
            failures.append({"variant": variant["key"], "reason": "surface study contains secondary-structure branches"})
        if variant["group"] == "compact" and (details["stem_count"] <= 0 or details["base_pair_bridge_count"] <= 0):
            failures.append({"variant": variant["key"], "reason": "missing compact RNA stems or base pairs"})
    report["validation_failures"] = failures
    if failures:
        raise RuntimeError(f"RNA variant validation failed: {failures}")

    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_PATH))
    print(f"Wrote {BLEND_PATH}")
    print(f"Wrote {ELONGATED_PATH}")
    print(f"Wrote {COMPACT_PATH}")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
