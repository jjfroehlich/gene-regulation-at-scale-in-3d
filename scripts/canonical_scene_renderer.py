#!/usr/bin/env python3
"""Shared renderer for the canonical gene-expression scene."""

from __future__ import annotations

import copy
import json
import math
import os
import shutil
import sys
from pathlib import Path

import bpy
from mathutils import Vector


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import blender_nucleic_meshes as direct_nucleic_meshes  # noqa: E402
import scene_core as base  # noqa: E402
import surface_assets as scene  # noqa: E402
import contact_validation as contact_helpers  # noqa: E402
import canonical_config as canonical  # noqa: E402
import procedural_nucleic_geometry as nucleic_geometry  # noqa: E402


OUTPUT_DIR = canonical.OUTPUT_DIR
BLEND_PATH = canonical.BLEND_PATH
PREVIEW_PATH = canonical.PREVIEW_PATH
REPORT_PATH = canonical.REPORT_PATH
DETAIL_PREVIEWS = {
    "full_overview": OUTPUT_DIR / "preview_gene_expression_surface_style_canonical_full_overview.png",
    "p53_dna": OUTPUT_DIR / "preview_gene_expression_surface_style_canonical_p53_dna.png",
    "polymerase_rna_start": OUTPUT_DIR / "preview_gene_expression_surface_style_canonical_polymerase_rna_start.png",
    "nucleosome_loop": OUTPUT_DIR / "preview_gene_expression_surface_style_canonical_nucleosome_loop.png",
    "ribosome_trna": OUTPUT_DIR / "preview_gene_expression_surface_style_canonical_ribosome_trna.png",
    "actin_product": OUTPUT_DIR / "preview_gene_expression_surface_style_canonical_actin_product.png",
    "cas9_dna": OUTPUT_DIR / "preview_gene_expression_surface_style_canonical_cas9_dna.png",
}
DETAIL_SHARED_ORTHO_SCALE_MM = 42.0
PROTEIN_AA_CONTOUR_NM = canonical.PROTEIN_AA_CONTOUR_NM
NANOMETER_SCALE_BAR_NM = 10.0
Canonical_MATERIAL_COLORS = {
    # DNA keeps distinct promoter/exon/intron semantics.
    "dna_promoter_blue": (0.12, 0.36, 0.78, 1.0),
    "dna_exon_orange": (0.95, 0.62, 0.16, 1.0),
    "dna_intron_olive": (0.42, 0.55, 0.34, 1.0),
    "dna_orange": (0.88, 0.56, 0.18, 1.0),
    "dna_dark": (0.18, 0.30, 0.52, 1.0),
    "dna_guide_amber": (0.91, 0.57, 0.19, 1.0),
    "guide_orange": (0.91, 0.57, 0.19, 1.0),
    # mRNA keeps separate 5' UTR, coding sequence, and 3' UTR colors.
    "olive": (0.16, 0.62, 0.56, 1.0),
    "orange": (0.95, 0.45, 0.16, 1.0),
    "yellow_olive": (0.82, 0.68, 0.26, 1.0),
    "rna_gold": (0.95, 0.45, 0.16, 1.0),
    "rna_red": (0.95, 0.45, 0.16, 1.0),
    "ribosome_red": (0.95, 0.45, 0.16, 1.0),
    "rna_guide_orange": (0.95, 0.45, 0.16, 1.0),
    # Context proteins are deliberately uniform; only actin product is highlighted.
    "protein_uniform_slate": (0.36, 0.44, 0.54, 1.0),
    "protein_product_coral": (0.94, 0.20, 0.15, 1.0),
    "actin_blue": (0.94, 0.20, 0.15, 1.0),
    "pol_green": (0.36, 0.44, 0.54, 1.0),
    "histone_green": (0.36, 0.44, 0.54, 1.0),
    "cas9_blue": (0.36, 0.44, 0.54, 1.0),
    "ribosome_blue": (0.36, 0.44, 0.54, 1.0),
    "tf_cyan": (0.36, 0.44, 0.54, 1.0),
    "tf_lavender": (0.36, 0.44, 0.54, 1.0),
    "tf_purple": (0.36, 0.44, 0.54, 1.0),
    "ago_pink": (0.36, 0.44, 0.54, 1.0),
    "rbp_purple": (0.36, 0.44, 0.54, 1.0),
    "mcp_blue": (0.36, 0.44, 0.54, 1.0),
    "rfp_red": (0.36, 0.44, 0.54, 1.0),
    "label_grey": (0.14, 0.15, 0.16, 1.0),
    "scale_grey": (0.24, 0.25, 0.26, 1.0),
    "black": (0.02, 0.022, 0.025, 1.0),
}
Canonical_WORLD_COLOR = (0.94, 0.955, 0.97, 1.0)
Canonical_BACKDROP_TOP_LEFT = (0.86, 0.925, 0.975, 1.0)
Canonical_BACKDROP_TOP_RIGHT = (0.96, 0.985, 0.995, 1.0)
Canonical_BACKDROP_BOTTOM_LEFT = (0.985, 0.945, 0.82, 1.0)
Canonical_BACKDROP_BOTTOM_RIGHT = (1.0, 0.885, 0.74, 1.0)
Canonical_LABEL_MATERIALS = {"black", "label_grey", "scale_grey"}
Canonical_CYCLES_SAMPLES = 64
OVERVIEW_LABEL_POSITIONS = {
    "label_Transcription factor 4": (0.36, 0.39),
    "label_Cas9": (0.46, 0.38),
    "label_Transcription factor 1": (0.48, 0.45),
    "label_MS2 coat protein MCP": (0.49, 0.30),
    "label_mCherry/RFP tag": (0.61, 0.29),
    "label_p53 tetramer bound to DNA": (0.61, 0.48),
    "label_Pumilio RBP": (0.52, 0.50),
    "label_RNA polymerase II elongation complex": (0.66, 0.33),
    "label_Nucleosome": (0.77, 0.50),
    "label_Transcription factor 3": (0.79, 0.44),
    "label_Argonaute": (0.73, 0.55),
    "label_Poly(A)-binding RBP": (0.60, 0.69),
    "label_HuR-like RBP": (0.73, 0.67),
    "label_Ribosome large subunit": (0.64, 0.73),
    "label_Ribosome small subunit": (0.71, 0.73),
    "label_Standalone tRNA": (0.78, 0.71),
}
PRIMARY_CALLOUTS = {
    "label_DNA_canonical": {"text": "Actb promoter + gene DNA 3954 bp", "view_position": (0.825, 0.29), "span": (0.16, 0.42)},
    "label_mRNA_canonical": {"text": "Actb mRNA 1852 nt", "view_position": (0.825, 0.63), "span": (0.44, 0.82)},
    "label_ACTB_primary_canonical": {"text": "ACTB protein 375 aa", "view_position": (0.825, 0.92), "span": (0.87, 0.97)},
}
COMPACT_CALLOUT = {"object": "label_compact_mrna_canonical", "text": "mRNA compact", "offset": (0.014, 0.014)}
DETAIL_TITLES = {
    "p53_dna": "p53 tetramer + DNA (3TS8)",
    "nucleosome_loop": "nucleosome core + wrapped DNA (1AOI)",
    "polymerase_rna_start": "RNA polymerase II + nascent RNA (2E2I)",
    "ribosome_trna": "ribosome + tRNA",
    "actin_product": "ACTB protein (1J6Z)",
    "cas9_dna": "Cas9 + guide/target DNA (4UN3)",
}


def add_source_path_curve(name: str, path, collection, material) -> dict:
    points = [(point.x, point.y, point.z + 0.16) for point in path.points]
    obj = base.create_curve(name, points, 0.028, material, collection, resolution=1)
    obj.hide_viewport = True
    obj.hide_render = True
    try:
        obj.hide_set(True)
    except RuntimeError:
        pass
    obj["source_path_for_rebake"] = True
    obj["sampled_points"] = len(points)
    obj["curve_role"] = "editable_source_path"
    return {
        "object": obj.name,
        "points": len(points),
        "hide_viewport": bool(obj.hide_viewport),
        "hide_render": bool(obj.hide_render),
    }


def build_canonical_dna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = direct_nucleic_meshes.build_dna_meshes(manifest, collections, materials)
    report["source_curve"] = add_source_path_curve("Canonical_DNA_source_path", path, collections["DNA"], materials["black"])
    start = path.points[0]
    base.create_text(
        "label_DNA_canonical",
        "Actb promoter + gene DNA 3954 bp",
        (start.x + 4.0, start.y - 3.0, start.z + 0.25),
        1.55,
        materials["black"],
        collections["Labels"],
    )
    return {"path": path, "report": report, "features": {"nucleosome_loop": report.get("nucleosome_loop")}}


def build_canonical_mrna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = direct_nucleic_meshes.build_mrna_meshes(manifest, collections, materials)
    report["source_curve"] = add_source_path_curve("Canonical_mRNA_source_path", path, collections["mRNA"], materials["black"])
    start = path.points[0]
    base.create_text(
        "label_mRNA_canonical",
        "Actb mRNA 1852 nt",
        (start.x + 3.0, start.y + 3.0, start.z + 0.3),
        1.55,
        materials["black"],
        collections["Labels"],
    )
    return {"path": path, "report": report}


def compact_positioned_manifest(manifest: dict, mrna_path) -> tuple[dict, dict]:
    adjusted = copy.deepcopy(manifest)
    mrna_end = Vector((mrna_path.points[-1].x, mrna_path.points[-1].y, mrna_path.points[-1].z))
    actin = canonical.asset_by_name(adjusted, "Actin protein")
    actin_location = Vector(actin["location_mm"])
    center = mrna_end.lerp(actin_location, 0.55)
    adjusted["procedural_nucleic_acids"]["mrna"]["compact_center_mm"] = [center.x, center.y, center.z]
    return adjusted, {
        "mode": "between_mrna_end_and_actin",
        "mrna_end_mm": [mrna_end.x, mrna_end.y, mrna_end.z],
        "actin_location_mm": [actin_location.x, actin_location.y, actin_location.z],
        "computed_center_mm": [center.x, center.y, center.z],
        "interpolation_from_mrna_end_to_actin": 0.55,
    }


def place_actin_from_mrna_endpoint(manifest: dict, mrna_path) -> dict:
    endpoint = Vector((mrna_path.points[-1].x, mrna_path.points[-1].y, mrna_path.points[-1].z))
    offset = Vector((-6.4, 0.6, 27.0))
    location = endpoint + offset
    actin = canonical.asset_by_name(manifest, "Actin protein")
    actin.pop("path_anchor", None)
    actin["location_mm"] = [location.x, location.y, location.z]
    actin["label_text"] = "actin (1J6Z)"
    return {
        "mode": "above_current_structured_mrna_endpoint",
        "mrna_end_mm": [endpoint.x, endpoint.y, endpoint.z],
        "offset_mm": [offset.x, offset.y, offset.z],
        "actin_location_mm": [location.x, location.y, location.z],
    }


def build_canonical_compact_mrna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = direct_nucleic_meshes.build_compact_mrna_meshes(manifest, collections, materials)
    compact_center = contact_helpers.path_center(path)
    base.create_text(
        "label_compact_mrna_canonical",
        "mRNA compact",
        (compact_center.x, compact_center.y + 3.2, compact_center.z + 0.35),
        1.45,
        materials["label_grey"],
        collections["Labels"],
    )
    return {"path": path, "report": report, "center_mm": [compact_center.x, compact_center.y, compact_center.z]}


def add_canonical_scale_bars(manifest: dict, collections: dict, materials: dict) -> dict:
    units = manifest["units"]
    nm_to_mm = float(units["nm_to_mm"])
    protein_aa_contour_nm = float(units.get("protein_aa_contour_nm", PROTEIN_AA_CONTOUR_NM))
    bars = [
        {
            "name": "scale_dna_100_bp",
            "label": "DNA 100 bp",
            "length_nm": 100.0 * float(units["dna_bp_rise_nm"]),
            "origin": (-82.0, -2.0, 0.0),
            "quantity": 100,
            "unit": "bp",
            "basis": "B-form DNA axial rise",
        },
        {
            "name": "scale_rna_100_nt",
            "label": "RNA 100 nt",
            "length_nm": 100.0 * float(units["mrna_nt_contour_nm"]),
            "origin": (-82.0, -7.0, 0.0),
            "quantity": 100,
            "unit": "nt",
            "basis": "single-stranded RNA contour length",
        },
        {
            "name": "scale_protein_33_aa",
            "label": "protein 33 aa",
            "length_nm": 33.0 * protein_aa_contour_nm,
            "origin": (-82.0, -12.0, 0.0),
            "quantity": 33,
            "unit": "aa",
            "basis": f"translation-equivalent unfolded polypeptide contour length at {protein_aa_contour_nm:g} nm/aa",
        },
        {
            "name": "scale_10_nm",
            "label": "10 nm",
            "length_nm": NANOMETER_SCALE_BAR_NM,
            "origin": (-82.0, -17.0, 0.0),
            "quantity": NANOMETER_SCALE_BAR_NM,
            "unit": "nm",
            "basis": "nanometer reference",
        },
    ]
    reports = {}
    for bar in bars:
        x, y, z = bar["origin"]
        length = bar["length_nm"] * nm_to_mm
        name = bar["name"]
        base.create_curve(name, [(x, y, z), (x + length, y, z)], 0.07, materials["scale_grey"], collections["Scale bars"])
        base.create_curve(f"{name}_left_tick", [(x, y - 0.8, z), (x, y + 0.8, z)], 0.055, materials["scale_grey"], collections["Scale bars"])
        base.create_curve(
            f"{name}_right_tick",
            [(x + length, y - 0.8, z), (x + length, y + 0.8, z)],
            0.055,
            materials["scale_grey"],
            collections["Scale bars"],
        )
        base.create_text(
            f"label_{name}",
            bar["label"],
            (x + length + 2.2, y, z),
            1.55,
            materials["label_grey"],
            collections["Labels"],
            align="LEFT",
        )
        reports[name] = {
            "label": bar["label"],
            "quantity": bar["quantity"],
            "unit": bar["unit"],
            "basis": bar["basis"],
            "length_nm": bar["length_nm"],
            "length_mm": length,
            "origin_mm": [x, y, z],
        }
    return reports


def create_camera(name: str, location: tuple[float, float, float], target: Vector, ortho_scale: float) -> str:
    return contact_helpers.create_camera(name, location, target, ortho_scale)


def set_bsdf_input(bsdf, names: tuple[str, ...], value) -> None:
    for name in names:
        if name in bsdf.inputs:
            bsdf.inputs[name].default_value = value
            return


def set_material_color(material: bpy.types.Material, color: tuple[float, float, float, float]) -> None:
    material.diffuse_color = color
    material.use_nodes = True
    bsdf = material.node_tree.nodes.get("Principled BSDF")
    if bsdf:
        bsdf.inputs["Base Color"].default_value = color
        label_material = material.name.split(".")[0] in Canonical_LABEL_MATERIALS
        bsdf.inputs["Roughness"].default_value = 0.74 if label_material else 0.56
        bsdf.inputs["Metallic"].default_value = 0.0
        set_bsdf_input(bsdf, ("Specular IOR Level", "Specular"), 0.44 if not label_material else 0.22)
        set_bsdf_input(bsdf, ("Coat Weight", "Coat"), 0.08 if not label_material else 0.0)
        set_bsdf_input(bsdf, ("Coat Roughness",), 0.72)
        if "Emission Color" in bsdf.inputs:
            bsdf.inputs["Emission Color"].default_value = color
        elif "Emission" in bsdf.inputs:
            bsdf.inputs["Emission"].default_value = color
        if "Emission Strength" in bsdf.inputs:
            bsdf.inputs["Emission Strength"].default_value = 0.035 if label_material else 0.0


def set_color_management(view_transform: str, look: str | None = None) -> None:
    view_settings = bpy.context.scene.view_settings
    transform_items = view_settings.bl_rna.properties["view_transform"].enum_items
    if view_transform in {item.identifier for item in transform_items}:
        view_settings.view_transform = view_transform
    if look is not None:
        look_items = view_settings.bl_rna.properties["look"].enum_items
        if look in {item.identifier for item in look_items}:
            view_settings.look = look


def configure_canonical_beauty_render() -> dict:
    scene_data = bpy.context.scene
    requested_engine = "CYCLES"
    engine_fallback_reason = None
    try:
        bpy.ops.preferences.addon_enable(module="cycles")
        scene_data.render.engine = requested_engine
    except Exception as exc:  # pragma: no cover - depends on local Blender build.
        scene_data.render.engine = "BLENDER_EEVEE"
        engine_fallback_reason = str(exc)
    scene_data.render.resolution_x = 3200
    scene_data.render.resolution_y = 1920
    scene_data.render.film_transparent = False
    cycles_settings = {}
    if scene_data.render.engine == "CYCLES":
        cycles = scene_data.cycles
        cycles.samples = Canonical_CYCLES_SAMPLES
        cycles.preview_samples = 24
        cycles.max_bounces = 6
        cycles.diffuse_bounces = 3
        cycles.glossy_bounces = 4
        cycles.transparent_max_bounces = 4
        if hasattr(cycles, "use_adaptive_sampling"):
            cycles.use_adaptive_sampling = True
        if hasattr(cycles, "adaptive_threshold"):
            cycles.adaptive_threshold = 0.035
        if hasattr(cycles, "use_denoising"):
            cycles.use_denoising = True
        if hasattr(cycles, "use_fast_gi"):
            cycles.use_fast_gi = True
        cycles_settings = {
            "samples": int(cycles.samples),
            "preview_samples": int(cycles.preview_samples),
            "max_bounces": int(cycles.max_bounces),
            "diffuse_bounces": int(cycles.diffuse_bounces),
            "glossy_bounces": int(cycles.glossy_bounces),
            "denoising": bool(getattr(cycles, "use_denoising", False)),
            "adaptive_sampling": bool(getattr(cycles, "use_adaptive_sampling", False)),
            "adaptive_threshold": float(getattr(cycles, "adaptive_threshold", 0.0)),
        }
    if hasattr(scene_data.eevee, "use_gtao"):
        scene_data.eevee.use_gtao = True
    if hasattr(scene_data.eevee, "gtao_distance"):
        scene_data.eevee.gtao_distance = 5.0
    if hasattr(scene_data.eevee, "gtao_factor"):
        scene_data.eevee.gtao_factor = 0.55
    scene_data.eevee.taa_render_samples = 96
    scene_data.eevee.taa_samples = 32
    scene_data.eevee.use_shadows = True
    scene_data.eevee.shadow_resolution_scale = 1.5
    if hasattr(scene_data.eevee, "use_fast_gi"):
        scene_data.eevee.use_fast_gi = True
    if hasattr(scene_data.eevee, "gi_diffuse_bounces"):
        scene_data.eevee.gi_diffuse_bounces = 3
    set_color_management("AgX", "Medium High Contrast")
    if scene_data.view_settings.view_transform == "Standard":
        scene_data.view_settings.look = "None"
    scene_data.view_settings.exposure = 0.55 if scene_data.render.engine == "CYCLES" else -0.12
    scene_data.view_settings.gamma = 1.0

    world = scene_data.world or bpy.data.worlds.new("World")
    scene_data.world = world
    world.color = Canonical_WORLD_COLOR[:3]
    world.use_nodes = True
    background = world.node_tree.nodes.get("Background")
    if background:
        background.inputs["Color"].default_value = Canonical_WORLD_COLOR
        background.inputs["Strength"].default_value = 0.52 if scene_data.render.engine == "CYCLES" else 0.48
    return {
        "requested_engine": requested_engine,
        "engine": scene_data.render.engine,
        "engine_fallback_reason": engine_fallback_reason,
        "resolution": [scene_data.render.resolution_x, scene_data.render.resolution_y],
        "world_color": list(Canonical_WORLD_COLOR),
        "world_strength": float(background.inputs["Strength"].default_value) if background else None,
        "cycles": cycles_settings,
        "taa_render_samples": int(scene_data.eevee.taa_render_samples),
        "view_transform": scene_data.view_settings.view_transform,
        "look": scene_data.view_settings.look,
        "exposure": float(scene_data.view_settings.exposure),
    }


def camera_view_dimensions(camera: bpy.types.Object) -> tuple[float, float]:
    aspect = bpy.context.scene.render.resolution_x / max(float(bpy.context.scene.render.resolution_y), 1.0)
    ortho_scale = float(camera.data.ortho_scale)
    if aspect >= 1.0:
        return ortho_scale, ortho_scale / aspect
    return ortho_scale * aspect, ortho_scale


def camera_space_renderable_bounds(camera_name: str, include_backdrop: bool = False) -> dict:
    camera = bpy.data.objects[camera_name]
    bpy.context.view_layer.update()
    inv = camera.matrix_world.inverted()
    xs = []
    ys = []
    zs = []
    object_count = 0
    for obj in bpy.data.objects:
        if obj.hide_render or obj.type in {"CAMERA", "LIGHT"}:
            continue
        if obj.get("canonical_beauty_backdrop") and not include_backdrop:
            continue
        object_count += 1
        for point in renderable_object_corners(obj):
            local = inv @ point
            xs.append(local.x)
            ys.append(local.y)
            zs.append(local.z)
    if not xs or not ys or not zs:
        return {"available": False, "object_count": object_count}
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)
    return {
        "available": True,
        "object_count": object_count,
        "min_x": min_x,
        "max_x": max_x,
        "min_y": min_y,
        "max_y": max_y,
        "min_z": min_z,
        "max_z": max_z,
        "width": max(max_x - min_x, 1e-6),
        "height": max(max_y - min_y, 1e-6),
        "depth": max(max_z - min_z, 1e-6),
    }


def canonical_beauty_collection(collections: dict[str, bpy.types.Collection]) -> bpy.types.Collection:
    collection = collections.get("Beauty") or bpy.data.collections.get("Beauty")
    if collection is None:
        collection = bpy.data.collections.new("Beauty")
    if collection.name not in {child.name for child in bpy.context.scene.collection.children}:
        bpy.context.scene.collection.children.link(collection)
    collections["Beauty"] = collection
    return collection


def remove_existing_canonical_beauty_objects() -> None:
    for obj in list(bpy.data.objects):
        if obj.name == "large_softbox" or obj.get("canonical_beauty_object"):
            data = obj.data
            bpy.data.objects.remove(obj, do_unlink=True)
            if isinstance(data, bpy.types.Light) and data.users == 0:
                bpy.data.lights.remove(data)
            elif isinstance(data, bpy.types.Mesh) and data.users == 0:
                bpy.data.meshes.remove(data)


def backdrop_material(name: str, color: tuple[float, float, float, float]) -> bpy.types.Material:
    material = bpy.data.materials.new(name)
    material.diffuse_color = color
    material.use_nodes = True
    nodes = material.node_tree.nodes
    bsdf = nodes.get("Principled BSDF")
    if bsdf:
        bsdf.inputs["Base Color"].default_value = color
        bsdf.inputs["Roughness"].default_value = 0.86
        bsdf.inputs["Metallic"].default_value = 0.0
        set_bsdf_input(bsdf, ("Specular IOR Level", "Specular"), 0.18)
        if "Emission Color" in bsdf.inputs:
            bsdf.inputs["Emission Color"].default_value = color
        elif "Emission" in bsdf.inputs:
            bsdf.inputs["Emission"].default_value = color
        if "Emission Strength" in bsdf.inputs:
            bsdf.inputs["Emission Strength"].default_value = 0.10
    return material


def blend_rgba(
    left: tuple[float, float, float, float],
    right: tuple[float, float, float, float],
    factor: float,
) -> tuple[float, float, float, float]:
    factor = max(0.0, min(1.0, factor))
    return tuple(float(a + (b - a) * factor) for a, b in zip(left, right))


def canonical_backdrop_color(x_norm: float, y_norm: float) -> tuple[float, float, float, float]:
    bottom = blend_rgba(Canonical_BACKDROP_BOTTOM_LEFT, Canonical_BACKDROP_BOTTOM_RIGHT, x_norm)
    top = blend_rgba(Canonical_BACKDROP_TOP_LEFT, Canonical_BACKDROP_TOP_RIGHT, x_norm)
    return blend_rgba(bottom, top, y_norm)


def create_canonical_backdrop(camera_name: str, collections: dict[str, bpy.types.Collection]) -> dict:
    camera = bpy.data.objects[camera_name]
    bounds = camera_space_renderable_bounds(camera_name)
    if not bounds.get("available"):
        return {"created": False, "reason": "no_renderable_objects"}

    collection = canonical_beauty_collection(collections)
    view_width, view_height = camera_view_dimensions(camera)
    backdrop_width = view_width * 1.34
    backdrop_height = view_height * 1.44
    backdrop_z = float(bounds["min_z"]) - 55.0
    camera.data.clip_end = max(float(camera.data.clip_end), abs(backdrop_z) + 120.0)

    columns = 24
    rows = 14
    vertices = []
    for row in range(rows + 1):
        y = (-0.5 + row / rows) * backdrop_height
        for column in range(columns + 1):
            x = (-0.5 + column / columns) * backdrop_width
            vertices.append(camera.matrix_world @ Vector((x, y, backdrop_z)))

    faces = []
    for row in range(rows):
        for column in range(columns):
            idx = row * (columns + 1) + column
            faces.append((idx, idx + 1, idx + columns + 2, idx + columns + 1))

    mesh = bpy.data.meshes.new("canonical_camera_gradient_backdrop_mesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    obj = bpy.data.objects.new("canonical_camera_gradient_backdrop", mesh)
    obj["canonical_beauty_object"] = True
    obj["canonical_beauty_backdrop"] = True
    obj["camera_aligned_to"] = camera_name
    obj.hide_select = True
    collection.objects.link(obj)

    face_index = 0
    for row in range(rows):
        for column in range(columns):
            x_norm = (column + 0.5) / columns
            y_norm = (row + 0.5) / rows
            mat = backdrop_material(f"canonical_backdrop_gradient_{row:02d}_{column:02d}", canonical_backdrop_color(x_norm, y_norm))
            mesh.materials.append(mat)
            mesh.polygons[face_index].material_index = face_index
            face_index += 1

    return {
        "created": True,
        "object": obj.name,
        "camera": camera_name,
        "grid": [columns, rows],
        "camera_local_z_mm": backdrop_z,
        "size_mm": [backdrop_width, backdrop_height],
        "palette": {
            "top_left": list(Canonical_BACKDROP_TOP_LEFT),
            "top_right": list(Canonical_BACKDROP_TOP_RIGHT),
            "bottom_left": list(Canonical_BACKDROP_BOTTOM_LEFT),
            "bottom_right": list(Canonical_BACKDROP_BOTTOM_RIGHT),
        },
    }


def point_object_at(obj: bpy.types.Object, target: Vector) -> None:
    direction = target - obj.location
    if direction.length > 1e-9:
        obj.rotation_euler = direction.to_track_quat("-Z", "Y").to_euler()


def create_area_light(
    name: str,
    collection: bpy.types.Collection,
    camera: bpy.types.Object,
    target: Vector,
    local_position: tuple[float, float, float],
    energy: float,
    size: float,
    color: tuple[float, float, float],
) -> dict:
    light_data = bpy.data.lights.new(name, "AREA")
    light_data.energy = energy
    light_data.size = size
    light_data.color = color
    if hasattr(light_data, "use_shadow"):
        light_data.use_shadow = True
    if hasattr(light_data, "use_contact_shadow"):
        light_data.use_contact_shadow = True
    light = bpy.data.objects.new(name, light_data)
    light.location = camera.matrix_world @ Vector(local_position)
    point_object_at(light, target)
    light["canonical_beauty_object"] = True
    collection.objects.link(light)
    return {
        "name": name,
        "type": "AREA",
        "energy": energy,
        "size_mm": size,
        "color": list(color),
        "camera_local_position_mm": list(local_position),
        "world_location_mm": [light.location.x, light.location.y, light.location.z],
    }


def add_canonical_beauty_environment(camera_name: str, collections: dict[str, bpy.types.Collection]) -> dict:
    remove_existing_canonical_beauty_objects()
    collection = canonical_beauty_collection(collections)
    camera = bpy.data.objects[camera_name]
    bounds = camera_space_renderable_bounds(camera_name)
    if not bounds.get("available"):
        return {"created": False, "reason": "no_renderable_objects"}

    view_width, view_height = camera_view_dimensions(camera)
    scene_center_local = Vector(
        (
            (float(bounds["min_x"]) + float(bounds["max_x"])) * 0.5,
            (float(bounds["min_y"]) + float(bounds["max_y"])) * 0.5,
            (float(bounds["min_z"]) + float(bounds["max_z"])) * 0.5,
        )
    )
    scene_center_world = camera.matrix_world @ scene_center_local
    front_z = float(bounds["max_z"]) + 42.0
    rim_z = float(bounds["min_z"]) - 24.0
    lights = [
        create_area_light(
            "canonical_key_softbox",
            collection,
            camera,
            scene_center_world,
            (-0.34 * view_width, 0.32 * view_height, front_z),
            10500.0,
            74.0,
            (1.0, 0.94, 0.84),
        ),
        create_area_light(
            "canonical_fill_softbox",
            collection,
            camera,
            scene_center_world,
            (0.42 * view_width, -0.14 * view_height, front_z + 14.0),
            2600.0,
            210.0,
            (0.84, 0.91, 1.0),
        ),
        create_area_light(
            "canonical_rim_softbox",
            collection,
            camera,
            scene_center_world,
            (0.26 * view_width, 0.48 * view_height, rim_z),
            6200.0,
            56.0,
            (0.78, 0.86, 1.0),
        ),
        create_area_light(
            "canonical_warm_wash",
            collection,
            camera,
            scene_center_world,
            (-0.08 * view_width, -0.46 * view_height, front_z + 26.0),
            2600.0,
            210.0,
            (1.0, 0.88, 0.74),
        ),
    ]
    backdrop = create_canonical_backdrop(camera_name, collections)
    bpy.context.view_layer.update()
    return {
        "created": True,
        "collection": collection.name,
        "camera": camera_name,
        "scene_center_world_mm": [scene_center_world.x, scene_center_world.y, scene_center_world.z],
        "lights": lights,
        "backdrop": backdrop,
    }


def polish_canonical_surface_rendering() -> dict:
    smoothed_meshes = []
    weighted_normal_meshes = []
    for obj in bpy.data.objects:
        if obj.type != "MESH" or obj.hide_render or obj.get("canonical_beauty_backdrop"):
            continue
        if obj.data and obj.data.polygons:
            for polygon in obj.data.polygons:
                polygon.use_smooth = True
            smoothed_meshes.append(obj.name)
        if obj.modifiers.get("canonical_weighted_normals") is None:
            modifier = obj.modifiers.new("canonical_weighted_normals", "WEIGHTED_NORMAL")
            modifier.keep_sharp = True
            modifier.weight = 50
            weighted_normal_meshes.append(obj.name)
    bpy.context.view_layer.update()
    return {
        "policy": "visible molecular meshes use smooth shading plus weighted normals for sciart-style surface form; geometry scale is unchanged",
        "smoothed_mesh_count": len(smoothed_meshes),
        "weighted_normal_modifier_count": len(weighted_normal_meshes),
        "sample_smoothed_meshes": smoothed_meshes[:12],
    }


def apply_canonical_color_palette(materials: dict[str, bpy.types.Material]) -> dict:
    applied = {}
    for name, color in Canonical_MATERIAL_COLORS.items():
        material = materials.get(name)
        if material is None:
            material = bpy.data.materials.new(name)
            materials[name] = material
        set_material_color(material, color)
        applied[name] = list(color)
    return {
        "policy": "DNA promoter/exons/introns remain distinct; mRNA UTRs/coding remain distinct; context proteins are uniform slate; actin product is highlighted coral",
        "dna_order_colors": {
            "promoter": "dna_promoter_blue",
            "exons": "dna_exon_orange",
            "introns": "dna_intron_olive",
        },
        "mrna_order_colors": {
            "5_utr": "olive",
            "coding_sequence": "orange",
            "3_utr": "yellow_olive",
        },
        "protein_colors": {
            "context_proteins": "protein_uniform_slate",
            "actin_product": "protein_product_coral",
        },
        "applied_material_rgba": applied,
    }


def renderable_object_corners(obj: bpy.types.Object) -> list[Vector]:
    if obj.type in {"MESH", "CURVE", "FONT"} and obj.bound_box:
        return [obj.matrix_world @ Vector(corner) for corner in obj.bound_box]
    return [obj.matrix_world.translation]


def orient_labels_to_camera(camera_name: str, *, reveal: bool = True) -> dict:
    camera = bpy.data.objects[camera_name]
    rows = []
    for obj in bpy.data.objects:
        if obj.type != "FONT":
            continue
        if reveal:
            obj.hide_render = False
        obj.rotation_euler = camera.rotation_euler
        obj["oriented_to_camera"] = camera_name
        rows.append({"object": obj.name, "text": obj.data.body})
    bpy.context.view_layer.update()
    return {
        "policy": "Canonical text labels are billboarded to the active camera",
        "camera": camera_name,
        "oriented_label_count": len(rows),
        "rows": rows,
    }


def fit_camera_to_renderables(camera_name: str, margin_fraction: float = 0.075) -> dict:
    camera = bpy.data.objects[camera_name]
    bounds = camera_space_renderable_bounds(camera_name)
    if not bounds.get("available"):
        return {"camera": camera_name, "fit_applied": False, "reason": "no_renderable_objects"}
    min_x = float(bounds["min_x"])
    max_x = float(bounds["max_x"])
    min_y = float(bounds["min_y"])
    max_y = float(bounds["max_y"])
    width = float(bounds["width"])
    height = float(bounds["height"])
    aspect = bpy.context.scene.render.resolution_x / max(float(bpy.context.scene.render.resolution_y), 1.0)
    usable = max(1.0 - 2.0 * margin_fraction, 0.1)
    old_location = [camera.location.x, camera.location.y, camera.location.z]
    old_ortho = float(camera.data.ortho_scale)
    local_center = Vector(((min_x + max_x) * 0.5, (min_y + max_y) * 0.5, 0.0))
    camera.location += camera.matrix_world.to_3x3() @ local_center
    if aspect >= 1.0:
        required_ortho = max(width, height * aspect)
    else:
        required_ortho = max(width / aspect, height)
    camera.data.ortho_scale = required_ortho / usable
    if camera_name == "Camera_canonical_full_overview":
        camera.data.ortho_scale *= 1.08
        gutter_shift = camera.matrix_world.to_3x3() @ Vector((float(camera.data.ortho_scale) * 0.065, 0.0, 0.0))
        camera.location += gutter_shift
    bpy.context.view_layer.update()
    return {
        "camera": camera_name,
        "fit_applied": True,
        "renderable_object_count": int(bounds["object_count"]),
        "margin_fraction": margin_fraction,
        "local_bounds_before_fit": {
            "min_x": min_x,
            "max_x": max_x,
            "min_y": min_y,
            "max_y": max_y,
            "min_z": float(bounds["min_z"]),
            "max_z": float(bounds["max_z"]),
            "width": width,
            "height": height,
            "depth": float(bounds["depth"]),
        },
        "old_location_mm": old_location,
        "new_location_mm": [camera.location.x, camera.location.y, camera.location.z],
        "old_ortho_scale_mm": old_ortho,
        "new_ortho_scale_mm": float(camera.data.ortho_scale),
    }


def _camera_overlay_point(camera: bpy.types.Object, x_norm: float, y_norm: float, z_local: float = -1.0) -> Vector:
    scene_render = bpy.context.scene.render
    aspect = scene_render.resolution_x / max(float(scene_render.resolution_y), 1.0)
    if aspect >= 1.0:
        local = Vector(
            (
                (x_norm - 0.5) * float(camera.data.ortho_scale),
                (y_norm - 0.5) * float(camera.data.ortho_scale) / aspect,
                z_local,
            )
        )
    else:
        local = Vector(
            (
                (x_norm - 0.5) * float(camera.data.ortho_scale) * aspect,
                (y_norm - 0.5) * float(camera.data.ortho_scale),
                z_local,
            )
        )
    return camera.matrix_world @ local


def _remove_annotation_guides() -> None:
    for obj in list(bpy.data.objects):
        if obj.get("canonical_annotation_leader") or obj.get("canonical_annotation_bracket") or obj.get("canonical_label_backing"):
            bpy.data.objects.remove(obj, do_unlink=True)


def _projected_object_bounds(camera: bpy.types.Object, label: bpy.types.Object) -> tuple[float, float, float, float]:
    inv = camera.matrix_world.inverted()
    view_width, view_height = camera_view_dimensions(camera)
    points = []
    pdb_id = label.get("pdb_id")
    for obj in bpy.data.objects:
        if obj.type != "MESH" or not pdb_id or obj.get("pdb_id") != pdb_id:
            continue
        points.extend(obj.matrix_world @ Vector(corner) for corner in obj.bound_box)
    if not points:
        anchor = Vector(label.get("molecule_anchor_mm", label.location))
        points = [anchor]
    normalized = []
    for point in points:
        local = inv @ point
        normalized.append((local.x / view_width + 0.5, local.y / view_height + 0.5))
    return (
        min(point[0] for point in normalized), max(point[0] for point in normalized),
        min(point[1] for point in normalized), max(point[1] for point in normalized),
    )


def _projected_objects_bounds(camera: bpy.types.Object, objects: list[bpy.types.Object]) -> tuple[float, float, float, float]:
    inv = camera.matrix_world.inverted()
    view_width, view_height = camera_view_dimensions(camera)
    normalized = []
    for obj in objects:
        for point in renderable_object_corners(obj):
            local = inv @ point
            normalized.append((local.x / view_width + 0.5, local.y / view_height + 0.5))
    if not normalized:
        raise ValueError("Cannot project bounds for an empty object collection")
    return (
        min(point[0] for point in normalized), max(point[0] for point in normalized),
        min(point[1] for point in normalized), max(point[1] for point in normalized),
    )


def _boxes_overlap(a: tuple[float, float, float, float], b: tuple[float, float, float, float], pad: float = 0.006) -> bool:
    return not (a[1] + pad < b[0] or b[1] + pad < a[0] or a[3] + pad < b[2] or b[3] + pad < a[2])


def _label_box(center: tuple[float, float], width: float, height: float) -> tuple[float, float, float, float]:
    return (center[0] - width * 0.5, center[0] + width * 0.5, center[1] - height * 0.5, center[1] + height * 0.5)


def _ensure_label_backing_material(materials: dict) -> bpy.types.Material:
    if "label_backing" in materials:
        return materials["label_backing"]
    mat = bpy.data.materials.new("label_backing")
    color = (0.975, 0.982, 0.985, 0.72)
    mat.diffuse_color = color
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf:
        bsdf.inputs["Base Color"].default_value = color
        bsdf.inputs["Roughness"].default_value = 0.92
        if "Alpha" in bsdf.inputs:
            bsdf.inputs["Alpha"].default_value = color[3]
    if hasattr(mat, "surface_render_method"):
        mat.surface_render_method = "DITHERED"
    materials["label_backing"] = mat
    return mat


def _create_label_backing(name: str, camera: bpy.types.Object, box, collections: dict, materials: dict) -> str:
    left, right, bottom, top = box
    margin_x, margin_y = 0.0025, 0.002
    corners = [
        _camera_overlay_point(camera, left - margin_x, bottom - margin_y, -1.08),
        _camera_overlay_point(camera, right + margin_x, bottom - margin_y, -1.08),
        _camera_overlay_point(camera, right + margin_x, top + margin_y, -1.08),
        _camera_overlay_point(camera, left - margin_x, top + margin_y, -1.08),
    ]
    mesh = bpy.data.meshes.new(f"{name}_mesh")
    mesh.from_pydata([tuple(point) for point in corners], [], [(0, 1, 2, 3)])
    obj = bpy.data.objects.new(name, mesh)
    collections["Labels"].objects.link(obj)
    obj.data.materials.append(_ensure_label_backing_material(materials))
    obj["canonical_label_backing"] = True
    return obj.name


def _create_group_bracket(
    name: str,
    camera: bpy.types.Object,
    span: tuple[float, float],
    collections: dict,
    materials: dict,
) -> str:
    bracket_x = 0.805
    tick = 0.018
    bottom, top = span
    points = [
        _camera_overlay_point(camera, bracket_x - tick, top),
        _camera_overlay_point(camera, bracket_x, top),
        _camera_overlay_point(camera, bracket_x, bottom),
        _camera_overlay_point(camera, bracket_x - tick, bottom),
    ]
    bracket = base.create_curve(
        name,
        [(point.x, point.y, point.z) for point in points],
        0.060,
        materials["label_grey"],
        collections["Labels"],
        resolution=1,
    )
    bracket["canonical_annotation_bracket"] = True
    return bracket.name


def place_overview_labels(camera_name: str, collections: dict, materials: dict) -> dict:
    camera = bpy.data.objects[camera_name]
    inv = camera.matrix_world.inverted()
    rows = []
    _remove_annotation_guides()
    if bpy.data.objects.get("label_ACTB_primary_canonical") is None:
        primary = base.create_text(
            "label_ACTB_primary_canonical", "ACTB protein 375 aa", (0.0, 0.0, 0.0), 1.55,
            materials["black"], collections["Labels"], align="LEFT"
        )
        primary.rotation_euler = camera.rotation_euler

    priority_labels = {
        "label_mCherry/RFP tag": 0,
        "label_Ribosome small subunit": 0,
        "label_Ribosome large subunit": 0,
    }
    labels = sorted(
        [obj for obj in bpy.data.objects if obj.type == "FONT" and obj.get("pdb_id")],
        key=lambda obj: (priority_labels.get(obj.name, 1), -_projected_object_bounds(camera, obj)[3]),
    )
    molecule_boxes = {obj.name: _projected_object_bounds(camera, obj) for obj in labels}
    occupied = []
    for obj in labels:
        molecule_box = molecule_boxes[obj.name]
        anchor_norm = ((molecule_box[0] + molecule_box[1]) * 0.5, (molecule_box[2] + molecule_box[3]) * 0.5)
        lines = obj.data.body.splitlines() or [obj.data.body]
        box_width = max(0.026, min(0.078, 0.006 + max(len(line) for line in lines) * 0.00215))
        box_height = max(0.017, len(lines) * 0.017)
        gap = 0.008 if obj.name == "label_HuR-like RBP" else 0.003
        half_width, half_height = box_width * 0.5, box_height * 0.5
        candidates = [
            (molecule_box[1] + gap + half_width, anchor_norm[1]),
            (anchor_norm[0], molecule_box[3] + gap + half_height),
            (molecule_box[0] - gap - half_width, anchor_norm[1]),
            (anchor_norm[0], molecule_box[2] - gap - half_height),
            (molecule_box[1] + gap + half_width, molecule_box[3] + gap + half_height),
            (molecule_box[0] - gap - half_width, molecule_box[3] + gap + half_height),
            (molecule_box[1] + gap + half_width, molecule_box[2] - gap - half_height),
            (molecule_box[0] - gap - half_width, molecule_box[2] - gap - half_height),
        ]
        for extra_gap in (gap + 0.012, gap + 0.024, gap + 0.036):
            candidates.extend(
                [
                    (molecule_box[1] + extra_gap + half_width, anchor_norm[1]),
                    (anchor_norm[0], molecule_box[3] + extra_gap + half_height),
                    (molecule_box[0] - extra_gap - half_width, anchor_norm[1]),
                    (anchor_norm[0], molecule_box[2] - extra_gap - half_height),
                    (molecule_box[1] + extra_gap + half_width, molecule_box[3] + extra_gap + half_height),
                    (molecule_box[0] - extra_gap - half_width, molecule_box[3] + extra_gap + half_height),
                    (molecule_box[1] + extra_gap + half_width, molecule_box[2] - extra_gap - half_height),
                    (molecule_box[0] - extra_gap - half_width, molecule_box[2] - extra_gap - half_height),
                ]
            )
        valid = []
        for candidate in candidates:
            box = _label_box(candidate, box_width, box_height)
            inside = box[0] >= 0.08 and box[1] <= 0.775 and box[2] >= 0.08 and box[3] <= 0.93
            ribosome_cluster = {"label_Ribosome small subunit", "label_Ribosome large subunit"}
            overlaps_other_molecule = any(
                other_name != obj.name
                and not ({obj.name, other_name} <= ribosome_cluster)
                and _boxes_overlap(box, other_box, pad=0.002)
                for other_name, other_box in molecule_boxes.items()
            )
            if inside and not overlaps_other_molecule and not any(_boxes_overlap(box, other) for other in occupied):
                valid.append((candidate, box))
        if not valid:
            for candidate in candidates:
                box = _label_box(candidate, box_width, box_height)
                inside = box[0] >= 0.08 and box[1] <= 0.775 and box[2] >= 0.08 and box[3] <= 0.93
                if inside and not any(_boxes_overlap(box, other) for other in occupied):
                    valid.append((candidate, box))
        if valid:
            (x_norm, y_norm), placed_box = min(valid, key=lambda item: abs(item[0][0] - anchor_norm[0]) + abs(item[0][1] - anchor_norm[1]))
        else:
            x_norm = max(0.08 + box_width * 0.5, min(0.775 - box_width * 0.5, candidates[0][0]))
            y_norm = max(0.08 + box_height * 0.5, min(0.93 - box_height * 0.5, candidates[0][1]))
            placed_box = _label_box((x_norm, y_norm), box_width, box_height)
            placed_box = min(
                (_label_box(candidate, box_width, box_height) for candidate in candidates),
                key=lambda box: sum(1 for other in occupied if _boxes_overlap(box, other)),
            )
            x_norm = (placed_box[0] + placed_box[1]) * 0.5
            y_norm = (placed_box[2] + placed_box[3]) * 0.5
        occupied.append(placed_box)
        molecule_anchor = Vector(obj.get("molecule_anchor_mm", obj.location))
        molecule_depth = (inv @ molecule_anchor).z
        scene_depth = molecule_depth + 1.2
        obj.location = _camera_overlay_point(camera, x_norm, y_norm, scene_depth)
        obj.data.size = min(float(obj.data.size), 1.00)
        obj.data.align_x = "CENTER"
        obj.data.align_y = "CENTER"
        displacement = abs(anchor_norm[0] - x_norm) + abs(anchor_norm[1] - y_norm)
        obj["overview_label_position"] = [x_norm, y_norm]
        rows.append(
            {
                "object": obj.name,
                "anchor_position": [anchor_norm[0], anchor_norm[1]],
                "view_position": [x_norm, y_norm],
                "estimated_view_box": [box_width, box_height],
                "molecule_projected_bounds": list(molecule_box),
                "placement_space": "world_scene_at_molecule_depth",
                "camera_local_depth_mm": scene_depth,
                "leader": None,
                "displacement": displacement,
            }
        )

    compact_obj = bpy.data.objects.get(COMPACT_CALLOUT["object"])
    compact_row = None
    if compact_obj is not None:
        compact_meshes = [
            obj
            for obj in bpy.data.objects
            if obj.type == "MESH" and obj.name.startswith("actin mRNA compact")
        ]
        compact_bounds = _projected_objects_bounds(camera, compact_meshes)
        compact_anchor = (compact_bounds[1], compact_bounds[3])
        compact_position = (
            compact_anchor[0] + COMPACT_CALLOUT["offset"][0],
            compact_anchor[1] + COMPACT_CALLOUT["offset"][1],
        )
        label_box = (
            compact_position[0], compact_position[0] + 0.060,
            compact_position[1] - 0.010, compact_position[1] + 0.010,
        )
        compact_local = inv @ Vector(compact_obj.get("molecule_anchor_mm", compact_obj.location))
        compact_obj.data.body = COMPACT_CALLOUT["text"]
        compact_obj.data.align_x = "LEFT"
        compact_obj.data.align_y = "CENTER"
        compact_obj.data.size = 1.34
        compact_obj.location = _camera_overlay_point(camera, *compact_position, compact_local.z + 1.2)
        compact_obj["overview_label_position"] = list(compact_position)
        compact_row = {
            "object": COMPACT_CALLOUT["object"],
            "text": COMPACT_CALLOUT["text"],
            "anchor_position": list(compact_anchor),
            "view_position": list(compact_position),
            "offset": list(COMPACT_CALLOUT["offset"]),
            "molecule_projected_bounds": list(compact_bounds),
            "estimated_view_box": list(label_box),
            "overlaps_molecule": _boxes_overlap(label_box, compact_bounds, pad=0.004),
            "leader": None,
            "placement_space": "world_scene_at_molecule_depth",
        }

    primary_rows = []
    for object_name, spec in PRIMARY_CALLOUTS.items():
        obj = bpy.data.objects.get(object_name)
        if obj is None:
            continue
        x_norm, y_norm = spec["view_position"]
        obj.data.body = spec["text"]
        obj.data.align_x = "LEFT"
        obj.data.align_y = "CENTER"
        obj.data.size = 1.55
        obj.location = _camera_overlay_point(camera, x_norm, y_norm)
        bracket_name = _create_group_bracket(
            f"overview_group_bracket_{object_name}", camera, spec["span"], collections, materials
        )
        obj["overview_label_position"] = [x_norm, y_norm]
        primary_rows.append(
            {
                "object": object_name,
                "text": spec["text"],
                "view_position": [x_norm, y_norm],
                "span": list(spec["span"]),
                "bracket": bracket_name,
                "leader": None,
            }
        )
    bpy.context.view_layer.update()
    collision_pairs = []
    molecule_overlap_pairs = []
    for left_index, left in enumerate(rows):
        left_box = _label_box(tuple(left["view_position"]), *left["estimated_view_box"])
        for right in rows[left_index + 1:]:
            right_box = _label_box(tuple(right["view_position"]), *right["estimated_view_box"])
            if _boxes_overlap(left_box, right_box, pad=0.002):
                collision_pairs.append([left["object"], right["object"]])
        for molecule_name, molecule_box in molecule_boxes.items():
            ribosome_cluster = {"label_Ribosome small subunit", "label_Ribosome large subunit"}
            if (
                molecule_name != left["object"]
                and not ({left["object"], molecule_name} <= ribosome_cluster)
                and _boxes_overlap(left_box, molecule_box, pad=0.001)
            ):
                molecule_overlap_pairs.append([left["object"], molecule_name])
    return {
        "policy": "three right-side camera-space brackets group protein, mRNA, and DNA; PDB labels remain in world space at molecule depth and are camera-packed immediately outside projected molecule bounds without leaders or backings",
        "camera": camera_name,
        "placed_label_count": len(rows) + len(primary_rows) + (1 if compact_row else 0),
        "primary_callouts": primary_rows,
        "compact_callout": compact_row,
        "rows": rows,
        "collision_pairs": collision_pairs,
        "molecule_overlap_pairs": molecule_overlap_pairs,
        "maximum_displacement": max((row["displacement"] for row in rows), default=0.0),
        "maximum_displacement_object": max(rows, key=lambda row: row["displacement"])["object"] if rows else None,
    }


def _asset_report_location(asset_reports: list[dict], name: str) -> Vector:
    for item in asset_reports:
        if item.get("name") == name:
            return Vector(item["location_mm"])
    raise KeyError(f"Missing asset report for camera target: {name}")


def add_canonical_cameras(
    manifest: dict,
    dna_path,
    mrna_path,
    compact_path,
    asset_reports: list[dict],
    collections: dict,
    materials: dict,
) -> dict[str, str]:
    base.add_lighting_and_camera(collections, materials)
    dna_center = contact_helpers.path_center(dna_path)
    p53 = _asset_report_location(asset_reports, "p53 tetramer bound to DNA")
    polymerase = _asset_report_location(asset_reports, "RNA polymerase II elongation complex")
    nucleosome = _asset_report_location(asset_reports, "Nucleosome")
    ribosome = (
        _asset_report_location(asset_reports, "Ribosome large subunit")
        + _asset_report_location(asset_reports, "Ribosome small subunit")
        + _asset_report_location(asset_reports, "Standalone tRNA")
    ) * (1.0 / 3.0)
    actin = _asset_report_location(asset_reports, "Actin protein")
    cas9 = _asset_report_location(asset_reports, "Cas9")

    camera_names = {
        "full_overview": create_camera(
            "Camera_canonical_full_overview",
            (dna_center.x - 48.0, dna_center.y - 142.0, 128.0),
            Vector((dna_center.x + 12.0, dna_center.y - 6.0, 36.0)),
            190.0,
        ),
        "p53_dna": create_camera(
            "Camera_canonical_p53_dna",
            (p53.x - 8.0, p53.y - 12.0, p53.z + 9.0),
            p53,
            9.0,
        ),
        "polymerase_rna_start": create_camera(
            "Camera_canonical_polymerase_rna_start",
            (polymerase.x - 10.0, polymerase.y - 14.0, polymerase.z + 10.0),
            polymerase,
            11.0,
        ),
        "nucleosome_loop": create_camera(
            "Camera_canonical_nucleosome_loop",
            (nucleosome.x - 9.0, nucleosome.y - 13.0, nucleosome.z + 10.0),
            nucleosome,
            10.0,
        ),
        "ribosome_trna": create_camera(
            "Camera_canonical_ribosome_trna",
            (ribosome.x - 13.0, ribosome.y - 18.0, ribosome.z + 15.0),
            ribosome,
            17.0,
        ),
        "actin_product": create_camera(
            "Camera_canonical_actin_product",
            (actin.x - 7.0, actin.y - 10.0, actin.z + 8.0),
            actin,
            7.0,
        ),
        "cas9_dna": create_camera(
            "Camera_canonical_cas9_dna",
            (cas9.x - 9.0, cas9.y - 13.0, cas9.z + 10.0),
            cas9,
            DETAIL_SHARED_ORTHO_SCALE_MM,
        ),
    }
    bpy.context.scene.camera = bpy.data.objects[camera_names["full_overview"]]
    return camera_names


def _detail_path_curve(name: str, path, start_fraction: float, end_fraction: float, material, collection) -> str:
    sample_count = 180
    start = path.length * start_fraction
    end = path.length * end_fraction
    points = [path.point_at_length(start + (end - start) * index / sample_count) for index in range(sample_count + 1)]
    obj = base.create_curve(
        name,
        [(point.x, point.y, point.z) for point in points],
        0.13,
        material,
        collection,
        resolution=2,
    )
    obj.hide_render = True
    obj["detail_context"] = True
    return obj.name


def add_detail_context_curves(mrna_path, collections: dict, materials: dict) -> dict:
    collection = collections.get("Detail context")
    if collection is None:
        collection = bpy.data.collections.new("Detail context")
        bpy.context.scene.collection.children.link(collection)
        collections["Detail context"] = collection
    return {
        "polymerase_rna_start": _detail_path_curve(
            "detail_polymerase_nascent_rna",
            mrna_path,
            0.0,
            0.055,
            materials["rna_gold"],
            collection,
        ),
        "ribosome_trna": _detail_path_curve(
            "detail_ribosome_mrna",
            mrna_path,
            0.855,
            0.935,
            materials["rna_gold"],
            collection,
        ),
    }


FOCUS_OBJECT_PATTERNS = {
    "p53_dna": ("p53 tetramer bound to DNA (3TS8)",),
    "nucleosome_loop": ("Nucleosome (1AOI)",),
    "polymerase_rna_start": ("RNA polymerase II elongation complex (2E2I)", "detail_polymerase_nascent_rna"),
    "ribosome_trna": (
        "Ribosome small subunit (1J5E)",
        "Ribosome large subunit (1JJ2)",
        "Standalone tRNA (4TNA)",
        "detail_ribosome_mrna",
    ),
    "actin_product": ("Actin protein (1J6Z)",),
    "cas9_dna": ("Cas9 (4UN3)",),
}


def snapshot_render_visibility() -> dict[str, bool]:
    return {obj.name: bool(obj.hide_render) for obj in bpy.data.objects}


def restore_render_visibility(snapshot: dict[str, bool]) -> None:
    for name, hidden in snapshot.items():
        obj = bpy.data.objects.get(name)
        if obj is not None:
            obj.hide_render = hidden


def apply_focus_visibility(key: str) -> dict:
    patterns = FOCUS_OBJECT_PATTERNS[key]
    visible = []
    for obj in bpy.data.objects:
        if obj.type in {"CAMERA", "LIGHT"} or obj.get("canonical_beauty_object"):
            continue
        should_show = any(pattern in obj.name for pattern in patterns)
        obj.hide_render = not should_show
        if should_show:
            visible.append(obj.name)
    return {"key": key, "patterns": list(patterns), "visible_objects": visible}


def create_detail_title(key: str, camera_name: str, collections: dict, materials: dict) -> str:
    for obj in list(bpy.data.objects):
        if obj.get("canonical_detail_title"):
            bpy.data.objects.remove(obj, do_unlink=True)
    camera = bpy.data.objects[camera_name]
    obj = base.create_text(
        f"detail_title_{key}",
        DETAIL_TITLES[key],
        (0.0, 0.0, 0.0),
        max(0.32, float(camera.data.ortho_scale) * 0.036),
        materials["black"],
        collections["Labels"],
        align="LEFT",
    )
    obj.location = _camera_overlay_point(camera, 0.07, 0.91)
    obj.rotation_euler = camera.rotation_euler
    obj.hide_render = False
    obj["canonical_detail_title"] = True
    obj["detail_view"] = key
    if key == "ribosome_trna":
        source = bpy.data.objects.get("label_Standalone tRNA")
        if source is not None:
            anchor_local = camera.matrix_world.inverted() @ Vector(source.get("molecule_anchor_mm", source.location))
            view_width, view_height = camera_view_dimensions(camera)
            anchor_norm = (anchor_local.x / view_width + 0.5, anchor_local.y / view_height + 0.5)
            label_position = (
                max(0.10, min(0.74, anchor_norm[0] + 0.045)),
                max(0.12, min(0.84, anchor_norm[1] + 0.025)),
            )
            trna = base.create_text(
                "detail_label_trna",
                "yeast tRNA-Phe (4TNA)",
                (0.0, 0.0, 0.0),
                max(0.24, float(camera.data.ortho_scale) * 0.023),
                materials["black"],
                collections["Labels"],
                align="LEFT",
            )
            trna.location = _camera_overlay_point(camera, *label_position)
            trna.rotation_euler = camera.rotation_euler
            trna.hide_render = False
            trna["canonical_detail_title"] = True
            trna["detail_view"] = key
    return obj.name


def path_bounds(path) -> dict:
    points = path.points
    min_v = Vector((min(point.x for point in points), min(point.y for point in points), min(point.z for point in points)))
    max_v = Vector((max(point.x for point in points), max(point.y for point in points), max(point.z for point in points)))
    return {
        "min_mm": [min_v.x, min_v.y, min_v.z],
        "max_mm": [max_v.x, max_v.y, max_v.z],
        "size_mm": [max_v.x - min_v.x, max_v.y - min_v.y, max_v.z - min_v.z],
    }


def source_curve_visibility_validation(report: dict) -> dict:
    rows = []
    failures = []
    for path_name, object_name in (("dna", "Canonical_DNA_source_path"), ("mrna", "Canonical_mRNA_source_path")):
        obj = bpy.data.objects.get(object_name)
        row = {
            "path": path_name,
            "object": object_name,
            "exists": obj is not None,
            "hide_viewport": bool(obj.hide_viewport) if obj else None,
            "hide_render": bool(obj.hide_render) if obj else None,
        }
        rows.append(row)
        if obj is None:
            failures.append({"object": object_name, "reason": "missing_source_curve"})
        elif not obj.hide_viewport or not obj.hide_render:
            failures.append({"object": object_name, "reason": "source_curve_not_hidden", "row": row})
    return {"policy": "source paths kept for editing but hidden in viewport and render", "rows": rows, "failures": failures}


def reader_order_validation(dna_path) -> dict:
    bounds = path_bounds(dna_path)
    start = dna_path.points[0]
    end = dna_path.points[-1]
    width = max(bounds["size_mm"][0], 1e-6)
    height = max(bounds["size_mm"][1], 1e-6)
    min_x, min_y, _min_z = bounds["min_mm"]
    max_x, max_y, _max_z = bounds["max_mm"]
    start_top_left = start.x <= min_x + width * 0.20 and start.y >= max_y - height * 0.20
    end_lower_right = end.x > start.x and end.y < start.y
    failures = []
    if not start_top_left:
        failures.append({"reason": "promoter_start_not_top_left"})
    if not end_lower_right:
        failures.append({"reason": "path_end_not_lower_right_than_start"})
    return {
        "path_mode": "canonical_reader_order_serpentine_with_nucleosome_loop",
        "bounds": bounds,
        "start_mm": [start.x, start.y, start.z],
        "end_mm": [end.x, end.y, end.z],
        "start_is_top_left": start_top_left,
        "end_is_lower_right_than_start": end_lower_right,
        "failures": failures,
    }


def dna_cluster_validation(asset_reports: list[dict]) -> dict:
    by_name = {item["name"]: item for item in asset_reports}
    rows = []
    failures = []
    for name in sorted(canonical.DNA_STRICT_CONTACT_NAMES):
        item = by_name.get(name)
        attachment = (item or {}).get("attachment_empty") or {}
        fraction = attachment.get("fraction")
        row = {"name": name, "path": attachment.get("path"), "fraction": fraction, "max_allowed_fraction": 0.205}
        rows.append(row)
        if not item:
            failures.append({"name": name, "reason": "missing_asset_report"})
        elif attachment.get("path") != "dna":
            failures.append({"name": name, "reason": "not_attached_to_dna", "actual_path": attachment.get("path")})
        elif fraction is None or float(fraction) > 0.2055:
            failures.append({"name": name, "reason": "outside_early_gene_cluster", "fraction": fraction, "max": 0.205})
    return {"policy": "all DNA regulatory binders cluster in the first 20.5 percent of Canonical DNA", "rows": rows, "failures": failures}


def compact_position_validation(positioning: dict, compact_center: list[float]) -> dict:
    start = Vector(positioning["mrna_end_mm"])
    end = Vector(positioning["actin_location_mm"])
    center = Vector(compact_center)
    segment = end - start
    length = segment.length
    projection = 0.0 if length == 0 else (center - start).dot(segment) / (length * length)
    closest = start + segment * max(0.0, min(1.0, projection))
    distance_to_segment = (center - closest).length
    z_low = min(start.z, end.z)
    z_high = max(start.z, end.z)
    failures = []
    if projection < -0.05 or projection > 1.05:
        failures.append({"reason": "compact_mrna_not_between_mrna_end_and_actin", "projection": projection})
    if center.z < z_low - 0.5 or center.z > z_high + 0.5:
        failures.append({"reason": "compact_mrna_z_outside_mrna_end_to_actin_range", "center_z": center.z, "range": [z_low, z_high]})
    if distance_to_segment > 2.5:
        failures.append({"reason": "compact_mrna_center_off_interpolation_line", "distance_mm": distance_to_segment, "max_mm": 2.5})
    return {
        **positioning,
        "actual_center_mm": [center.x, center.y, center.z],
        "projection_from_mrna_end_to_actin": projection,
        "distance_to_mrna_end_actin_segment_mm": distance_to_segment,
        "failures": failures,
    }


def trna_ribosome_separation_report(asset_reports: list[dict]) -> dict:
    by_name = {item["name"]: item for item in asset_reports}
    trna = by_name["Standalone tRNA"]
    ribosomes = [by_name["Ribosome small subunit"], by_name["Ribosome large subunit"]]
    trna_location = Vector(trna["location_mm"])
    distances = {item["name"]: (trna_location - Vector(item["location_mm"])).length for item in ribosomes}
    fractions = {
        item["name"]: (item.get("attachment_empty") or {}).get("fraction")
        for item in [trna, *ribosomes]
    }
    minimum = min(distances.values())
    def component_bounds(item: dict) -> tuple[Vector, Vector]:
        components = item.get("components", [])
        return (
            Vector(tuple(min(component["min_mm"][axis] for component in components) for axis in range(3))),
            Vector(tuple(max(component["max_mm"][axis] for component in components) for axis in range(3))),
        )

    def aabb_gap(left: tuple[Vector, Vector], right: tuple[Vector, Vector]) -> float:
        delta = Vector(tuple(max(0.0, left[0][axis] - right[1][axis], right[0][axis] - left[1][axis]) for axis in range(3)))
        return delta.length

    trna_bounds = component_bounds(trna)
    surface_gaps = {item["name"]: aabb_gap(trna_bounds, component_bounds(item)) for item in ribosomes}
    minimum_surface_gap = min(surface_gaps.values())
    failures = []
    if minimum < 4.0:
        failures.append({"reason": "trna_not_visibly_separated_from_ribosome", "minimum_center_distance_mm": minimum, "minimum_required_mm": 4.0})
    if minimum_surface_gap < 0.8:
        failures.append({"reason": "trna_surface_silhouette_not_separated_from_ribosome", "minimum_surface_bbox_gap_mm": minimum_surface_gap, "minimum_required_mm": 0.8})
    return {
        "policy": "incoming tRNA is displayed beside the ribosome without a misleading visual connection to the mRNA path",
        "center_distances_mm": distances,
        "minimum_center_distance_mm": minimum,
        "surface_bbox_gaps_mm": surface_gaps,
        "minimum_surface_bbox_gap_mm": minimum_surface_gap,
        "path_fractions": fractions,
        "failures": failures,
    }


def validate_canonical_report(report: dict) -> None:
    failures = []
    source_manifest = report.get("source_manifest", "").replace("\\", "/")
    if not source_manifest.endswith("config/scene_manifest.json"):
        failures.append({"reason": "wrong_source_manifest", "source_manifest": report.get("source_manifest")})
    if report.get("dna", {}).get("represented_bp") != 3954:
        failures.append({"reason": "wrong_dna_bp", "represented_bp": report.get("dna", {}).get("represented_bp")})
    if not math.isclose(float(report.get("dna", {}).get("axis_length_mm", 0.0)), 537.744, abs_tol=0.01):
        failures.append({"reason": "wrong_dna_axis_length", "axis_length_mm": report.get("dna", {}).get("axis_length_mm")})
    mrna_nt = sum(int(segment.get("nt", 0)) for segment in report.get("mrna", {}).get("segments", []))
    if mrna_nt != 1852:
        failures.append({"reason": "wrong_mrna_nt", "mrna_nt": mrna_nt})
    mrna_segment_order = [segment.get("name") for segment in report.get("mrna", {}).get("segments", [])]
    if mrna_segment_order != ["3' UTR", "coding sequence", "5' UTR"]:
        failures.append(
            {
                "reason": "wrong_mrna_order_from_polymerase",
                "expected": ["3' UTR", "coding sequence", "5' UTR"],
                "actual": mrna_segment_order,
            }
        )
    if report.get("mrna", {}).get("path_origin") != "3_prime_at_polymerase_ii":
        failures.append(
            {
                "reason": "wrong_mrna_path_origin",
                "actual": report.get("mrna", {}).get("path_origin"),
            }
        )
    for key in ("mrna", "compact_mrna"):
        rna_report = report.get(key, {})
        if not math.isclose(float(rna_report.get("total_measured_mm", 0.0)), 222.24, rel_tol=0.001, abs_tol=0.02):
            failures.append({"reason": "wrong_rna_contour_length", "rna": key, "measured_mm": rna_report.get("total_measured_mm")})
    mrna_report = report.get("mrna", {})
    if mrna_report.get("surface_mode") != "twisted_groove":
        failures.append({"reason": "wrong_elongated_rna_surface", "actual": mrna_report.get("surface_mode")})
    if mrna_report.get("show_radial_nucleotide_detail") is not False:
        failures.append({"reason": "elongated_rna_radial_detail_enabled"})
    compact_report = report.get("compact_mrna", {})
    if compact_report.get("variant") != "compact_rosette":
        failures.append({"reason": "wrong_compact_rna_variant", "actual": compact_report.get("variant")})
    if compact_report.get("stem_count") != 38 or compact_report.get("paired_stem_bridge_count", 0) <= 0:
        failures.append({"reason": "wrong_compact_rosette_secondary_structure", "stem_count": compact_report.get("stem_count"), "bridges": compact_report.get("paired_stem_bridge_count")})
    expected_detail_keys = {"p53_dna", "nucleosome_loop", "polymerase_rna_start", "ribosome_trna", "actin_product", "cas9_dna"}
    actual_detail_keys = set(report.get("detail_render_specs", {}))
    if actual_detail_keys != expected_detail_keys:
        failures.append({"reason": "wrong_detail_render_set", "actual": sorted(actual_detail_keys), "expected": sorted(expected_detail_keys)})
    label_policy = report.get("render_label_policy", {})
    callout_text = {row.get("text") for row in label_policy.get("overview_label_placement", {}).get("primary_callouts", [])}
    required_callout_text = {spec["text"] for spec in PRIMARY_CALLOUTS.values()}
    if callout_text != required_callout_text:
        failures.append({"reason": "wrong_primary_callouts", "actual": sorted(callout_text), "expected": sorted(required_callout_text)})
    overview_labels = label_policy.get("overview_label_placement", {})
    if len(overview_labels.get("primary_callouts", [])) != 3:
        failures.append({"reason": "wrong_primary_callout_count", "actual": len(overview_labels.get("primary_callouts", []))})
    if overview_labels.get("collision_pairs"):
        failures.append({"reason": "overview_label_collisions", "pairs": overview_labels.get("collision_pairs")})
    if overview_labels.get("molecule_overlap_pairs"):
        failures.append({"reason": "overview_label_overlaps_unrelated_molecule", "pairs": overview_labels.get("molecule_overlap_pairs")})
    if float(overview_labels.get("maximum_displacement", 1.0)) > 0.11:
        failures.append({"reason": "overview_label_too_far_from_molecule", "object": overview_labels.get("maximum_displacement_object"), "maximum_displacement": overview_labels.get("maximum_displacement"), "maximum_allowed": 0.11})
    if len(overview_labels.get("rows", [])) != len(report.get("pdb_assets", [])):
        failures.append({"reason": "missing_pdb_overview_labels", "actual": len(overview_labels.get("rows", [])), "expected": len(report.get("pdb_assets", []))})
    compact_callout = overview_labels.get("compact_callout") or {}
    if compact_callout.get("overlaps_molecule") is not False:
        failures.append({"reason": "compact_mrna_label_overlaps_molecule", "callout": compact_callout})
    if compact_callout.get("offset") != list(COMPACT_CALLOUT["offset"]):
        failures.append({"reason": "wrong_compact_mrna_label_offset", "actual": compact_callout.get("offset"), "expected": list(COMPACT_CALLOUT["offset"])})
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
    failures.extend(report.get("source_curve_visibility", {}).get("failures", []))
    failures.extend(report.get("reader_order_validation", {}).get("failures", []))
    failures.extend(report.get("dna_cluster_validation", {}).get("failures", []))
    failures.extend(report.get("compact_positioning", {}).get("failures", []))
    failures.extend(report.get("trna_ribosome_separation", {}).get("failures", []))
    failures.extend(report.get("contact_validation", {}).get("strict_contact_failures", []))
    scale_bars = report.get("scale_bars", {})
    units = report.get("units", {})
    nm_to_mm = float(units.get("nm_to_mm", 0.4))
    protein_aa_contour_nm = float(units.get("protein_aa_contour_nm", PROTEIN_AA_CONTOUR_NM))
    required_scale_bars = {
        "scale_dna_100_bp": 100.0 * float(units.get("dna_bp_rise_nm", 0.34)) * nm_to_mm,
        "scale_rna_100_nt": 100.0 * float(units.get("mrna_nt_contour_nm", 0.3)) * nm_to_mm,
        "scale_protein_33_aa": 33.0 * protein_aa_contour_nm * nm_to_mm,
        "scale_10_nm": NANOMETER_SCALE_BAR_NM * nm_to_mm,
    }
    if set(scale_bars) != set(required_scale_bars):
        failures.append(
            {
                "reason": "wrong_canonical_scale_bars",
                "actual": sorted(scale_bars),
                "expected": sorted(required_scale_bars),
            }
        )
    actual_scale_bar_order = list(scale_bars)
    expected_scale_bar_order = list(required_scale_bars)
    if actual_scale_bar_order != expected_scale_bar_order:
        failures.append(
            {
                "reason": "wrong_scale_bar_display_order",
                "actual": actual_scale_bar_order,
                "expected": expected_scale_bar_order,
            }
        )
    for name, expected_length in required_scale_bars.items():
        actual_length = (scale_bars.get(name) or {}).get("length_mm")
        if actual_length is None or not math.isclose(float(actual_length), expected_length, abs_tol=1e-6):
            failures.append(
                {
                    "reason": "wrong_scale_bar_length",
                    "name": name,
                    "actual_length_mm": actual_length,
                    "expected_length_mm": expected_length,
                }
            )
    asset_reports = {item.get("name"): item for item in report.get("pdb_assets", [])}
    requested_fractions = {
        "Pumilio RBP": 0.20,
        "Poly(A)-binding RBP": 0.24,
        "Argonaute": 0.54,
        "HuR-like RBP": 0.60,
        "Ribosome large subunit": 0.89,
        "Transcription factor 3": 0.115,
        "Nucleosome": 0.165,
    }
    fractions = {}
    for name, expected_fraction in requested_fractions.items():
        attachment = (asset_reports.get(name) or {}).get("attachment_empty") or {}
        fraction = attachment.get("fraction")
        fractions[name] = fraction
        if fraction is None or not math.isclose(float(fraction), expected_fraction, abs_tol=1e-6):
            failures.append({"reason": "wrong_requested_path_anchor", "name": name, "actual": fraction, "expected": expected_fraction})
    pabp = fractions.get("Poly(A)-binding RBP")
    pumilio = fractions.get("Pumilio RBP")
    hur = fractions.get("HuR-like RBP")
    argonaute = fractions.get("Argonaute")
    ribosome = fractions.get("Ribosome large subunit")
    foxm1 = fractions.get("Transcription factor 3")
    nucleosome = fractions.get("Nucleosome")
    promoter_fraction = 500.0 / 3954.0
    if pabp is not None and pumilio is not None and not float(pumilio) < float(pabp):
        failures.append({"reason": "pabp_not_after_pumilio", "pumilio": pumilio, "pabp": pabp})
    if hur is not None and argonaute is not None and ribosome is not None and not float(argonaute) < float(hur) < float(ribosome):
        failures.append({"reason": "hur_not_between_argonaute_and_ribosome", "argonaute": argonaute, "hur": hur, "ribosome": ribosome})
    if foxm1 is not None and not float(foxm1) < promoter_fraction:
        failures.append({"reason": "foxm1_outside_promoter", "fraction": foxm1, "promoter_end_fraction": promoter_fraction})
    if foxm1 is not None and nucleosome is not None and abs(float(nucleosome) - float(foxm1)) < 0.04:
        failures.append({"reason": "foxm1_too_close_to_nucleosome", "foxm1": foxm1, "nucleosome": nucleosome, "minimum_fraction_gap": 0.04})
    if failures:
        raise RuntimeError(f"Canonical Canonical validation failed: {failures}")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = canonical.write_manifest()

    scene.clean_scene()
    scene.configure_scene()
    beauty_rendering = configure_canonical_beauty_render()
    collections = base.make_collections(manifest["collections"])
    materials = base.make_materials()
    scene.soften_materials(materials)
    color_palette = apply_canonical_color_palette(materials)

    dna = build_canonical_dna(manifest, collections, materials)
    mrna = build_canonical_mrna(manifest, collections, materials)
    actin_positioning = place_actin_from_mrna_endpoint(manifest, mrna["path"])
    compact_manifest, compact_positioning = compact_positioned_manifest(manifest, mrna["path"])
    compact_mrna = build_canonical_compact_mrna(compact_manifest, collections, materials)
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
        "canonical_version": "canonical",
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
        "color_palette": color_palette,
        "beauty_rendering": beauty_rendering,
        "layout_intent": manifest.get("layout_intent", {}),
        "actin_positioning": actin_positioning,
    }
    report["dna"] = dna["report"]
    report["mrna"] = mrna["report"]
    report["compact_mrna"] = compact_mrna["report"]
    report["source_curve_visibility"] = source_curve_visibility_validation(report)
    report["reader_order_validation"] = reader_order_validation(dna["path"])
    report["compact_positioning"] = compact_position_validation(compact_positioning, compact_mrna["center_mm"])
    report["scale_bars"] = add_canonical_scale_bars(manifest, collections, materials)
    report["pdb_assets"] = scene.build_assets(manifest, collections, materials, path_context)
    report["trna_ribosome_separation"] = trna_ribosome_separation_report(report["pdb_assets"])
    report["attachments"] = contact_helpers.attachment_table(report["pdb_assets"])
    report["dna_cluster_validation"] = dna_cluster_validation(report["pdb_assets"])
    report["polymerase_rna_origin"] = contact_helpers.polymerase_rna_origin_report(report["pdb_assets"], mrna["path"])
    report["contact_validation"] = contact_helpers.contact_validation(report["pdb_assets"])
    final_contact = contact_helpers.final_surface_contact_validation(
        report["pdb_assets"],
        path_context,
        contact_radii,
        {"dna": dna["report"], "mrna": mrna["report"]},
    )
    report["contact_validation"]["final_surface_contact"] = final_contact
    report["contact_validation"]["strict_contact_failures"].extend(final_contact["failures"])
    camera_names = add_canonical_cameras(
        manifest,
        dna["path"],
        mrna["path"],
        compact_mrna["path"],
        report["pdb_assets"],
        collections,
        materials,
    )
    report["detail_context"] = add_detail_context_curves(mrna["path"], collections, materials)
    for obj in bpy.data.objects:
        if obj.type == "FONT":
            obj.hide_render = True
    full_overview_camera_fit = fit_camera_to_renderables(camera_names["full_overview"])
    label_orientation = orient_labels_to_camera(camera_names["full_overview"])
    overview_label_placement = place_overview_labels(camera_names["full_overview"], collections, materials)
    report["beauty_rendering"]["environment"] = add_canonical_beauty_environment(camera_names["full_overview"], collections)
    report["beauty_rendering"]["surface_polish"] = polish_canonical_surface_rendering()
    report["render_label_policy"] = {
        "hidden_non_scale_label_count": 0,
        "scale_bar_labels_visible": True,
        "pdb_asset_labels_visible": True,
        "label_orientation": label_orientation,
        "overview_label_placement": overview_label_placement,
        "full_overview_camera_fit": full_overview_camera_fit,
    }
    report["detail_render_specs"] = {
        key: {
            "camera": camera_names[key],
            "title": DETAIL_TITLES[key],
            "output": str(DETAIL_PREVIEWS[key]),
            "object_patterns": list(FOCUS_OBJECT_PATTERNS[key]),
            "margin_fraction": 0.12,
            "shared_ortho_scale_mm": DETAIL_SHARED_ORTHO_SCALE_MM,
        }
        for key in DETAIL_TITLES
    }
    scene.validate_scene(report)
    validate_canonical_report(report)

    overview_only = os.environ.get("CANONICAL_OVERVIEW_ONLY", "") == "1"
    primary_camera = bpy.data.objects[camera_names["full_overview"]]
    bpy.context.scene.camera = primary_camera
    bpy.context.scene.render.filepath = str(DETAIL_PREVIEWS["full_overview"])
    bpy.ops.render.render(write_still=True)
    shutil.copyfile(DETAIL_PREVIEWS["full_overview"], PREVIEW_PATH)

    overview_visibility = snapshot_render_visibility()
    detail_runs = {}
    if overview_only and REPORT_PATH.exists():
        detail_runs = json.loads(REPORT_PATH.read_text(encoding="utf-8")).get("detail_rendering", {})
    for key in (() if overview_only else DETAIL_TITLES):
        restore_render_visibility(overview_visibility)
        visibility = apply_focus_visibility(key)
        camera = bpy.data.objects[camera_names[key]]
        bpy.context.scene.camera = camera
        camera_fit = fit_camera_to_renderables(camera_names[key], margin_fraction=0.12)
        fitted_scale = float(camera.data.ortho_scale)
        if fitted_scale > DETAIL_SHARED_ORTHO_SCALE_MM + 1e-6:
            raise RuntimeError(
                f"Detail target {key} needs {fitted_scale:.3f} mm, exceeding shared "
                f"orthographic scale {DETAIL_SHARED_ORTHO_SCALE_MM:.3f} mm"
            )
        camera.data.ortho_scale = DETAIL_SHARED_ORTHO_SCALE_MM
        camera_fit["auto_fit_ortho_scale_mm"] = fitted_scale
        camera_fit["new_ortho_scale_mm"] = DETAIL_SHARED_ORTHO_SCALE_MM
        camera_fit["scale_policy"] = "identical_orthographic_scale_for_cross_panel_size_comparison"
        environment = add_canonical_beauty_environment(camera_names[key], collections)
        title_object = create_detail_title(key, camera_names[key], collections, materials)
        bpy.context.scene.render.filepath = str(DETAIL_PREVIEWS[key])
        bpy.ops.render.render(write_still=True)
        detail_runs[key] = {
            "visibility": visibility,
            "camera_fit": camera_fit,
            "environment": environment,
            "title_object": title_object,
        }

    restore_render_visibility(overview_visibility)
    for obj in list(bpy.data.objects):
        if obj.get("canonical_detail_title"):
            bpy.data.objects.remove(obj, do_unlink=True)
    report["detail_rendering"] = detail_runs
    bpy.context.scene.camera = primary_camera
    report["beauty_rendering"]["environment"] = add_canonical_beauty_environment(camera_names["full_overview"], collections)
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_PATH))
    print(f"Wrote {BLEND_PATH}")
    print(f"Wrote {PREVIEW_PATH}")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
