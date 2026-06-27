#!/usr/bin/env python3
"""Build the canonical V5 reader-order gene-expression scene."""

from __future__ import annotations

import copy
import json
import math
import sys
from pathlib import Path

import bpy
from mathutils import Vector


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import blender_nucleic_meshes as direct_nucleic_meshes  # noqa: E402
import build_gene_expression_scene as base  # noqa: E402
import build_gene_expression_surface_scene as scene  # noqa: E402
import build_gene_expression_surface_scene_v4 as v4_builder  # noqa: E402
import canonical_v5_config as v5  # noqa: E402
import procedural_nucleic_geometry as nucleic_geometry  # noqa: E402


OUTPUT_DIR = v5.OUTPUT_DIR
BLEND_PATH = v5.BLEND_PATH
PREVIEW_PATH = v5.PREVIEW_PATH
REPORT_PATH = v5.REPORT_PATH
DETAIL_PREVIEWS = {
    "full_overview": OUTPUT_DIR / "preview_gene_expression_surface_style_v5_full_overview.png",
    "polymerase_rna_start": OUTPUT_DIR / "preview_gene_expression_surface_style_v5_polymerase_rna_start.png",
    "nucleosome_loop": OUTPUT_DIR / "preview_gene_expression_surface_style_v5_nucleosome_loop.png",
    "mrna_spiral": OUTPUT_DIR / "preview_gene_expression_surface_style_v5_mrna_spiral.png",
    "ribosome_top_translation": OUTPUT_DIR / "preview_gene_expression_surface_style_v5_ribosome_top_translation.png",
    "compact_mrna": OUTPUT_DIR / "preview_gene_expression_surface_style_v5_compact_mrna.png",
}
PROTEIN_AA_CONTOUR_NM = v5.PROTEIN_AA_CONTOUR_NM
NANOMETER_SCALE_BAR_NM = 10.0
V5_MATERIAL_COLORS = {
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
V5_WORLD_COLOR = (0.94, 0.955, 0.97, 1.0)
V5_BACKDROP_TOP_LEFT = (0.86, 0.925, 0.975, 1.0)
V5_BACKDROP_TOP_RIGHT = (0.96, 0.985, 0.995, 1.0)
V5_BACKDROP_BOTTOM_LEFT = (0.985, 0.945, 0.82, 1.0)
V5_BACKDROP_BOTTOM_RIGHT = (1.0, 0.885, 0.74, 1.0)
V5_LABEL_MATERIALS = {"black", "label_grey", "scale_grey"}
V5_CYCLES_SAMPLES = 64
OVERVIEW_LABEL_POSITIONS = {
    "label_scene_scale": (0.17, 0.50),
    "label_DNA_v5": (0.38, 0.34),
    "label_Transcription factor 4": (0.36, 0.39),
    "label_Cas9": (0.46, 0.38),
    "label_Transcription factor 1": (0.48, 0.45),
    "label_MS2 coat protein MCP": (0.49, 0.30),
    "label_mCherry/RFP tag": (0.61, 0.29),
    "label_p53 tetramer bound to DNA": (0.61, 0.48),
    "label_Pumilio RBP": (0.52, 0.50),
    "label_RNA polymerase II elongation complex": (0.66, 0.33),
    "label_mRNA_v5": (0.70, 0.46),
    "label_Nucleosome": (0.77, 0.50),
    "label_Transcription factor 3": (0.79, 0.44),
    "label_Argonaute": (0.73, 0.55),
    "label_Poly(A)-binding RBP": (0.60, 0.69),
    "label_HuR-like RBP": (0.73, 0.67),
    "label_Ribosome large subunit": (0.64, 0.73),
    "label_Ribosome small subunit": (0.71, 0.73),
    "label_Standalone tRNA": (0.78, 0.71),
    "label_compact_mrna_v5": (0.72, 0.86),
    "label_Actin protein": (0.69, 0.92),
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


def build_v5_dna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = direct_nucleic_meshes.build_dna_meshes(manifest, collections, materials)
    report["source_curve"] = add_source_path_curve("V5_DNA_source_path", path, collections["DNA"], materials["black"])
    start = path.points[0]
    base.create_text(
        "label_DNA_v5",
        "ACTB promoter + gene DNA\n3954 bp",
        (start.x + 4.0, start.y - 3.0, start.z + 0.25),
        1.55,
        materials["black"],
        collections["Labels"],
    )
    return {"path": path, "report": report, "features": {"nucleosome_loop": report.get("nucleosome_loop")}}


def build_v5_mrna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = direct_nucleic_meshes.build_mrna_meshes(manifest, collections, materials)
    report["source_curve"] = add_source_path_curve("V5_mRNA_source_path", path, collections["mRNA"], materials["black"])
    start = path.points[0]
    base.create_text(
        "label_mRNA_v5",
        "actin mRNA\n1852 nt",
        (start.x + 3.0, start.y + 3.0, start.z + 0.3),
        1.55,
        materials["black"],
        collections["Labels"],
    )
    return {"path": path, "report": report}


def compact_positioned_manifest(manifest: dict, mrna_path) -> tuple[dict, dict]:
    adjusted = copy.deepcopy(manifest)
    mrna_end = Vector((mrna_path.points[-1].x, mrna_path.points[-1].y, mrna_path.points[-1].z))
    actin = v5.asset_by_name(adjusted, "Actin protein")
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


def build_v5_compact_mrna(manifest: dict, collections: dict, materials: dict) -> dict:
    path, report = direct_nucleic_meshes.build_compact_mrna_meshes(manifest, collections, materials)
    compact_center = v4_builder.path_center(path)
    base.create_text(
        "label_compact_mrna_v5",
        "compact mRNP-like\nmRNA reference",
        (compact_center.x, compact_center.y + 3.2, compact_center.z + 0.35),
        1.45,
        materials["label_grey"],
        collections["Labels"],
    )
    return {"path": path, "report": report, "center_mm": [compact_center.x, compact_center.y, compact_center.z]}


def add_v5_scale_bars(manifest: dict, collections: dict, materials: dict) -> dict:
    units = manifest["units"]
    nm_to_mm = float(units["nm_to_mm"])
    protein_aa_contour_nm = float(units.get("protein_aa_contour_nm", PROTEIN_AA_CONTOUR_NM))
    bars = [
        {
            "name": "scale_rna_100_nt",
            "label": "RNA 100 nt",
            "length_nm": 100.0 * float(units["mrna_nt_contour_nm"]),
            "origin": (-82.0, -2.0, 0.0),
            "quantity": 100,
            "unit": "nt",
            "basis": "single-stranded RNA contour length",
        },
        {
            "name": "scale_dna_100_bp",
            "label": "DNA 100 bp",
            "length_nm": 100.0 * float(units["dna_bp_rise_nm"]),
            "origin": (-82.0, -7.0, 0.0),
            "quantity": 100,
            "unit": "bp",
            "basis": "B-form DNA axial rise",
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
    return v4_builder.create_camera(name, location, target, ortho_scale)


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
        label_material = material.name.split(".")[0] in V5_LABEL_MATERIALS
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


def configure_v5_beauty_render() -> dict:
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
        cycles.samples = V5_CYCLES_SAMPLES
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
    world.color = V5_WORLD_COLOR[:3]
    world.use_nodes = True
    background = world.node_tree.nodes.get("Background")
    if background:
        background.inputs["Color"].default_value = V5_WORLD_COLOR
        background.inputs["Strength"].default_value = 0.52 if scene_data.render.engine == "CYCLES" else 0.48
    return {
        "requested_engine": requested_engine,
        "engine": scene_data.render.engine,
        "engine_fallback_reason": engine_fallback_reason,
        "resolution": [scene_data.render.resolution_x, scene_data.render.resolution_y],
        "world_color": list(V5_WORLD_COLOR),
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
        if obj.get("v5_beauty_backdrop") and not include_backdrop:
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


def v5_beauty_collection(collections: dict[str, bpy.types.Collection]) -> bpy.types.Collection:
    collection = collections.get("Beauty") or bpy.data.collections.get("Beauty")
    if collection is None:
        collection = bpy.data.collections.new("Beauty")
    if collection.name not in {child.name for child in bpy.context.scene.collection.children}:
        bpy.context.scene.collection.children.link(collection)
    collections["Beauty"] = collection
    return collection


def remove_existing_v5_beauty_objects() -> None:
    for obj in list(bpy.data.objects):
        if obj.name == "large_softbox" or obj.get("v5_beauty_object"):
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


def v5_backdrop_color(x_norm: float, y_norm: float) -> tuple[float, float, float, float]:
    bottom = blend_rgba(V5_BACKDROP_BOTTOM_LEFT, V5_BACKDROP_BOTTOM_RIGHT, x_norm)
    top = blend_rgba(V5_BACKDROP_TOP_LEFT, V5_BACKDROP_TOP_RIGHT, x_norm)
    return blend_rgba(bottom, top, y_norm)


def create_v5_backdrop(camera_name: str, collections: dict[str, bpy.types.Collection]) -> dict:
    camera = bpy.data.objects[camera_name]
    bounds = camera_space_renderable_bounds(camera_name)
    if not bounds.get("available"):
        return {"created": False, "reason": "no_renderable_objects"}

    collection = v5_beauty_collection(collections)
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

    mesh = bpy.data.meshes.new("v5_camera_gradient_backdrop_mesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    obj = bpy.data.objects.new("v5_camera_gradient_backdrop", mesh)
    obj["v5_beauty_object"] = True
    obj["v5_beauty_backdrop"] = True
    obj["camera_aligned_to"] = camera_name
    obj.hide_select = True
    collection.objects.link(obj)

    face_index = 0
    for row in range(rows):
        for column in range(columns):
            x_norm = (column + 0.5) / columns
            y_norm = (row + 0.5) / rows
            mat = backdrop_material(f"v5_backdrop_gradient_{row:02d}_{column:02d}", v5_backdrop_color(x_norm, y_norm))
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
            "top_left": list(V5_BACKDROP_TOP_LEFT),
            "top_right": list(V5_BACKDROP_TOP_RIGHT),
            "bottom_left": list(V5_BACKDROP_BOTTOM_LEFT),
            "bottom_right": list(V5_BACKDROP_BOTTOM_RIGHT),
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
    light["v5_beauty_object"] = True
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


def add_v5_beauty_environment(camera_name: str, collections: dict[str, bpy.types.Collection]) -> dict:
    remove_existing_v5_beauty_objects()
    collection = v5_beauty_collection(collections)
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
            "v5_key_softbox",
            collection,
            camera,
            scene_center_world,
            (-0.34 * view_width, 0.32 * view_height, front_z),
            10500.0,
            74.0,
            (1.0, 0.94, 0.84),
        ),
        create_area_light(
            "v5_fill_softbox",
            collection,
            camera,
            scene_center_world,
            (0.42 * view_width, -0.14 * view_height, front_z + 14.0),
            2600.0,
            210.0,
            (0.84, 0.91, 1.0),
        ),
        create_area_light(
            "v5_rim_softbox",
            collection,
            camera,
            scene_center_world,
            (0.26 * view_width, 0.48 * view_height, rim_z),
            6200.0,
            56.0,
            (0.78, 0.86, 1.0),
        ),
        create_area_light(
            "v5_warm_wash",
            collection,
            camera,
            scene_center_world,
            (-0.08 * view_width, -0.46 * view_height, front_z + 26.0),
            2600.0,
            210.0,
            (1.0, 0.88, 0.74),
        ),
    ]
    backdrop = create_v5_backdrop(camera_name, collections)
    bpy.context.view_layer.update()
    return {
        "created": True,
        "collection": collection.name,
        "camera": camera_name,
        "scene_center_world_mm": [scene_center_world.x, scene_center_world.y, scene_center_world.z],
        "lights": lights,
        "backdrop": backdrop,
    }


def polish_v5_surface_rendering() -> dict:
    smoothed_meshes = []
    weighted_normal_meshes = []
    for obj in bpy.data.objects:
        if obj.type != "MESH" or obj.hide_render or obj.get("v5_beauty_backdrop"):
            continue
        if obj.data and obj.data.polygons:
            for polygon in obj.data.polygons:
                polygon.use_smooth = True
            smoothed_meshes.append(obj.name)
        if obj.modifiers.get("v5_weighted_normals") is None:
            modifier = obj.modifiers.new("v5_weighted_normals", "WEIGHTED_NORMAL")
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


def apply_v5_color_palette(materials: dict[str, bpy.types.Material]) -> dict:
    applied = {}
    for name, color in V5_MATERIAL_COLORS.items():
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


def orient_labels_to_camera(camera_name: str) -> dict:
    camera = bpy.data.objects[camera_name]
    rows = []
    for obj in bpy.data.objects:
        if obj.type != "FONT":
            continue
        obj.hide_render = False
        obj.rotation_euler = camera.rotation_euler
        obj["oriented_to_camera"] = camera_name
        rows.append({"object": obj.name, "text": obj.data.body})
    bpy.context.view_layer.update()
    return {
        "policy": "all V5 text labels are visible and billboarded to the full-overview camera",
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


def place_overview_labels(camera_name: str) -> dict:
    camera = bpy.data.objects[camera_name]
    scene_render = bpy.context.scene.render
    aspect = scene_render.resolution_x / max(float(scene_render.resolution_y), 1.0)
    inv = camera.matrix_world.inverted()
    rows = []
    for object_name, (x_norm, y_norm) in OVERVIEW_LABEL_POSITIONS.items():
        obj = bpy.data.objects.get(object_name)
        if obj is None:
            continue
        local = inv @ obj.location
        if aspect >= 1.0:
            local.x = (x_norm - 0.5) * float(camera.data.ortho_scale)
            local.y = (y_norm - 0.5) * float(camera.data.ortho_scale) / aspect
        else:
            local.x = (x_norm - 0.5) * float(camera.data.ortho_scale) * aspect
            local.y = (y_norm - 0.5) * float(camera.data.ortho_scale)
        obj.location = camera.matrix_world @ local
        obj["overview_label_position"] = [x_norm, y_norm]
        rows.append({"object": object_name, "view_position": [x_norm, y_norm]})
    bpy.context.view_layer.update()
    return {
        "policy": "selected overview labels are spread in camera space after fitting the full-scene camera",
        "camera": camera_name,
        "placed_label_count": len(rows),
        "rows": rows,
    }


def add_v5_cameras(manifest: dict, dna_path, mrna_path, compact_path, collections: dict, materials: dict) -> dict[str, str]:
    base.add_lighting_and_camera(collections, materials)
    dna_center = v4_builder.path_center(dna_path)
    mrna_center = v4_builder.path_center(mrna_path)
    compact_center = v4_builder.path_center(compact_path)
    ribosome_point = mrna_path.point_at_length(mrna_path.length * 0.892)
    actin_point = mrna_path.point_at_length(mrna_path.length)
    loop_report = nucleic_geometry.dna_nucleosome_loop_report(manifest) or {}
    loop_center = loop_report.get("center_mm") or [0.0, -9.0, 0.0]

    camera_names = {
        "full_overview": create_camera(
            "Camera_v5_full_overview",
            (dna_center.x - 48.0, dna_center.y - 142.0, 128.0),
            Vector((dna_center.x + 12.0, dna_center.y - 6.0, 36.0)),
            190.0,
        ),
        "polymerase_rna_start": create_camera(
            "Camera_v5_polymerase_rna_start",
            (mrna_path.points[0].x - 24.0, mrna_path.points[0].y - 38.0, mrna_path.points[0].z + 30.0),
            Vector((mrna_path.points[0].x, mrna_path.points[0].y, mrna_path.points[0].z + 4.0)),
            26.0,
        ),
        "nucleosome_loop": create_camera(
            "Camera_v5_nucleosome_loop",
            (loop_center[0] - 8.0, loop_center[1] - 12.0, (loop_center[2] if len(loop_center) > 2 else 0.0) + 9.0),
            Vector((loop_center[0], loop_center[1], loop_center[2] if len(loop_center) > 2 else 0.0)),
            8.0,
        ),
        "mrna_spiral": create_camera(
            "Camera_v5_mrna_spiral",
            (mrna_center.x - 48.0, mrna_center.y - 72.0, mrna_center.z + 58.0),
            mrna_center,
            92.0,
        ),
        "ribosome_top_translation": create_camera(
            "Camera_v5_ribosome_top_translation",
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
            "Camera_v5_compact_mrna",
            (compact_center.x - 12.0, compact_center.y - 30.0, compact_center.z + 32.0),
            compact_center,
            24.0,
        ),
    }
    bpy.context.scene.camera = bpy.data.objects[camera_names["full_overview"]]
    return camera_names


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
    for path_name, object_name in (("dna", "V5_DNA_source_path"), ("mrna", "V5_mRNA_source_path")):
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
        "path_mode": "v5_reader_order_serpentine_with_nucleosome_loop",
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
    for name in sorted(v5.DNA_STRICT_CONTACT_NAMES):
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
    return {"policy": "all DNA regulatory binders cluster in the first 20.5 percent of V5 DNA", "rows": rows, "failures": failures}


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


def validate_v5_report(report: dict) -> None:
    failures = []
    source_manifest = report.get("source_manifest", "").replace("\\", "/")
    if not source_manifest.endswith("config/scene_manifest_v5.json"):
        failures.append({"reason": "wrong_source_manifest", "source_manifest": report.get("source_manifest")})
    if report.get("dna", {}).get("represented_bp") != 3954:
        failures.append({"reason": "wrong_dna_bp", "represented_bp": report.get("dna", {}).get("represented_bp")})
    if not math.isclose(float(report.get("dna", {}).get("axis_length_mm", 0.0)), 537.744, abs_tol=0.01):
        failures.append({"reason": "wrong_dna_axis_length", "axis_length_mm": report.get("dna", {}).get("axis_length_mm")})
    mrna_nt = sum(int(segment.get("nt", 0)) for segment in report.get("mrna", {}).get("segments", []))
    if mrna_nt != 1852:
        failures.append({"reason": "wrong_mrna_nt", "mrna_nt": mrna_nt})
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
    failures.extend(report.get("contact_validation", {}).get("strict_contact_failures", []))
    scale_bars = report.get("scale_bars", {})
    units = report.get("units", {})
    nm_to_mm = float(units.get("nm_to_mm", 0.4))
    protein_aa_contour_nm = float(units.get("protein_aa_contour_nm", PROTEIN_AA_CONTOUR_NM))
    required_scale_bars = {
        "scale_rna_100_nt": 100.0 * float(units.get("mrna_nt_contour_nm", 0.3)) * nm_to_mm,
        "scale_dna_100_bp": 100.0 * float(units.get("dna_bp_rise_nm", 0.34)) * nm_to_mm,
        "scale_protein_33_aa": 33.0 * protein_aa_contour_nm * nm_to_mm,
        "scale_10_nm": NANOMETER_SCALE_BAR_NM * nm_to_mm,
    }
    if set(scale_bars) != set(required_scale_bars):
        failures.append(
            {
                "reason": "wrong_v5_scale_bars",
                "actual": sorted(scale_bars),
                "expected": sorted(required_scale_bars),
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
    if failures:
        raise RuntimeError(f"Canonical V5 validation failed: {failures}")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if not v5.MANIFEST_PATH.exists():
        v5.write_manifest()
    manifest = json.loads(v5.MANIFEST_PATH.read_text(encoding="utf-8"))

    scene.clean_scene()
    scene.configure_scene()
    beauty_rendering = configure_v5_beauty_render()
    collections = base.make_collections(manifest["collections"])
    materials = base.make_materials()
    scene.soften_materials(materials)
    color_palette = apply_v5_color_palette(materials)

    dna = build_v5_dna(manifest, collections, materials)
    mrna = build_v5_mrna(manifest, collections, materials)
    compact_manifest, compact_positioning = compact_positioned_manifest(manifest, mrna["path"])
    compact_mrna = build_v5_compact_mrna(compact_manifest, collections, materials)
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
        "kind": v5.REPORT_KIND,
        "canonical_version": "v5",
        "units": manifest["units"],
        "source_manifest": str(v5.MANIFEST_PATH),
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
    }
    report["dna"] = dna["report"]
    report["mrna"] = mrna["report"]
    report["compact_mrna"] = compact_mrna["report"]
    report["source_curve_visibility"] = source_curve_visibility_validation(report)
    report["reader_order_validation"] = reader_order_validation(dna["path"])
    report["compact_positioning"] = compact_position_validation(compact_positioning, compact_mrna["center_mm"])
    report["scale_bars"] = add_v5_scale_bars(manifest, collections, materials)
    report["pdb_assets"] = scene.build_assets(manifest, collections, materials, path_context)
    report["attachments"] = v4_builder.attachment_table(report["pdb_assets"])
    report["dna_cluster_validation"] = dna_cluster_validation(report["pdb_assets"])
    report["polymerase_rna_origin"] = v4_builder.polymerase_rna_origin_report(report["pdb_assets"], mrna["path"])
    report["contact_validation"] = v4_builder.contact_validation(report["pdb_assets"])
    final_contact = v4_builder.final_surface_contact_validation(
        report["pdb_assets"],
        path_context,
        contact_radii,
        {"dna": dna["report"], "mrna": mrna["report"]},
    )
    report["contact_validation"]["final_surface_contact"] = final_contact
    report["contact_validation"]["strict_contact_failures"].extend(final_contact["failures"])
    camera_names = add_v5_cameras(manifest, dna["path"], mrna["path"], compact_mrna["path"], collections, materials)
    label_orientation = orient_labels_to_camera(camera_names["full_overview"])
    full_overview_camera_fit = fit_camera_to_renderables(camera_names["full_overview"])
    overview_label_placement = place_overview_labels(camera_names["full_overview"])
    report["beauty_rendering"]["environment"] = add_v5_beauty_environment(camera_names["full_overview"], collections)
    report["beauty_rendering"]["surface_polish"] = polish_v5_surface_rendering()
    report["render_label_policy"] = {
        "hidden_non_scale_label_count": 0,
        "scale_bar_labels_visible": True,
        "pdb_asset_labels_visible": True,
        "label_orientation": label_orientation,
        "overview_label_placement": overview_label_placement,
        "full_overview_camera_fit": full_overview_camera_fit,
    }
    scene.validate_scene(report)
    validate_v5_report(report)

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
    main()
