#!/usr/bin/env python3
"""Build and optionally render the V5 educational flythrough animation."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import bpy
from mathutils import Vector


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE_BLEND = ROOT / "outputs" / "canonical" / "gene_expression_surface_style_v5.blend"
DEFAULT_REPORT = ROOT / "outputs" / "canonical" / "gene_expression_surface_scene_v5_report.json"
DEFAULT_OUTPUT_DIR = ROOT / "experiments" / "v5_flythrough_animation" / "outputs"


def parse_args() -> argparse.Namespace:
    argv = sys.argv
    if "--" in argv:
        argv = argv[argv.index("--") + 1 :]
    else:
        argv = []
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-blend", type=Path, default=DEFAULT_SOURCE_BLEND)
    parser.add_argument("--source-report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--output-blend", type=Path)
    parser.add_argument("--output-mp4", type=Path)
    parser.add_argument("--duration-seconds", type=float, default=60.0)
    parser.add_argument("--fps", type=int, default=24)
    parser.add_argument("--resolution-x", type=int, default=1920)
    parser.add_argument("--resolution-y", type=int, default=1080)
    parser.add_argument("--smoke-test", action="store_true")
    parser.add_argument("--skip-video-render", action="store_true")
    return parser.parse_args(argv)


def as_path(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def scene_collection(name: str) -> bpy.types.Collection:
    existing = bpy.data.collections.get(name)
    if existing:
        return existing
    collection = bpy.data.collections.new(name)
    bpy.context.scene.collection.children.link(collection)
    return collection


def link_to_collection(obj: bpy.types.Object, collection: bpy.types.Collection) -> None:
    collection.objects.link(obj)
    for user_collection in list(obj.users_collection):
        if user_collection != collection:
            user_collection.objects.unlink(obj)


def material(name: str, color: tuple[float, float, float, float], emission_strength: float = 0.0) -> bpy.types.Material:
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = color
    mat.use_nodes = True
    mat.blend_method = "BLEND"
    mat.show_transparent_back = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf:
        if "Base Color" in bsdf.inputs:
            bsdf.inputs["Base Color"].default_value = color
        if "Alpha" in bsdf.inputs:
            bsdf.inputs["Alpha"].default_value = color[3]
        if "Emission Color" in bsdf.inputs:
            bsdf.inputs["Emission Color"].default_value = color
        elif "Emission" in bsdf.inputs:
            bsdf.inputs["Emission"].default_value = color
        if "Emission Strength" in bsdf.inputs:
            bsdf.inputs["Emission Strength"].default_value = emission_strength
        if "Roughness" in bsdf.inputs:
            bsdf.inputs["Roughness"].default_value = 0.6
    return mat


def configure_view_settings() -> None:
    scene = bpy.context.scene
    for transform in ("AgX", "Filmic"):
        try:
            scene.view_settings.view_transform = transform
            break
        except TypeError:
            continue
    for look in ("Medium High Contrast", "High Contrast", "None"):
        try:
            scene.view_settings.look = look
            break
        except TypeError:
            continue
    scene.view_settings.exposure = 0.0
    scene.view_settings.gamma = 1.0
    if scene.world is None:
        scene.world = bpy.data.worlds.new("Animation_Cinematic_World")
    scene.world.color = (0.012, 0.014, 0.021)
    scene.world.use_nodes = True
    background = scene.world.node_tree.nodes.get("Background") if scene.world.node_tree else None
    if background:
        if "Color" in background.inputs:
            background.inputs["Color"].default_value = (0.012, 0.014, 0.021, 1.0)
        if "Strength" in background.inputs:
            background.inputs["Strength"].default_value = 0.055
    if hasattr(scene.world, "mist_settings"):
        scene.world.mist_settings.use_mist = False


def set_node_input(node: bpy.types.Node, names: tuple[str, ...], value) -> None:
    for socket in node.inputs:
        if socket.name in names and hasattr(socket, "default_value"):
            try:
                socket.default_value = value
                return
            except Exception:
                continue


def set_node_property(node: bpy.types.Node, name: str, value) -> None:
    if hasattr(node, name):
        try:
            setattr(node, name, value)
        except Exception:
            return


def configure_compositor() -> list[str]:
    scene = bpy.context.scene
    if hasattr(scene, "use_nodes"):
        scene.use_nodes = True
    if hasattr(scene.render, "use_compositing"):
        scene.render.use_compositing = True
    if hasattr(scene.render, "use_sequencer"):
        scene.render.use_sequencer = False

    if hasattr(scene, "compositing_node_group"):
        tree = bpy.data.node_groups.new("Animation_Cinematic_Compositor", "CompositorNodeTree")
        scene.compositing_node_group = tree
        output_node_type = "NodeGroupOutput"
        tree.interface.new_socket(name="Image", in_out="OUTPUT", socket_type="NodeSocketColor")
    else:
        tree = scene.node_tree
        tree.nodes.clear()
        output_node_type = "CompositorNodeComposite"

    render_layers = tree.nodes.new("CompositorNodeRLayers")
    render_layers.name = "Animation_Compositor_Render_Layers"
    glare = tree.nodes.new("CompositorNodeGlare")
    glare.name = "Animation_Compositor_Fog_Glow"
    set_node_property(glare, "glare_type", "FOG_GLOW")
    set_node_property(glare, "quality", "MEDIUM")
    set_node_property(glare, "threshold", 1.15)
    set_node_property(glare, "size", 6)
    set_node_property(glare, "mix", -0.72)
    set_node_input(glare, ("Type",), "Fog Glow")
    set_node_input(glare, ("Quality",), "Medium")
    set_node_input(glare, ("Threshold",), 1.15)
    set_node_input(glare, ("Strength",), 0.18)
    set_node_input(glare, ("Saturation",), 0.82)
    set_node_input(glare, ("Size",), 0.38)

    color_balance = tree.nodes.new("CompositorNodeColorBalance")
    color_balance.name = "Animation_Compositor_Cool_Shadows_Warm_Highlights"
    set_node_input(color_balance, ("Factor",), 0.62)
    set_node_input(color_balance, ("Type",), "Lift/Gamma/Gain")
    set_node_input(color_balance, ("Lift",), (0.94, 0.97, 1.03, 1.0))
    set_node_input(color_balance, ("Gamma",), (1.0, 1.0, 1.0, 1.0))
    set_node_input(color_balance, ("Gain",), (1.04, 1.02, 0.97, 1.0))
    try:
        color_balance.lift = (0.90, 0.95, 1.05)
        color_balance.gamma = (1.0, 1.0, 1.0)
        color_balance.gain = (1.06, 1.02, 0.94)
    except Exception:
        pass

    hue_sat = tree.nodes.new("CompositorNodeHueSat")
    hue_sat.name = "Animation_Compositor_Gentle_Saturation"
    set_node_input(hue_sat, ("Saturation", "Sat"), 1.025)
    set_node_input(hue_sat, ("Value", "Val"), 1.015)

    composite = tree.nodes.new(output_node_type)
    composite.name = "Animation_Compositor_Output"
    viewer = tree.nodes.new("CompositorNodeViewer")
    viewer.name = "Animation_Compositor_Viewer"

    render_layers.location = (-620, 120)
    glare.location = (-390, 120)
    color_balance.location = (-160, 120)
    hue_sat.location = (70, 120)
    composite.location = (320, 160)
    viewer.location = (320, -20)

    tree.links.new(render_layers.outputs["Image"], glare.inputs["Image"])
    tree.links.new(glare.outputs["Image"], color_balance.inputs["Image"])
    tree.links.new(color_balance.outputs["Image"], hue_sat.inputs["Image"])
    tree.links.new(hue_sat.outputs["Image"], composite.inputs["Image"])
    tree.links.new(hue_sat.outputs["Image"], viewer.inputs["Image"])

    return [node.name for node in (render_layers, glare, color_balance, hue_sat, composite, viewer)]


def volume_material(name: str, color: tuple[float, float, float, float], density: float) -> bpy.types.Material:
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = color
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    nodes.clear()
    output = nodes.new("ShaderNodeOutputMaterial")
    output.location = (240, 0)
    try:
        volume = nodes.new("ShaderNodeVolumePrincipled")
        set_node_input(volume, ("Color",), color)
        set_node_input(volume, ("Density",), density)
        set_node_input(volume, ("Anisotropy",), 0.18)
    except RuntimeError:
        volume = nodes.new("ShaderNodeVolumeScatter")
        set_node_input(volume, ("Color",), color)
        set_node_input(volume, ("Density",), density)
        set_node_input(volume, ("Anisotropy",), 0.12)
    volume.location = (0, 0)
    if "Volume" in output.inputs and "Volume" in volume.outputs:
        mat.node_tree.links.new(volume.outputs["Volume"], output.inputs["Volume"])
    return mat


def create_box_mesh(name: str, center: Vector, dimensions: Vector) -> bpy.types.Mesh:
    half = dimensions * 0.5
    verts = [
        (center.x - half.x, center.y - half.y, center.z - half.z),
        (center.x + half.x, center.y - half.y, center.z - half.z),
        (center.x + half.x, center.y + half.y, center.z - half.z),
        (center.x - half.x, center.y + half.y, center.z - half.z),
        (center.x - half.x, center.y - half.y, center.z + half.z),
        (center.x + half.x, center.y - half.y, center.z + half.z),
        (center.x + half.x, center.y + half.y, center.z + half.z),
        (center.x - half.x, center.y + half.y, center.z + half.z),
    ]
    faces = [
        (0, 1, 2, 3),
        (4, 7, 6, 5),
        (0, 4, 5, 1),
        (1, 5, 6, 2),
        (2, 6, 7, 3),
        (3, 7, 4, 0),
    ]
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, [], faces)
    mesh.update()
    return mesh


def create_volumetric_atmosphere(
    collection: bpy.types.Collection,
    points: list[Vector],
    density: float = 0.003,
) -> bpy.types.Object:
    min_v = Vector((min(point.x for point in points), min(point.y for point in points), min(point.z for point in points)))
    max_v = Vector((max(point.x for point in points), max(point.y for point in points), max(point.z for point in points)))
    min_v -= Vector((90.0, 90.0, 55.0))
    max_v += Vector((90.0, 90.0, 110.0))
    center = (min_v + max_v) * 0.5
    dimensions = max_v - min_v
    mesh = create_box_mesh("Animation_Atmosphere_Box_Mesh", center, dimensions)
    obj = bpy.data.objects.new("Animation_Atmosphere_Subtle_Volume", mesh)
    obj.data.materials.append(volume_material("anim_atmosphere_cool_subtle_volume", (0.30, 0.34, 0.42, 1.0), density))
    obj.display_type = "WIRE"
    obj.hide_select = True
    obj["atmosphere_density"] = density
    obj["purpose"] = "subtle volumetric depth for cinematic flythrough"
    link_to_collection(obj, collection)
    return obj


def configure_light_quality(light_obj: bpy.types.Object) -> None:
    data = light_obj.data
    for attr, value in [
        ("use_shadow", True),
        ("use_contact_shadow", True),
        ("contact_shadow_bias", 0.02),
        ("contact_shadow_distance", 35.0),
        ("shadow_soft_size", getattr(data, "shadow_soft_size", 1.0)),
    ]:
        if hasattr(data, attr):
            try:
                setattr(data, attr, value)
            except TypeError:
                pass


def create_area_light(
    name: str,
    location: Vector,
    energy: float,
    size: float,
    color: tuple[float, float, float],
    collection: bpy.types.Collection,
    target: bpy.types.Object | None = None,
    parent: bpy.types.Object | None = None,
) -> bpy.types.Object:
    light = bpy.data.lights.new(name, "AREA")
    light.energy = energy
    light.size = size
    light.color = color
    obj = bpy.data.objects.new(name, light)
    obj.location = location
    obj.parent = parent
    if hasattr(obj, "visible_camera"):
        obj.visible_camera = False
    configure_light_quality(obj)
    link_to_collection(obj, collection)
    if target is not None:
        constraint = obj.constraints.new(type="TRACK_TO")
        constraint.target = target
        constraint.track_axis = "TRACK_NEGATIVE_Z"
        constraint.up_axis = "UP_Y"
    return obj


def create_point_light(
    name: str,
    location: Vector,
    energy: float,
    color: tuple[float, float, float],
    collection: bpy.types.Collection,
    parent: bpy.types.Object | None = None,
    soft_size: float = 18.0,
) -> bpy.types.Object:
    light = bpy.data.lights.new(name, "POINT")
    light.energy = energy
    light.color = color
    light.shadow_soft_size = soft_size
    obj = bpy.data.objects.new(name, light)
    obj.parent = parent
    obj.location = location
    configure_light_quality(obj)
    link_to_collection(obj, collection)
    return obj


def create_sun_light(
    name: str,
    location: Vector,
    energy: float,
    angle: float,
    color: tuple[float, float, float],
    collection: bpy.types.Collection,
    target: bpy.types.Object,
) -> bpy.types.Object:
    light = bpy.data.lights.new(name, "SUN")
    light.energy = energy
    light.angle = angle
    light.color = color
    obj = bpy.data.objects.new(name, light)
    obj.location = location
    configure_light_quality(obj)
    link_to_collection(obj, collection)
    constraint = obj.constraints.new(type="TRACK_TO")
    constraint.target = target
    constraint.track_axis = "TRACK_NEGATIVE_Z"
    constraint.up_axis = "UP_Y"
    return obj


def configure_cinematic_lighting(collection: bpy.types.Collection, target: bpy.types.Object, scene_center: Vector) -> list[str]:
    light_target = bpy.data.objects.new("Animation_Light_Target", None)
    light_target.empty_display_type = "PLAIN_AXES"
    light_target.empty_display_size = 7.0
    light_target.location = scene_center
    link_to_collection(light_target, collection)

    lights = [
        create_area_light(
            "Animation_Soft_Key_Area",
            Vector((-66.0, -92.0, 126.0)),
            5800.0,
            106.0,
            (1.0, 0.88, 0.72),
            collection,
            light_target,
        ),
        create_sun_light(
            "Animation_Key_Sun",
            Vector((-70.0, -95.0, 155.0)),
            1.25,
            0.16,
            (1.0, 0.91, 0.78),
            collection,
            light_target,
        ),
        create_point_light(
            "Animation_Cool_Fill",
            Vector((86.0, 28.0, 92.0)),
            12500.0,
            (0.58, 0.72, 1.0),
            collection,
            None,
            95.0,
        ),
        create_area_light(
            "Animation_Long_Rim_Area",
            Vector((60.0, 74.0, 118.0)),
            5600.0,
            48.0,
            (0.62, 0.88, 1.0),
            collection,
            light_target,
        ),
        create_point_light(
            "Animation_Rim_Point",
            Vector((24.0, 52.0, 175.0)),
            30000.0,
            (0.72, 0.93, 1.0),
            collection,
            None,
            42.0,
        ),
        create_point_light(
            "Animation_Camera_Eye_Light",
            Vector((0.0, 0.0, -18.0)),
            360.0,
            (0.92, 0.96, 1.0),
            collection,
            target,
            24.0,
        ),
    ]
    return [light_target.name, *(light.name for light in lights)]


def set_material_alpha(mat: bpy.types.Material, frame: int, alpha: float, emission_strength: float | None = None) -> None:
    mat.diffuse_color = (mat.diffuse_color[0], mat.diffuse_color[1], mat.diffuse_color[2], alpha)
    mat.keyframe_insert(data_path="diffuse_color", frame=frame)
    bsdf = mat.node_tree.nodes.get("Principled BSDF") if mat.use_nodes else None
    if bsdf:
        if "Alpha" in bsdf.inputs:
            bsdf.inputs["Alpha"].default_value = alpha
            bsdf.inputs["Alpha"].keyframe_insert(data_path="default_value", frame=frame)
        if emission_strength is not None and "Emission Strength" in bsdf.inputs:
            bsdf.inputs["Emission Strength"].default_value = emission_strength
            bsdf.inputs["Emission Strength"].keyframe_insert(data_path="default_value", frame=frame)


def set_light_energy(light_obj: bpy.types.Object, frame: int, energy: float) -> None:
    light_obj.data.energy = energy
    light_obj.data.keyframe_insert(data_path="energy", frame=frame)


def animate_light_window(
    light_obj: bpy.types.Object,
    start: int,
    end: int,
    peak_energy: float,
    fade: int = 10,
) -> None:
    start = max(1, start)
    end = max(start + 1, end)
    for frame, energy in [
        (max(1, start - fade), 0.0),
        (start, 0.0),
        (min(end, start + fade), peak_energy),
        (max(start, end - fade), peak_energy),
        (end, 0.0),
        (end + fade, 0.0),
    ]:
        set_light_energy(light_obj, frame, energy)


def animate_material_window(
    mat: bpy.types.Material,
    start: int,
    end: int,
    peak_alpha: float,
    peak_emission: float = 0.7,
    fade: int = 8,
) -> None:
    start = max(1, start)
    end = max(start + 1, end)
    for frame, alpha, emission in [
        (max(1, start - fade), 0.0, 0.0),
        (start, 0.0, 0.0),
        (min(end, start + fade), peak_alpha, peak_emission),
        (max(start, end - fade), peak_alpha, peak_emission),
        (end, 0.0, 0.0),
        (end + fade, 0.0, 0.0),
    ]:
        set_material_alpha(mat, frame, alpha, emission)


def curve_points(obj_name: str) -> list[Vector]:
    obj = bpy.data.objects.get(obj_name)
    if obj is None or obj.type != "CURVE":
        raise RuntimeError(f"Expected curve object {obj_name!r}")
    points: list[Vector] = []
    for spline in obj.data.splines:
        if spline.type == "POLY":
            for point in spline.points:
                points.append(obj.matrix_world @ Vector((point.co.x, point.co.y, point.co.z)))
    if len(points) < 2:
        raise RuntimeError(f"Curve {obj_name!r} does not contain enough points")
    return points


def polyline_length(points: list[Vector]) -> float:
    return sum((points[index] - points[index - 1]).length for index in range(1, len(points)))


def point_at_length(points: list[Vector], distance: float) -> Vector:
    if distance <= 0:
        return points[0].copy()
    remaining = distance
    for index in range(1, len(points)):
        segment = points[index] - points[index - 1]
        length = segment.length
        if length >= remaining:
            return points[index - 1].lerp(points[index], remaining / max(length, 1e-9))
        remaining -= length
    return points[-1].copy()


def point_at_fraction(points: list[Vector], fraction: float) -> Vector:
    return point_at_length(points, polyline_length(points) * max(0.0, min(1.0, fraction)))


def create_curve_object(
    name: str,
    points: list[Vector],
    bevel_depth: float,
    mat: bpy.types.Material,
    collection: bpy.types.Collection,
) -> bpy.types.Object:
    curve = bpy.data.curves.new(name, "CURVE")
    curve.dimensions = "3D"
    curve.resolution_u = 2
    curve.bevel_depth = bevel_depth
    curve.bevel_resolution = 4
    spline = curve.splines.new("POLY")
    spline.points.add(len(points) - 1)
    for spline_point, point in zip(spline.points, points):
        spline_point.co = (point.x, point.y, point.z, 1.0)
    curve.materials.append(mat)
    obj = bpy.data.objects.new(name, curve)
    link_to_collection(obj, collection)
    return obj


def circle_points(center: Vector, radius: float, plane: str, count: int = 96) -> list[Vector]:
    points = []
    for index in range(count + 1):
        theta = 2.0 * math.pi * index / count
        c = math.cos(theta) * radius
        s = math.sin(theta) * radius
        if plane == "xy":
            points.append(center + Vector((c, s, 0.0)))
        elif plane == "xz":
            points.append(center + Vector((c, 0.0, s)))
        else:
            points.append(center + Vector((0.0, c, s)))
    return points


def create_target_rings(
    name: str,
    center: Vector,
    radius: float,
    mat: bpy.types.Material,
    collection: bpy.types.Collection,
) -> list[bpy.types.Object]:
    rings = []
    for plane in ("xy", "xz", "yz"):
        rings.append(create_curve_object(f"{name}_{plane}", circle_points(center, radius, plane), 0.035, mat, collection))
    return rings


def mesh_objects_for_asset(asset_name: str) -> list[bpy.types.Object]:
    return [
        obj
        for obj in bpy.data.objects
        if obj.type == "MESH" and obj.name.startswith(asset_name)
    ]


def create_mesh_highlight_overlays(
    asset_name: str,
    name_prefix: str,
    mat: bpy.types.Material,
    collection: bpy.types.Collection,
    start: int,
    end: int,
    scale: float = 1.012,
) -> list[bpy.types.Object]:
    overlays = []
    for source in mesh_objects_for_asset(asset_name):
        overlay = source.copy()
        overlay.data = source.data.copy()
        overlay.animation_data_clear()
        overlay.name = f"{name_prefix}_{source.name}"
        overlay.scale = source.scale * scale
        overlay.hide_viewport = True
        overlay.hide_render = True
        overlay.data.materials.clear()
        overlay.data.materials.append(mat)
        link_to_collection(overlay, collection)
        animate_visibility_window(overlay, start, end)
        overlays.append(overlay)
    return overlays


def set_object_visibility(obj: bpy.types.Object, frame: int, visible: bool) -> None:
    obj.hide_viewport = not visible
    obj.hide_render = not visible
    obj.keyframe_insert(data_path="hide_viewport", frame=frame)
    obj.keyframe_insert(data_path="hide_render", frame=frame)


def animate_visibility_window(obj: bpy.types.Object, start: int, end: int) -> None:
    for frame, visible in [(1, False), (max(1, start - 1), False), (start, True), (end, True), (end + 1, False)]:
        set_object_visibility(obj, frame, visible)


def hide_existing_text_labels() -> list[str]:
    hidden = []
    for obj in bpy.data.objects:
        if obj.type != "FONT" or obj.name.startswith("Animation_"):
            continue
        obj.hide_viewport = True
        obj.hide_render = True
        hidden.append(obj.name)
    return hidden


def hide_existing_backdrops() -> list[str]:
    hidden = []
    for obj in bpy.data.objects:
        name = obj.name.lower()
        if "backdrop" not in name and "background" not in name:
            continue
        obj.hide_viewport = True
        obj.hide_render = True
        hidden.append(obj.name)
    return hidden


def create_label(
    name: str,
    text: str,
    location: Vector,
    size: float,
    mat: bpy.types.Material,
    collection: bpy.types.Collection,
    camera: bpy.types.Object,
    start: int,
    end: int,
) -> bpy.types.Object:
    curve = bpy.data.curves.new(name, "FONT")
    curve.body = text
    curve.align_x = "CENTER"
    curve.align_y = "CENTER"
    curve.size = size
    curve.materials.append(mat)
    obj = bpy.data.objects.new(name, curve)
    obj.location = location
    link_to_collection(obj, collection)
    constraint = obj.constraints.new(type="COPY_ROTATION")
    constraint.target = camera
    animate_visibility_window(obj, start, end)
    return obj


def create_camera_label(
    name: str,
    text: str,
    camera: bpy.types.Object,
    local_location: tuple[float, float, float],
    size: float,
    mat: bpy.types.Material,
    collection: bpy.types.Collection,
    start: int,
    end: int,
) -> bpy.types.Object:
    curve = bpy.data.curves.new(name, "FONT")
    curve.body = text
    curve.align_x = "CENTER"
    curve.align_y = "CENTER"
    curve.size = size
    curve.materials.append(mat)
    obj = bpy.data.objects.new(name, curve)
    obj.parent = camera
    obj.location = local_location
    link_to_collection(obj, collection)
    animate_visibility_window(obj, start, end)
    return obj


def create_camera_panel(
    name: str,
    camera: bpy.types.Object,
    local_location: tuple[float, float, float],
    width: float,
    height: float,
    mat: bpy.types.Material,
    collection: bpy.types.Collection,
    start: int,
    end: int,
) -> bpy.types.Object:
    mesh = bpy.data.meshes.new(name)
    half_w = width * 0.5
    half_h = height * 0.5
    mesh.from_pydata(
        [(-half_w, -half_h, 0.0), (half_w, -half_h, 0.0), (half_w, half_h, 0.0), (-half_w, half_h, 0.0)],
        [],
        [(0, 1, 2, 3)],
    )
    mesh.materials.append(mat)
    obj = bpy.data.objects.new(name, mesh)
    obj.parent = camera
    obj.location = local_location
    link_to_collection(obj, collection)
    animate_visibility_window(obj, start, end)
    return obj


def vector_from_report(report: dict, name: str) -> Vector:
    for asset in report.get("pdb_assets", []):
        if asset.get("name") == name:
            return Vector(asset["location_mm"])
    raise KeyError(name)


def bbox_radius(report: dict, name: str, default: float = 3.0) -> float:
    for asset in report.get("pdb_assets", []):
        if asset.get("name") == name:
            bbox = asset.get("bbox_mm") or []
            if bbox:
                return max(max(float(v) for v in bbox) * 0.65, default)
    return default


def bbox_center(components: list[dict]) -> Vector:
    mins = []
    maxs = []
    for component in components:
        if "min_mm" in component and "max_mm" in component:
            mins.append(Vector(component["min_mm"]))
            maxs.append(Vector(component["max_mm"]))
    min_v = Vector((min(v.x for v in mins), min(v.y for v in mins), min(v.z for v in mins)))
    max_v = Vector((max(v.x for v in maxs), max(v.y for v in maxs), max(v.z for v in maxs)))
    return (min_v + max_v) * 0.5


def frame_at(seconds: float, duration_seconds: float, fps: int) -> int:
    scaled = seconds * duration_seconds / 60.0
    return max(1, int(round(1 + scaled * fps)))


def configure_render(args: argparse.Namespace, frame_end: int, frames_dir: Path) -> None:
    scene = bpy.context.scene
    scene.frame_start = 1
    scene.frame_end = frame_end
    scene.frame_set(1)
    scene.render.fps = args.fps
    scene.render.resolution_x = args.resolution_x
    scene.render.resolution_y = args.resolution_y
    scene.render.film_transparent = False
    try:
        scene.render.engine = "BLENDER_EEVEE_NEXT"
    except TypeError:
        scene.render.engine = "BLENDER_EEVEE"
    if hasattr(scene, "eevee"):
        scene.eevee.taa_render_samples = 20 if args.smoke_test else 96
        for attr, value in [
            ("use_gtao", True),
            ("gtao_distance", 24.0),
            ("gtao_factor", 1.55),
            ("use_bloom", True),
            ("bloom_intensity", 0.032),
            ("bloom_radius", 6.0),
            ("use_motion_blur", True),
            ("motion_blur_shutter", 0.18),
            ("use_raytracing", True),
            ("use_volumetric_lights", True),
            ("use_volumetric_shadows", True),
            ("volumetric_samples", 32 if args.smoke_test else 80),
            ("volumetric_sample_distribution", 0.62),
            ("volumetric_start", 0.1),
            ("volumetric_end", 700.0),
            ("volumetric_tile_size", "4"),
        ]:
            if hasattr(scene.eevee, attr):
                try:
                    setattr(scene.eevee, attr, value)
                except TypeError:
                    pass
    scene.render.image_settings.file_format = "PNG"
    scene.render.image_settings.color_mode = "RGB"
    scene.render.filepath = str(frames_dir / "frame_")


def create_animation(args: argparse.Namespace, report: dict, output_blend: Path, output_mp4: Path, frames_dir: Path) -> dict:
    dna_points = curve_points("V5_DNA_source_path")
    mrna_points = curve_points("V5_mRNA_source_path")
    dna_length = polyline_length(dna_points)
    mrna_length = polyline_length(mrna_points)

    collections = {
        "camera": scene_collection("Animation Camera"),
        "labels": scene_collection("Animation Labels"),
        "highlights": scene_collection("Animation Highlights"),
        "lights": scene_collection("Animation Lights"),
        "caption_panels": scene_collection("Animation Caption Panels"),
        "atmosphere": scene_collection("Animation Atmosphere"),
    }
    mats = {
        "dna_highlight": material("anim_highlight_dna_blue", (0.05, 0.42, 1.0, 0.0), 0.0),
        "rna_highlight": material("anim_highlight_rna_orange", (1.0, 0.44, 0.08, 0.0), 0.0),
        "protein_highlight": material("anim_highlight_protein_cyan", (0.05, 0.86, 1.0, 0.0), 0.0),
        "actin_highlight": material("anim_highlight_actin_red", (1.0, 0.15, 0.08, 0.0), 0.0),
        "compact_highlight": material("anim_highlight_compact_gold", (1.0, 0.78, 0.05, 0.0), 0.0),
        "tf_pulse": material("anim_pulse_tf_green", (0.16, 1.0, 0.38, 0.0), 0.0),
        "pol_pulse": material("anim_pulse_pol_cyan", (0.10, 0.88, 1.0, 0.0), 0.0),
        "nucleosome_pulse": material("anim_pulse_nucleosome_violet", (0.74, 0.45, 1.0, 0.0), 0.0),
        "mrna_pulse": material("anim_pulse_mrna_orange", (1.0, 0.52, 0.08, 0.0), 0.0),
        "ribosome_pulse": material("anim_pulse_ribosome_gold", (1.0, 0.76, 0.20, 0.0), 0.0),
        "trna_pulse": material("anim_pulse_trna_green", (0.38, 1.0, 0.68, 0.0), 0.0),
        "actin_pulse": material("anim_pulse_actin_red", (1.0, 0.20, 0.10, 0.0), 0.0),
        "label": material("anim_label_soft_white", (0.94, 0.965, 1.0, 1.0), 0.34),
        "caption_panel": material("anim_caption_panel_smoke", (0.015, 0.019, 0.026, 0.42), 0.0),
    }
    hidden_original_labels = hide_existing_text_labels()
    hidden_original_backdrops = hide_existing_backdrops()
    configure_view_settings()
    compositor_nodes = configure_compositor()

    duration = max(4.0, args.duration_seconds)
    fps = max(1, args.fps)
    frame_end = frame_at(60.0, duration, fps)

    camera_data = bpy.data.cameras.new("Animation_Flythrough_Camera")
    camera_data.type = "PERSP"
    camera_data.lens = 22.0
    camera_data.sensor_width = 32.0
    camera_data.clip_end = 900.0
    camera_data.dof.use_dof = True
    camera_data.dof.aperture_fstop = 9.5
    camera_data.dof.aperture_blades = 7
    camera = bpy.data.objects.new("Animation_Flythrough_Camera", camera_data)
    link_to_collection(camera, collections["camera"])
    target = bpy.data.objects.new("Animation_Flythrough_Target", None)
    target.empty_display_type = "SPHERE"
    target.empty_display_size = 2.5
    link_to_collection(target, collections["camera"])
    constraint = camera.constraints.new(type="TRACK_TO")
    constraint.target = target
    constraint.track_axis = "TRACK_NEGATIVE_Z"
    constraint.up_axis = "UP_Y"
    bpy.context.scene.camera = camera
    camera.data.dof.focus_object = target

    dna_center = bbox_center(report["dna"]["components"])
    compact_center = bbox_center(report["compact_mrna"]["components"])
    actin = vector_from_report(report, "Actin protein")
    ribosome = (vector_from_report(report, "Ribosome small subunit") + vector_from_report(report, "Ribosome large subunit")) * 0.5
    pol = vector_from_report(report, "RNA polymerase II elongation complex")
    nucleosome = vector_from_report(report, "Nucleosome")
    tf_focus = vector_from_report(report, "p53 tetramer bound to DNA")
    mrna_mid = point_at_fraction(mrna_points, 0.52)
    rbp_focus = (
        vector_from_report(report, "Pumilio RBP")
        + vector_from_report(report, "MS2 coat protein MCP")
        + vector_from_report(report, "Argonaute")
        + vector_from_report(report, "Poly(A)-binding RBP")
        + vector_from_report(report, "HuR-like RBP")
    ) / 5.0
    scene_center = (dna_center + compact_center + actin + ribosome + mrna_mid) / 5.0
    lighting_objects = configure_cinematic_lighting(collections["lights"], camera, scene_center)
    atmosphere = create_volumetric_atmosphere(
        collections["atmosphere"],
        [*dna_points, *mrna_points, dna_center, compact_center, actin, ribosome, pol, nucleosome, tf_focus, mrna_mid, rbp_focus],
        density=0.003,
    )

    top_down = (scene_center + Vector((-18.0, -44.0, 205.0)), scene_center + Vector((0.0, -2.0, -5.0)), 20.0, 11.0)
    overview = (dna_center + Vector((-52.0, -118.0, 126.0)), scene_center + Vector((2.0, -7.0, 8.0)), 24.0, 8.0)
    dna_entry = (point_at_fraction(dna_points, 0.06) + Vector((-31.0, -43.0, 34.0)), point_at_fraction(dna_points, 0.06), 34.0, 5.6)
    tf_shot = (tf_focus + Vector((-11.5, -17.0, 9.8)), tf_focus, 84.0, 1.45)
    pol_shot = (pol + Vector((-16.0, -21.0, 13.0)), pol, 84.0, 1.55)
    nuc_rna = (nucleosome + Vector((-15.0, -22.0, 15.5)), (nucleosome + point_at_fraction(mrna_points, 0.12)) * 0.5, 42.0, 3.2)
    mrna_wide = (mrna_mid + Vector((-78.0, -112.0, 88.0)), mrna_mid, 22.0, 8.5)
    rbp_ribosome = ((rbp_focus + ribosome) * 0.5 + Vector((-24.0, -48.0, 38.0)), (rbp_focus + ribosome) * 0.5, 36.0, 3.5)
    ribosome_shot = (ribosome + Vector((9.0, -30.0, 31.0)), ribosome, 82.0, 1.55)
    actin_rise = ((ribosome + actin) * 0.5 + Vector((-11.0, -32.0, 46.0)), (ribosome + actin) * 0.5, 34.0, 3.0)
    actin_shot = (actin + Vector((-9.5, -16.0, 11.0)), actin, 90.0, 1.35)

    storyboard = [
        (0.0, "top_down_scale", *top_down),
        (6.0, "top_down_hold", *top_down),
        (10.0, "overview_oblique", *overview),
        (14.0, "dna_entry", *dna_entry),
        (18.0, "tf_binding_closeup", *tf_shot),
        (23.0, "tf_binding_hold", *tf_shot),
        (27.0, "polymerase_closeup", *pol_shot),
        (32.0, "polymerase_hold", *pol_shot),
        (36.0, "nucleosome_to_rna", *nuc_rna),
        (41.0, "full_mrna_zoomout", *mrna_wide),
        (46.0, "full_mrna_hold", *mrna_wide),
        (48.5, "rbp_to_ribosome", *rbp_ribosome),
        (50.5, "ribosome_closeup", *ribosome_shot),
        (52.5, "ribosome_hold", *ribosome_shot),
        (53.5, "rise_to_actin", *actin_rise),
        (55.0, "actin_closeup", *actin_shot),
        (60.0, "actin_close_hold", *actin_shot),
    ]

    for seconds, _name, location, focus, lens, fstop in storyboard:
        frame = frame_at(seconds, duration, fps)
        camera.location = location
        camera.data.lens = lens
        camera.data.dof.aperture_fstop = fstop
        target.location = focus
        camera.keyframe_insert(data_path="location", frame=frame)
        camera.data.keyframe_insert(data_path="lens", frame=frame)
        camera.data.dof.keyframe_insert(data_path="aperture_fstop", frame=frame)
        target.keyframe_insert(data_path="location", frame=frame)

    dna_highlight = create_curve_object("Animation_highlight_DNA_3954bp_path", dna_points, 0.18, mats["dna_highlight"], collections["highlights"])
    mrna_highlight = create_curve_object("Animation_highlight_mRNA_1852nt_path", mrna_points, 0.13, mats["rna_highlight"], collections["highlights"])
    dna_highlight["educational_callout"] = "ACTB promoter + gene DNA; 3954 bp; 537.7 mm"
    mrna_highlight["educational_callout"] = "actin mRNA; 1852 nt; 222.2 mm"
    animate_material_window(mats["dna_highlight"], frame_at(5, duration, fps), frame_at(34, duration, fps), 0.72, 1.0)
    animate_material_window(mats["rna_highlight"], frame_at(31, duration, fps), frame_at(58, duration, fps), 0.76, 1.0)

    highlight_objects = [dna_highlight.name, mrna_highlight.name]

    pulse_specs = [
        ("tf", "p53 tetramer bound to DNA", mats["tf_pulse"], 17, 24, tf_focus + Vector((-3.0, -6.0, 7.0)), 1150.0, 0.30, 1.8),
        ("pol", "RNA polymerase II elongation complex", mats["pol_pulse"], 26, 33, pol + Vector((-5.0, -8.0, 6.0)), 1850.0, 0.22, 1.35),
        ("nucleosome", "Nucleosome", mats["nucleosome_pulse"], 34, 38, nucleosome + Vector((-4.0, -6.0, 5.0)), 1100.0, 0.20, 1.25),
        ("ribosome", "Ribosome large subunit", mats["ribosome_pulse"], 49, 54, ribosome + Vector((4.0, -10.0, 9.0)), 950.0, 0.045, 0.22),
        ("trna", "Standalone tRNA", mats["trna_pulse"], 49, 54, vector_from_report(report, "Standalone tRNA") + Vector((1.0, -3.0, 4.0)), 420.0, 0.16, 0.7),
        ("actin", "Actin protein", mats["actin_pulse"], 54, 60, actin + Vector((-2.0, -5.0, 4.0)), 1900.0, 0.24, 1.35),
    ]
    pulse_overlay_objects = []
    pulse_lights = []
    for name, asset_name, mat, start, end, light_location, light_energy, peak_alpha, peak_emission in pulse_specs:
        start_frame = frame_at(start, duration, fps)
        end_frame = frame_at(end, duration, fps)
        animate_material_window(mat, start_frame, end_frame, peak_alpha, peak_emission, fade=6)
        for overlay in create_mesh_highlight_overlays(
            asset_name,
            f"Animation_pulse_{name}",
            mat,
            collections["highlights"],
            start_frame,
            end_frame,
        ):
            overlay["highlight_target"] = asset_name
            pulse_overlay_objects.append(overlay.name)
        beat_light = create_point_light(
            f"Animation_Beat_Light_{name}",
            light_location,
            0.0,
            tuple(mat.diffuse_color[:3]),
            collections["lights"],
            None,
            16.0,
        )
        animate_light_window(beat_light, start_frame, end_frame, light_energy, fade=8)
        pulse_lights.append(beat_light.name)

    animate_material_window(mats["mrna_pulse"], frame_at(39, duration, fps), frame_at(47, duration, fps), 0.32, 1.35, fade=8)
    mrna_pulse = create_curve_object("Animation_pulse_full_mRNA_color_hold", mrna_points, 0.22, mats["mrna_pulse"], collections["highlights"])
    animate_visibility_window(mrna_pulse, frame_at(39, duration, fps), frame_at(47, duration, fps))
    pulse_overlay_objects.append(mrna_pulse.name)
    highlight_objects.extend(pulse_overlay_objects)
    lighting_objects.extend(pulse_lights)

    label_specs = [
        ("label_dna", "ACTB DNA\n3,954 bp / 537.7 mm", point_at_fraction(dna_points, 0.05) + Vector((-4.0, 5.0, 2.5)), 1.05, 8, 18),
        ("label_tf", "p53 tetramer\nbound to DNA", tf_focus + Vector((2.6, 3.4, 2.1)), 0.48, 16, 25),
        ("label_pol", "Pol II", pol + Vector((-6.0, -6.0, 4.2)), 0.30, 25, 35),
        ("label_nucleosome", "nucleosome loop\nwrapped DNA", nucleosome + Vector((3.6, 3.6, 2.8)), 0.5, 33, 40),
        ("label_mrna", "full actin mRNA\n1,852 nt / 222.2 mm", point_at_fraction(mrna_points, 0.48) + Vector((8.0, 6.0, 7.0)), 0.9, 38, 48),
        ("label_rbps", "RNA-binding proteins\npost-transcriptional control", rbp_focus + Vector((8.0, 7.0, 5.0)), 0.78, 43, 51),
        ("label_ribosome", "ribosome + tRNA\ntranslation", ribosome + Vector((6.0, 6.0, 5.0)), 0.52, 48, 55),
        ("label_compact", "compact mRNP-like RNA\n10.5 x 8.0 x 3.4 mm envelope", compact_center + Vector((0.0, 7.0, 4.0)), 0.48, 53, 55),
        ("label_actin", "actin protein\n2.7 x 1.8 x 2.1 mm", actin + Vector((3.4, 3.8, 2.4)), 0.36, 54, 60),
    ]
    labels = []
    for name, text, location, size, start, end in label_specs:
        labels.append(
            create_label(
                f"Animation_{name}",
                text,
                location,
                size,
                mats["label"],
                collections["labels"],
                camera,
                frame_at(start, duration, fps),
                frame_at(end, duration, fps),
            ).name
        )
    caption_specs = [
        ("overview", "overall scale: DNA, mRNA, ribosome, actin", 0, 10),
        ("dna", "ACTB DNA: 3,954 bp -> 537.7 mm", 8, 18),
        ("tf", "transcription factor binds DNA", 16, 25),
        ("pol", "RNA Pol II starts mRNA", 25, 35),
        ("mrna", "full mRNA: 1,852 nt -> 222.2 mm", 38, 48),
        ("ribosome", "translation: ribosome + tRNA", 48, 54),
        ("actin", "actin endpoint: 2.7 x 1.8 x 2.1 mm", 54, 60),
    ]
    caption_panels = []
    for name, text, start, end in caption_specs:
        start_frame = frame_at(start, duration, fps)
        end_frame = frame_at(end, duration, fps)
        caption_panels.append(
            create_camera_panel(
                f"Animation_caption_panel_{name}",
                camera,
                (0.0, -1.36, -18.0),
                7.6,
                0.95,
                mats["caption_panel"],
                collections["caption_panels"],
                start_frame,
                end_frame,
            ).name
        )
        labels.append(
            create_camera_label(
                f"Animation_caption_{name}",
                text,
                camera,
                (0.0, -1.40, -17.8),
                0.21,
                mats["label"],
                collections["labels"],
                start_frame,
                end_frame,
            ).name
        )

    configure_render(args, frame_end, frames_dir)
    for action in bpy.data.actions:
        fcurves = getattr(action, "fcurves", None)
        if fcurves is None:
            continue
        for fcurve in fcurves:
            for keyframe in fcurve.keyframe_points:
                keyframe.interpolation = "BEZIER"

    output_blend.parent.mkdir(parents=True, exist_ok=True)
    bpy.ops.wm.save_as_mainfile(filepath=str(output_blend))

    rendered_frames = 0
    if not args.skip_video_render:
        frames_dir.mkdir(parents=True, exist_ok=True)
        for old_frame in frames_dir.glob("frame_*.png"):
            old_frame.unlink()
        bpy.ops.render.render(animation=True)
        rendered_frames = len(list(frames_dir.glob("frame_*.png")))

    return {
        "source_blend": str(args.source_blend),
        "source_report": str(args.source_report),
        "experiment_blend": str(output_blend),
        "mp4": str(output_mp4),
        "frames_dir": str(frames_dir),
        "rendered": rendered_frames > 0,
        "rendered_frames": rendered_frames,
        "duration_seconds": duration,
        "fps": fps,
        "frame_start": 1,
        "frame_end": frame_end,
        "resolution": [args.resolution_x, args.resolution_y],
        "smoke_test": bool(args.smoke_test),
        "camera": camera.name,
        "camera_type": camera.data.type,
        "target": target.name,
        "source_path_lengths_mm": {"dna": dna_length, "mrna": mrna_length},
        "storyboard": [
            {
                "story_seconds": seconds,
                "timeline_seconds": round(seconds * duration / 60.0, 3),
                "frame": frame_at(seconds, duration, fps),
                "name": name,
                "camera_location_mm": [location.x, location.y, location.z],
                "focus_mm": [focus.x, focus.y, focus.z],
                "lens_mm": lens,
                "dof_aperture_fstop": fstop,
            }
            for seconds, name, location, focus, lens, fstop in storyboard
        ],
        "labels": labels,
        "caption_panels": caption_panels,
        "lighting_objects": lighting_objects,
        "atmosphere_object": atmosphere.name,
        "atmosphere_density": atmosphere["atmosphere_density"],
        "compositor_nodes": compositor_nodes,
        "pulse_overlay_objects": pulse_overlay_objects,
        "pulse_lights": pulse_lights,
        "hidden_original_labels": hidden_original_labels,
        "hidden_original_backdrops": hidden_original_backdrops,
        "highlight_objects": highlight_objects,
    }


def main() -> None:
    args = parse_args()
    args.source_blend = as_path(args.source_blend)
    args.source_report = as_path(args.source_report)
    args.output_dir = as_path(args.output_dir)
    if args.smoke_test:
        args.duration_seconds = min(args.duration_seconds, 8.0)
        args.fps = min(args.fps, 6)
        args.resolution_x = min(args.resolution_x, 640)
        args.resolution_y = min(args.resolution_y, 360)

    output_blend = as_path(args.output_blend) if args.output_blend else args.output_dir / "v5_flythrough_animation.blend"
    output_mp4 = as_path(args.output_mp4) if args.output_mp4 else args.output_dir / ("v5_flythrough_animation_smoke.mp4" if args.smoke_test else "v5_flythrough_animation_1080p.mp4")
    report_path = args.output_dir / "v5_flythrough_animation_report.json"
    frames_dir = args.output_dir / ("frames_smoke" if args.smoke_test else "frames")

    if not args.source_blend.exists():
        raise FileNotFoundError(f"Canonical V5 blend not found: {args.source_blend}")
    if not args.source_report.exists():
        raise FileNotFoundError(f"Canonical V5 report not found: {args.source_report}")

    bpy.ops.wm.open_mainfile(filepath=str(args.source_blend))
    report = load_json(args.source_report)
    animation_report = create_animation(args, report, output_blend, output_mp4, frames_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(animation_report, indent=2), encoding="utf-8")
    print(f"Wrote {output_blend}")
    print(f"Wrote {report_path}")
    if animation_report["rendered"]:
        print(f"Rendered {animation_report['rendered_frames']} frames for {output_mp4}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        import traceback

        traceback.print_exc()
        sys.exit(1)
