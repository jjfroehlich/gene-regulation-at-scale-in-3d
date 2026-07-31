#!/usr/bin/env python3
"""Build and optionally render the educational flythrough animation."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import bpy
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Vector


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE_BLEND = ROOT / "outputs" / "canonical" / "gene_expression_surface_style.blend"
DEFAULT_REPORT = ROOT / "outputs" / "canonical" / "gene_expression_surface_scene_report.json"
DEFAULT_OUTPUT_DIR = ROOT / "experiments" / "flythrough_animation" / "outputs"
STORY_DURATION_SECONDS = 66.0


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
    parser.add_argument("--frames-dir", type=Path)
    parser.add_argument("--render-profile", choices=("final", "review", "smoke"), default="final")
    parser.add_argument("--duration-seconds", type=float, default=STORY_DURATION_SECONDS)
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
            1800.0,
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


def hide_existing_overview_guides() -> list[str]:
    hidden = []
    for obj in bpy.data.objects:
        if not obj.name.startswith("overview_group_"):
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


def create_label_leader(
    name: str,
    anchor: Vector,
    label_location: Vector,
    mat: bpy.types.Material,
    collection: bpy.types.Collection,
    start: int,
    end: int,
) -> bpy.types.Object:
    direction = label_location - anchor
    endpoint = label_location - direction.normalized() * min(0.75, direction.length * 0.18)
    obj = create_curve_object(name, [anchor, endpoint], 0.014, mat, collection)
    animate_visibility_window(obj, start, end)
    return obj


def camera_plane_label_location(
    anchor: Vector,
    camera_location: Vector,
    focus: Vector,
    screen_x: float,
    screen_y: float,
) -> Vector:
    """Offset a molecule-attached label in the intended shot's screen plane."""
    direction = (focus - camera_location).normalized()
    right = direction.cross(Vector((0.0, 0.0, 1.0)))
    if right.length < 1e-6:
        right = Vector((1.0, 0.0, 0.0))
    else:
        right.normalize()
    up = right.cross(direction).normalized()
    return anchor + right * screen_x + up * screen_y


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
    scaled = seconds * duration_seconds / STORY_DURATION_SECONDS
    return max(1, int(round(1 + scaled * fps)))


def create_bezier_camera_path(
    name: str,
    points: list[Vector],
    collection: bpy.types.Collection,
) -> bpy.types.Object:
    curve = bpy.data.curves.new(name, "CURVE")
    curve.dimensions = "3D"
    curve.resolution_u = 18
    curve.render_resolution_u = 24
    spline = curve.splines.new("BEZIER")
    spline.bezier_points.add(len(points) - 1)
    for bezier_point, point in zip(spline.bezier_points, points):
        bezier_point.co = point
        try:
            bezier_point.handle_left_type = "AUTO_CLAMPED"
            bezier_point.handle_right_type = "AUTO_CLAMPED"
        except TypeError:
            bezier_point.handle_left_type = "AUTO"
            bezier_point.handle_right_type = "AUTO"
    obj = bpy.data.objects.new(name, curve)
    link_to_collection(obj, collection)
    obj.hide_viewport = True
    obj.hide_render = True
    obj["continuous_camera_path"] = True
    return obj


def bezier_path_control_progress(path: bpy.types.Object, samples_per_segment: int = 48) -> list[float]:
    """Return arc-length-normalized Follow Path factors for each Bézier control point."""
    bpy.context.view_layer.update()
    spline = path.data.splines[0]
    points = spline.bezier_points
    if len(points) < 2:
        return [0.0]
    cumulative = [0.0]
    total = 0.0
    for index in range(len(points) - 1):
        current = points[index]
        following = points[index + 1]
        p0 = current.co.copy()
        p1 = current.handle_right.copy()
        p2 = following.handle_left.copy()
        p3 = following.co.copy()
        previous = p0
        segment_length = 0.0
        for sample in range(1, samples_per_segment + 1):
            t = sample / samples_per_segment
            omt = 1.0 - t
            point = omt**3 * p0 + 3.0 * omt**2 * t * p1 + 3.0 * omt * t**2 * p2 + t**3 * p3
            segment_length += (point - previous).length
            previous = point
        total += segment_length
        cumulative.append(total)
    if total <= 1e-9:
        return [index / (len(points) - 1) for index in range(len(points))]
    return [distance / total for distance in cumulative]


def create_camera_bracket(
    name: str,
    camera: bpy.types.Object,
    x: float,
    y_bottom: float,
    y_top: float,
    z: float,
    mat: bpy.types.Material,
    collection: bpy.types.Collection,
    start: int,
    end: int,
) -> bpy.types.Object:
    tick = 0.42
    points = [
        Vector((x - tick, y_top, z)),
        Vector((x, y_top, z)),
        Vector((x, y_bottom, z)),
        Vector((x - tick, y_bottom, z)),
    ]
    obj = create_curve_object(name, points, 0.018, mat, collection)
    obj.parent = camera
    animate_visibility_window(obj, start, end)
    return obj


def camera_plane_projection(camera: bpy.types.Object, world_point: Vector, plane_z: float) -> Vector:
    local = camera.matrix_world.inverted() @ world_point
    if abs(local.z) < 1e-8:
        return Vector((0.0, 0.0, plane_z))
    scale = plane_z / local.z
    return Vector((local.x * scale, local.y * scale, plane_z))


def create_animated_camera_pointer(
    name: str,
    camera: bpy.types.Object,
    world_target: Vector,
    label_endpoint: Vector,
    plane_z: float,
    mat: bpy.types.Material,
    collection: bpy.types.Collection,
    start: int,
    end: int,
) -> bpy.types.Object:
    scene = bpy.context.scene
    scene.frame_set(start)
    target_endpoint = camera_plane_projection(camera, world_target, plane_z)
    obj = create_curve_object(name, [label_endpoint, target_endpoint], 0.018, mat, collection)
    obj.parent = camera
    animate_visibility_window(obj, start, end)
    target_point = obj.data.splines[0].points[1]
    for frame in range(start, end + 1):
        scene.frame_set(frame)
        bpy.context.view_layer.update()
        projected = camera_plane_projection(camera, world_target, plane_z)
        target_point.co = (*projected, 1.0)
        target_point.keyframe_insert(data_path="co", frame=frame)
    scene.frame_set(start)
    return obj


def camera_pointer_projection_validation(
    pointer: bpy.types.Object,
    camera: bpy.types.Object,
    world_target: Vector,
    start: int,
    end: int,
) -> dict:
    scene = bpy.context.scene
    point = pointer.data.splines[0].points[1]
    max_error = 0.0
    max_error_frame = start
    for frame in range(start, end + 1):
        scene.frame_set(frame)
        bpy.context.view_layer.update()
        actual_world = pointer.matrix_world @ Vector(point.co[:3])
        actual = world_to_camera_view(scene, camera, actual_world)
        expected = world_to_camera_view(scene, camera, world_target)
        error = math.hypot(actual.x - expected.x, actual.y - expected.y)
        if error > max_error:
            max_error = error
            max_error_frame = frame
    scene.frame_set(start)
    return {
        "frames_sampled": end - start + 1,
        "maximum_normalized_viewport_error": max_error,
        "maximum_error_frame": max_error_frame,
        "threshold": 0.01,
        "passed": max_error <= 0.01,
    }


def set_material_emission(mat: bpy.types.Material, frame: int, strength: float) -> bool:
    bsdf = mat.node_tree.nodes.get("Principled BSDF") if mat.use_nodes and mat.node_tree else None
    if bsdf is None or "Emission Strength" not in bsdf.inputs:
        return False
    emission = bsdf.inputs["Emission Strength"]
    emission.default_value = strength
    emission.keyframe_insert(data_path="default_value", frame=frame)
    return True


def animate_emission_window(mat: bpy.types.Material, start: int, end: int, peak: float, fade: int = 8) -> bool:
    start = max(1, start)
    end = max(start + 2, end)
    midpoint = start + (end - start) // 2
    keyed = False
    for frame, strength in [
        (max(1, start - fade), 0.0),
        (start, 0.0),
        (midpoint, peak),
        (end, 0.0),
        (end + fade, 0.0),
    ]:
        keyed = set_material_emission(mat, frame, strength) or keyed
    return keyed


def pulse_asset_emission(
    asset_name: str,
    component: str,
    start: int,
    end: int,
    peak: float,
) -> dict:
    object_names = []
    material_names = []
    for obj in mesh_objects_for_asset(asset_name):
        if obj.get("component") != component:
            continue
        object_names.append(obj.name)
        for index, source in enumerate(list(obj.data.materials)):
            if source is None:
                continue
            clone = source.copy()
            clone.name = f"Animation_emission_{asset_name}_{component}_{index}"
            clone.animation_data_clear()
            bsdf = clone.node_tree.nodes.get("Principled BSDF") if clone.use_nodes and clone.node_tree else None
            if bsdf is not None:
                base_color = bsdf.inputs.get("Base Color")
                emission_color = bsdf.inputs.get("Emission Color") or bsdf.inputs.get("Emission")
                if base_color is not None and emission_color is not None:
                    emission_color.default_value = base_color.default_value
            obj.data.materials[index] = clone
            if animate_emission_window(clone, start, end, peak):
                material_names.append(clone.name)
    if not object_names or not material_names:
        raise RuntimeError(f"Could not create source-mesh emission pulse for {asset_name!r} component {component!r}")
    return {
        "asset": asset_name,
        "component": component,
        "objects": object_names,
        "materials": material_names,
        "start_frame": start,
        "end_frame": end,
        "peak_emission": peak,
    }


def set_animation_interpolation() -> None:
    for action in bpy.data.actions:
        fcurves = getattr(action, "fcurves", None)
        if fcurves is None:
            continue
        for fcurve in fcurves:
            discrete = "hide_viewport" in fcurve.data_path or "hide_render" in fcurve.data_path
            for keyframe in fcurve.keyframe_points:
                if discrete:
                    keyframe.interpolation = "CONSTANT"
                    continue
                keyframe.interpolation = "BEZIER"
                keyframe.handle_left_type = "AUTO_CLAMPED"
                keyframe.handle_right_type = "AUTO_CLAMPED"


def camera_motion_continuity(
    camera: bpy.types.Object,
    target: bpy.types.Object,
    frame_start: int,
    frame_end: int,
    duration_seconds: float,
    fps: int,
    moving_end_frame: int | None = None,
) -> dict:
    scene = bpy.context.scene
    previous_location = None
    previous_target = None
    previous_direction = None
    max_camera_step = 0.0
    max_target_step = 0.0
    max_angular_step = 0.0
    max_camera_step_frame = frame_start
    max_target_step_frame = frame_start
    max_angular_step_frame = frame_start
    max_angular_step_camera = None
    max_angular_step_target = None
    max_angular_step_distance = None
    non_finite_frames = []
    stationary_runs = []
    stationary_start = None
    stationary_epsilon = 1e-5
    for frame in range(frame_start, frame_end + 1):
        scene.frame_set(frame)
        bpy.context.view_layer.update()
        location = camera.matrix_world.translation.copy()
        focus = target.matrix_world.translation.copy()
        values = (*location, *focus)
        if not all(math.isfinite(float(value)) for value in values):
            non_finite_frames.append(frame)
            continue
        direction = (focus - location).normalized()
        if previous_location is not None:
            camera_step = (location - previous_location).length
            target_step = (focus - previous_target).length
            angular_step = math.degrees(previous_direction.angle(direction))
            if camera_step > max_camera_step:
                max_camera_step = camera_step
                max_camera_step_frame = frame
            if target_step > max_target_step:
                max_target_step = target_step
                max_target_step_frame = frame
            if angular_step > max_angular_step:
                max_angular_step = angular_step
                max_angular_step_frame = frame
                max_angular_step_camera = [location.x, location.y, location.z]
                max_angular_step_target = [focus.x, focus.y, focus.z]
                max_angular_step_distance = (focus - location).length
            if (moving_end_frame is None or frame <= moving_end_frame) and camera_step <= stationary_epsilon:
                stationary_start = stationary_start or frame - 1
            elif stationary_start is not None:
                if frame - stationary_start >= 3:
                    stationary_runs.append({"start_frame": stationary_start, "end_frame": frame - 1})
                stationary_start = None
        previous_location = location
        previous_target = focus
        previous_direction = direction
    if stationary_start is not None and frame_end + 1 - stationary_start >= 3:
        stationary_runs.append({"start_frame": stationary_start, "end_frame": min(frame_end, moving_end_frame or frame_end)})
    scene.frame_set(frame_start)
    time_scale = STORY_DURATION_SECONDS / max(duration_seconds, 1e-6)
    fps_scale = 24.0 / max(float(fps), 1.0)
    camera_threshold = 6.0 * time_scale * fps_scale
    target_threshold = 4.0 * time_scale * fps_scale
    angular_threshold = 9.0 * time_scale * fps_scale
    constraint_types = [constraint.type for constraint in camera.constraints]
    failures = []
    if non_finite_frames:
        failures.append({"reason": "non_finite_camera_transform", "frames": non_finite_frames[:20]})
    if max_camera_step > camera_threshold:
        failures.append({"reason": "camera_position_discontinuity", "frame": max_camera_step_frame, "actual": max_camera_step, "maximum": camera_threshold})
    if max_target_step > target_threshold:
        failures.append({"reason": "camera_target_discontinuity", "frame": max_target_step_frame, "actual": max_target_step, "maximum": target_threshold})
    if max_angular_step > angular_threshold:
        failures.append({
            "reason": "camera_rotation_discontinuity",
            "frame": max_angular_step_frame,
            "actual_degrees": max_angular_step,
            "maximum_degrees": angular_threshold,
            "camera_location_mm": max_angular_step_camera,
            "target_location_mm": max_angular_step_target,
            "camera_target_distance_mm": max_angular_step_distance,
        })
    if constraint_types != ["FOLLOW_PATH", "TRACK_TO"]:
        failures.append({"reason": "unexpected_camera_constraints", "actual": constraint_types, "expected": ["FOLLOW_PATH", "TRACK_TO"]})
    if stationary_runs:
        failures.append({"reason": "stationary_camera_before_final_hold", "runs": stationary_runs})
    return {
        "sampled_frames": frame_end - frame_start + 1,
        "single_camera": True,
        "constraint_types": constraint_types,
        "max_camera_step_mm": max_camera_step,
        "max_camera_step_frame": max_camera_step_frame,
        "max_target_step_mm": max_target_step,
        "max_target_step_frame": max_target_step_frame,
        "max_angular_step_degrees": max_angular_step,
        "max_angular_step_frame": max_angular_step_frame,
        "stationary_epsilon_mm": stationary_epsilon,
        "stationary_runs_before_final_hold": stationary_runs,
        "thresholds": {
            "camera_step_mm": camera_threshold,
            "target_step_mm": target_threshold,
            "angular_step_degrees": angular_threshold,
        },
        "failures": failures,
    }


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
        render_samples = 20 if args.smoke_test else (40 if args.render_profile == "review" else 96)
        volumetric_samples = 32 if args.smoke_test else (48 if args.render_profile == "review" else 80)
        scene.eevee.taa_render_samples = render_samples
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
            ("volumetric_samples", volumetric_samples),
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
    dna_points = curve_points("Canonical_DNA_source_path")
    mrna_points = curve_points("Canonical_mRNA_source_path")
    dna_length = polyline_length(dna_points)
    mrna_length = polyline_length(mrna_points)

    collections = {
        "camera": scene_collection("Animation Camera"),
        "labels": scene_collection("Animation Labels"),
        "highlights": scene_collection("Animation Highlights"),
        "lights": scene_collection("Animation Lights"),
        "atmosphere": scene_collection("Animation Atmosphere"),
    }
    mats = {
        "dna_highlight": material("anim_highlight_dna_blue", (0.05, 0.42, 1.0, 0.0), 0.0),
        "rna_highlight": material("anim_highlight_rna_orange", (1.0, 0.44, 0.08, 0.0), 0.0),
        "label": material("anim_label_soft_white", (0.94, 0.965, 1.0, 1.0), 0.34),
    }
    hidden_original_labels = hide_existing_text_labels()
    hidden_original_backdrops = hide_existing_backdrops()
    hidden_original_overview_guides = hide_existing_overview_guides()
    brightened_scale_objects = []
    for obj in bpy.data.objects:
        if obj.type != "CURVE" or not obj.name.startswith("scale_"):
            continue
        obj.data.materials.clear()
        obj.data.materials.append(mats["label"])
        brightened_scale_objects.append(obj.name)
    configure_view_settings()
    compositor_nodes = configure_compositor()

    duration = max(4.0, args.duration_seconds)
    fps = max(1, args.fps)
    frame_end = frame_at(STORY_DURATION_SECONDS, duration, fps)

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

    asset_locations = {asset["name"]: Vector(asset["location_mm"]) for asset in report.get("pdb_assets", [])}
    expected_assets = [asset["name"] for asset in report.get("pdb_assets", [])]
    dna_center = bbox_center(report["dna"]["components"])
    compact_center = bbox_center(report["compact_mrna"]["components"])
    actin = asset_locations["Actin protein"]
    ribosome = (asset_locations["Ribosome small subunit"] + asset_locations["Ribosome large subunit"]) * 0.5
    mrna_mid = point_at_fraction(mrna_points, 0.52)
    mrna_origin_focus = point_at_fraction(mrna_points, 0.06)
    ribosome_approach_focus = point_at_fraction(mrna_points, 0.70)
    ribosome_entry_focus = point_at_fraction(mrna_points, 0.82)
    ribosome_start_focus = point_at_fraction(mrna_points, 0.87)
    ribosome_exit_focus = point_at_fraction(mrna_points, 0.96)
    scene_center = (dna_center + compact_center + actin + ribosome + mrna_mid) / 5.0
    lighting_objects = configure_cinematic_lighting(collections["lights"], camera, scene_center)
    atmosphere = create_volumetric_atmosphere(
        collections["atmosphere"],
        [
            *dna_points,
            *mrna_points,
            dna_center,
            compact_center,
            actin,
            ribosome,
            *asset_locations.values(),
            mrna_mid,
            ribosome_approach_focus,
            ribosome_entry_focus,
            ribosome_start_focus,
            ribosome_exit_focus,
        ],
        density=0.003,
    )

    dna_traversal_offset = Vector((-20.0, -32.0, 24.0))
    mrna_traversal_offset = Vector((-22.0, -30.0, 24.0))
    gene_38 = point_at_fraction(dna_points, 0.38)
    gene_62 = point_at_fraction(dna_points, 0.62)
    gene_84 = point_at_fraction(dna_points, 0.84)
    shots = {
        "overview_start": (scene_center + Vector((-78.0, -168.0, 82.0)), scene_center + Vector((18.0, -2.0, 2.0)), 22.0, 10.0),
        "overview_end": (scene_center + Vector((-74.0, -160.0, 78.0)), scene_center + Vector((16.0, -1.0, 2.0)), 23.0, 9.5),
        "r2r3_myb": (asset_locations["Transcription factor 4"] + dna_traversal_offset, asset_locations["Transcription factor 4"], 50.0, 3.4),
        "cas9": (asset_locations["Cas9"] + dna_traversal_offset, asset_locations["Cas9"], 50.0, 3.2),
        "zbtb24": (asset_locations["Transcription factor 1"] + dna_traversal_offset, asset_locations["Transcription factor 1"], 52.0, 3.0),
        "p53": (asset_locations["p53 tetramer bound to DNA"] + dna_traversal_offset, asset_locations["p53 tetramer bound to DNA"], 64.0, 2.4),
        "foxm1": (asset_locations["Transcription factor 3"] + dna_traversal_offset, asset_locations["Transcription factor 3"], 54.0, 2.8),
        "nucleosome": (asset_locations["Nucleosome"] + dna_traversal_offset, asset_locations["Nucleosome"], 52.0, 3.2),
        "gene_38": (gene_38 + dna_traversal_offset, gene_38, 44.0, 4.0),
        "gene_62": (gene_62 + dna_traversal_offset, gene_62, 44.0, 4.0),
        "gene_84": (gene_84 + dna_traversal_offset, gene_84, 44.0, 4.0),
        "pol": (asset_locations["RNA polymerase II elongation complex"] + dna_traversal_offset, asset_locations["RNA polymerase II elongation complex"], 64.0, 2.4),
        "mrna_origin": (mrna_origin_focus + Vector((-18.0, -28.0, 24.0)), mrna_origin_focus, 56.0, 3.4),
        "pum2": (asset_locations["Pumilio RBP"] + mrna_traversal_offset, asset_locations["Pumilio RBP"], 44.0, 3.0),
        "pabp": (asset_locations["Poly(A)-binding RBP"] + mrna_traversal_offset, asset_locations["Poly(A)-binding RBP"], 44.0, 3.0),
        "ms2": (asset_locations["MS2 coat protein MCP"] + mrna_traversal_offset, asset_locations["MS2 coat protein MCP"], 44.0, 3.0),
        "mcherry": (asset_locations["mCherry/RFP tag"] + mrna_traversal_offset, asset_locations["mCherry/RFP tag"], 44.0, 3.0),
        "argonaute": (asset_locations["Argonaute"] + mrna_traversal_offset, asset_locations["Argonaute"], 44.0, 3.0),
        "hur": (asset_locations["HuR-like RBP"] + mrna_traversal_offset, asset_locations["HuR-like RBP"], 44.0, 3.0),
        "ribosome_approach": (ribosome_approach_focus + mrna_traversal_offset, ribosome_approach_focus, 48.0, 3.2),
        "ribosome_entry": (ribosome_entry_focus + mrna_traversal_offset, ribosome_entry_focus, 48.0, 3.2),
        "ribosome_start": (ribosome_start_focus + mrna_traversal_offset, ribosome_start_focus, 48.0, 3.0),
        "ribosome": (ribosome + mrna_traversal_offset, ribosome, 50.0, 3.0),
        "ribosome_exit": (ribosome_exit_focus + mrna_traversal_offset, ribosome_exit_focus, 48.0, 3.2),
        "compact": (compact_center + Vector((-8.0, -24.0, 20.0)), compact_center, 68.0, 2.4),
        "actin": (actin + Vector((-9.5, -16.0, 11.0)), actin, 90.0, 1.6),
    }
    path_order = [
        "overview_start",
        "overview_end",
        "r2r3_myb",
        "cas9",
        "zbtb24",
        "p53",
        "foxm1",
        "nucleosome",
        "gene_38",
        "gene_62",
        "gene_84",
        "pol",
        "mrna_origin",
        "pum2",
        "pabp",
        "ms2",
        "mcherry",
        "argonaute",
        "hur",
        "ribosome_approach",
        "ribosome_entry",
        "ribosome_start",
        "ribosome",
        "ribosome_exit",
        "compact",
        "actin",
    ]
    camera_path = create_bezier_camera_path(
        "Animation_Flythrough_Bezier_Path",
        [shots[name][0] for name in path_order],
        collections["camera"],
    )
    path_progresses = bezier_path_control_progress(camera_path)
    camera.location = (0.0, 0.0, 0.0)
    follow = camera.constraints.new(type="FOLLOW_PATH")
    follow.target = camera_path
    follow.use_fixed_location = True
    follow.use_curve_follow = False
    track = camera.constraints.new(type="TRACK_TO")
    track.target = target
    track.track_axis = "TRACK_NEGATIVE_Z"
    track.up_axis = "UP_Y"
    bpy.context.scene.camera = camera
    camera.data.dof.focus_object = target

    storyboard_specs = [
        (0.0, "opening_overview_start", 0, "overview_start"),
        (5.0, "opening_overview_glide", 1, "overview_end"),
        (8.0, "r2r3_myb", 2, "r2r3_myb"),
        (9.5, "cas9", 3, "cas9"),
        (11.0, "zbtb24", 4, "zbtb24"),
        (13.2, "p53_slow_pass", 5, "p53"),
        (15.2, "foxm1_dbd", 6, "foxm1"),
        (18.0, "nucleosome", 7, "nucleosome"),
        (20.5, "downstream_gene_38", 8, "gene_38"),
        (22.5, "downstream_gene_62", 9, "gene_62"),
        (24.5, "downstream_gene_84", 10, "gene_84"),
        (27.0, "rna_pol_ii_gene_end", 11, "pol"),
        (29.0, "mrna_origin", 12, "mrna_origin"),
        (31.0, "pum2", 13, "pum2"),
        (33.0, "pabp", 14, "pabp"),
        (35.0, "ms2", 15, "ms2"),
        (37.0, "mcherry", 16, "mcherry"),
        (40.0, "argonaute", 17, "argonaute"),
        (44.0, "hur", 18, "hur"),
        (46.0, "ribosome_approach", 19, "ribosome_approach"),
        (48.0, "ribosome_entry", 20, "ribosome_entry"),
        (50.0, "ribosome_slow_pass_start", 21, "ribosome_start"),
        (52.0, "ribosome_center", 22, "ribosome"),
        (54.0, "ribosome_slow_pass_exit", 23, "ribosome_exit"),
        (57.0, "compact_mrna_pass", 24, "compact"),
        (62.0, "actin_arrival", 25, "actin"),
        (66.0, "actin_final_hold", 25, "actin"),
    ]
    storyboard = []
    for seconds, name, path_index, shot_name in storyboard_specs:
        location, focus, lens, fstop = shots[shot_name]
        frame = frame_at(seconds, duration, fps)
        progress = path_progresses[path_index]
        follow.offset_factor = progress
        camera.data.lens = lens
        camera.data.dof.aperture_fstop = fstop
        target.location = focus
        follow.keyframe_insert(data_path="offset_factor", frame=frame)
        camera.data.keyframe_insert(data_path="lens", frame=frame)
        camera.data.dof.keyframe_insert(data_path="aperture_fstop", frame=frame)
        target.keyframe_insert(data_path="location", frame=frame)
        storyboard.append((seconds, name, progress, location, focus, lens, fstop))

    dna_highlight = create_curve_object("Animation_highlight_DNA_3954bp_path", dna_points, 0.18, mats["dna_highlight"], collections["highlights"])
    mrna_highlight = create_curve_object("Animation_highlight_mRNA_1852nt_path", mrna_points, 0.13, mats["rna_highlight"], collections["highlights"])
    dna_highlight["educational_callout"] = "ACTB promoter + gene DNA; 3954 bp"
    mrna_highlight["educational_callout"] = "actin mRNA; 1852 nt"
    animate_material_window(mats["dna_highlight"], frame_at(5, duration, fps), frame_at(29, duration, fps), 0.64, 0.75)
    animate_material_window(mats["rna_highlight"], frame_at(27, duration, fps), frame_at(57, duration, fps), 0.68, 0.85)

    highlight_objects = [dna_highlight.name, mrna_highlight.name]
    emission_specs = [
        ("p53 tetramer bound to DNA", "protein", 12.0, 14.2, 0.45),
        ("Nucleosome", "protein", 16.8, 19.4, 0.30),
        ("RNA polymerase II elongation complex", "protein", 25.5, 29.2, 0.35),
        ("Ribosome small subunit", "protein", 50.0, 53.5, 0.10),
        ("Ribosome large subunit", "protein", 50.0, 53.5, 0.10),
        ("Standalone tRNA", "nucleic", 50.0, 53.5, 0.16),
        ("Actin protein", "protein", 58.0, 62.0, 0.40),
    ]
    emission_targets = [
        pulse_asset_emission(
            asset_name,
            component,
            frame_at(start, duration, fps),
            frame_at(end, duration, fps),
            peak,
        )
        for asset_name, component, start, end, peak in emission_specs
    ]

    asset_label_specs = [
        ("Transcription factor 4", "r2r3_myb", "R2R3 MYB", "r2r3_myb", -4.5, 3.4, 0.56, 7.3, 9.2),
        ("Cas9", "cas9", "Cas9", "cas9", 4.5, 3.2, 0.56, 8.7, 10.5),
        ("Transcription factor 1", "zbtb24", "ZBTB24", "zbtb24", -4.5, 3.4, 0.56, 10.0, 12.2),
        ("p53 tetramer bound to DNA", "p53", "p53 tetramer", "p53", 4.5, 3.0, 0.56, 11.7, 14.4),
        ("Transcription factor 3", "foxm1", "FOXM1-DBD", "foxm1", -4.5, 3.4, 0.56, 14.0, 16.5),
        ("Nucleosome", "nucleosome", "nucleosome", "nucleosome", -5.0, 3.4, 0.58, 16.5, 19.5),
        ("RNA polymerase II elongation complex", "rna_pol_ii", "RNA Pol II at gene end", "pol", 5.0, 3.4, 0.58, 25.5, 29.3),
        ("Pumilio RBP", "pum2", "PUM2", "pum2", -5.0, 4.0, 0.60, 29.5, 32.2),
        ("Poly(A)-binding RBP", "pabp", "PABP", "pabp", 5.0, 4.0, 0.60, 31.5, 34.2),
        ("MS2 coat protein MCP", "ms2", "MS2 coat protein", "ms2", -5.0, 3.8, 0.56, 33.5, 36.2),
        ("mCherry/RFP tag", "mcherry", "mCherry", "mcherry", 5.0, 3.8, 0.56, 35.5, 38.5),
        ("Argonaute", "argonaute", "Argonaute", "argonaute", -5.0, 3.8, 0.58, 38.0, 42.2),
        ("HuR-like RBP", "hur", "HuR", "hur", 5.0, 3.8, 0.58, 42.0, 45.8),
        ("Ribosome small subunit", "ribosome_small", "ribosome small subunit", "ribosome", -6.5, 5.2, 0.50, 50.0, 53.0),
        ("Ribosome large subunit", "ribosome_large", "ribosome large subunit", "ribosome", 6.0, 1.0, 0.50, 50.0, 53.0),
        ("Standalone tRNA", "trna", "tRNA", "ribosome", 9.0, -5.0, 0.54, 50.0, 53.0),
        ("Actin protein", "actin", "ACTB protein\n375 aa", "actin", 1.7, 1.2, 0.46, 59.0, 66.0),
    ]
    molecule_label_specs = [
        ("dna", "ACTB promoter + gene DNA\n3,954 bp", point_at_fraction(dna_points, 0.03), "r2r3_myb", -6.0, 5.5, 0.82, 5.0, 8.0),
        ("mrna", "Actin mRNA\n1,852 nt", mrna_origin_focus, "mrna_origin", 5.0, 4.0, 0.72, 27.0, 30.0),
        ("compact", "compact mRNA", compact_center, "compact", -5.0, 2.2, 0.58, 54.0, 59.0),
    ]
    label_records = []
    label_texts = []
    asset_label_records = []
    for asset_name, name, text, shot_name, screen_x, screen_y, size, start, end in asset_label_specs:
        anchor = asset_locations[asset_name]
        shot = shots[shot_name]
        location = camera_plane_label_location(anchor, shot[0], shot[1], screen_x, screen_y)
        start_frame = frame_at(start, duration, fps)
        end_frame = max(start_frame, frame_at(end, duration, fps) - (0 if end >= STORY_DURATION_SECONDS else 1))
        obj = create_label(
            f"Animation_{name}",
            text,
            location,
            size,
            mats["label"],
            collections["labels"],
            camera,
            start_frame,
            end_frame,
        )
        leader = create_label_leader(
            f"Animation_{name}_leader",
            anchor,
            location,
            mats["label"],
            collections["labels"],
            start_frame,
            end_frame,
        )
        record = {
            "object": obj.name,
            "leader": leader.name,
            "asset": asset_name,
            "text": text,
            "start_frame": start_frame,
            "end_frame": end_frame,
            "kind": "asset",
        }
        asset_label_records.append(record)
        label_records.append(record)
        label_texts.append(text)

    for name, text, anchor, shot_name, screen_x, screen_y, size, start, end in molecule_label_specs:
        shot = shots[shot_name]
        location = camera_plane_label_location(anchor, shot[0], shot[1], screen_x, screen_y)
        start_frame = frame_at(start, duration, fps)
        end_frame = max(start_frame, frame_at(end, duration, fps) - (0 if end >= STORY_DURATION_SECONDS else 1))
        obj = create_label(
            f"Animation_label_{name}",
            text,
            location,
            size,
            mats["label"],
            collections["labels"],
            camera,
            start_frame,
            end_frame,
        )
        record = {"object": obj.name, "text": text, "start_frame": start_frame, "end_frame": end_frame, "kind": "molecule"}
        label_records.append(record)
        label_texts.append(text)

    scale_label_records = []
    scale_label_end = max(1, frame_at(5.0, duration, fps) - 1)
    for name, row in report.get("scale_bars", {}).items():
        origin = Vector(row["origin_mm"])
        location = origin + Vector((float(row["length_mm"]) + 2.4, 0.0, 0.4))
        text = row["label"]
        obj = create_label(
            f"Animation_label_{name}",
            text,
            location,
            1.90,
            mats["label"],
            collections["labels"],
            camera,
            1,
            scale_label_end,
        )
        obj.data.align_x = "LEFT"
        record = {"object": obj.name, "text": text, "start_frame": 1, "end_frame": scale_label_end, "kind": "scale"}
        scale_label_records.append(record)
        label_records.append(record)
        label_texts.append(text)

    set_animation_interpolation()
    overview_start = 1
    overview_end = max(overview_start, frame_at(5.0, duration, fps) - 1)
    overview_specs = [
        ("mrna", "Actin mRNA — 1,852 nt", (5.25, 1.15, -18.0), 0.31, (4.75, -0.35, 2.65)),
        ("dna", "ACTB promoter + gene DNA — 3,954 bp", (5.25, -3.15, -18.0), 0.27, (4.75, -4.45, -2.00)),
    ]
    overview_labels = []
    protein_label = create_camera_label(
        "Animation_overview_label_protein",
        "ACTB protein — 375 aa",
        camera,
        (5.25, 4.05, -18.0),
        0.31,
        mats["label"],
        collections["labels"],
        overview_start,
        overview_end,
    )
    protein_label.data.align_x = "LEFT"
    protein_pointer = create_animated_camera_pointer(
        "Animation_overview_pointer_protein",
        camera,
        actin,
        Vector((4.85, 4.05, -18.0)),
        -18.0,
        mats["label"],
        collections["labels"],
        overview_start,
        overview_end,
    )
    protein_row = {
        "label": protein_label.name,
        "pointer": protein_pointer.name,
        "bracket": None,
        "text": "ACTB protein — 375 aa",
        "start_frame": overview_start,
        "end_frame": overview_end,
    }
    overview_labels.append(protein_row)
    label_records.append({"object": protein_label.name, "text": protein_row["text"], "start_frame": overview_start, "end_frame": overview_end, "kind": "overview"})
    label_texts.append(protein_row["text"])
    for name, text, local_location, size, bracket in overview_specs:
        label = create_camera_label(
            f"Animation_overview_label_{name}",
            text,
            camera,
            local_location,
            size,
            mats["label"],
            collections["labels"],
            overview_start,
            overview_end,
        )
        label.data.align_x = "LEFT"
        bracket_obj = create_camera_bracket(
            f"Animation_overview_bracket_{name}",
            camera,
            *bracket,
            -18.0,
            mats["label"],
            collections["labels"],
            overview_start,
            overview_end,
        )
        row = {
            "label": label.name,
            "bracket": bracket_obj.name,
            "pointer": None,
            "text": text,
            "start_frame": overview_start,
            "end_frame": overview_end,
        }
        overview_labels.append(row)
        label_records.append({"object": label.name, "text": text, "start_frame": overview_start, "end_frame": overview_end, "kind": "overview"})
        label_texts.append(text)

    configure_render(args, frame_end, frames_dir)
    set_animation_interpolation()
    actin_arrival_frame = frame_at(62.0, duration, fps)
    motion_validation = camera_motion_continuity(camera, target, 1, frame_end, duration, fps, moving_end_frame=actin_arrival_frame)
    if motion_validation["failures"]:
        raise RuntimeError(f"Animation camera continuity validation failed: {motion_validation['failures']}")
    moving_progresses = [progress for seconds, _name, progress, *_rest in storyboard if seconds <= 62.0]
    non_increasing_progress = [
        {"index": index, "previous": moving_progresses[index - 1], "actual": moving_progresses[index]}
        for index in range(1, len(moving_progresses))
        if moving_progresses[index] <= moving_progresses[index - 1]
    ]
    if non_increasing_progress:
        raise RuntimeError(f"Camera path progress is not strictly increasing before the actin hold: {non_increasing_progress}")
    labeled_assets = [record["asset"] for record in asset_label_records]
    missing_asset_labels = sorted(set(expected_assets) - set(labeled_assets))
    duplicate_asset_labels = sorted({name for name in labeled_assets if labeled_assets.count(name) > 1})
    if missing_asset_labels or duplicate_asset_labels or len(labeled_assets) != 17:
        raise RuntimeError(
            f"Individual asset-label coverage failed: missing={missing_asset_labels}, duplicates={duplicate_asset_labels}, count={len(labeled_assets)}"
        )
    forbidden_grouped_labels = [text for text in label_texts if any(token in text for token in ("PUM2 + PABP", "Argonaute + HuR", "ribosome + tRNA"))]
    if forbidden_grouped_labels:
        raise RuntimeError(f"Grouped animation labels remain: {forbidden_grouped_labels}")
    pointer_validation = camera_pointer_projection_validation(protein_pointer, camera, actin, overview_start, overview_end)
    if not pointer_validation["passed"]:
        raise RuntimeError(f"Actin overview pointer is not aligned: {pointer_validation}")
    mm_text = [text for text in label_texts if "mm" in text.lower()]
    if mm_text:
        raise RuntimeError(f"Animation labels still contain millimeter dimensions: {mm_text}")
    duplicate_glow_objects = [obj.name for obj in bpy.data.objects if obj.name.startswith("Animation_pulse_")]
    if duplicate_glow_objects:
        raise RuntimeError(f"Duplicate glow geometry remains in animation: {duplicate_glow_objects}")

    output_blend.parent.mkdir(parents=True, exist_ok=True)
    bpy.context.preferences.filepaths.save_version = 0
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
        "story_duration_seconds": STORY_DURATION_SECONDS,
        "fps": fps,
        "frame_start": 1,
        "frame_end": frame_end,
        "resolution": [args.resolution_x, args.resolution_y],
        "render_profile": args.render_profile,
        "smoke_test": bool(args.smoke_test),
        "camera": camera.name,
        "camera_type": camera.data.type,
        "target": target.name,
        "camera_path": camera_path.name,
        "camera_path_points": [[point.x, point.y, point.z] for point in [shots[name][0] for name in path_order]],
        "camera_motion_continuity": motion_validation,
        "source_path_lengths_mm": {"dna": dna_length, "mrna": mrna_length},
        "storyboard": [
            {
                "story_seconds": seconds,
                "timeline_seconds": round(seconds * duration / STORY_DURATION_SECONDS, 3),
                "frame": frame_at(seconds, duration, fps),
                "name": name,
                "path_progress": progress,
                "camera_location_mm": [location.x, location.y, location.z],
                "focus_mm": [focus.x, focus.y, focus.z],
                "lens_mm": lens,
                "dof_aperture_fstop": fstop,
            }
            for seconds, name, progress, location, focus, lens, fstop in storyboard
        ],
        "labels": label_records,
        "scale_bar_labels": scale_label_records,
        "brightened_scale_objects": brightened_scale_objects,
        "overview_labels": overview_labels,
        "overview_actin_pointer_validation": pointer_validation,
        "individual_asset_labels": asset_label_records,
        "individual_asset_label_coverage": {
            "expected_count": len(expected_assets),
            "actual_count": len(labeled_assets),
            "expected_assets": expected_assets,
            "labeled_assets": labeled_assets,
            "missing": missing_asset_labels,
            "duplicates": duplicate_asset_labels,
            "passed": not missing_asset_labels and not duplicate_asset_labels and len(labeled_assets) == 17,
        },
        "path_progress_validation": {
            "strictly_increasing_until_story_second": 62.0,
            "non_increasing_steps": non_increasing_progress,
            "passed": not non_increasing_progress,
        },
        "text_validation": {"millimeter_labels": mm_text, "passed": not mm_text},
        "actin_hold_seconds": round(4.0 * duration / STORY_DURATION_SECONDS, 3),
        "emission_targets": emission_targets,
        "lighting_objects": lighting_objects,
        "atmosphere_object": atmosphere.name,
        "atmosphere_density": atmosphere["atmosphere_density"],
        "compositor_nodes": compositor_nodes,
        "hidden_original_labels": hidden_original_labels,
        "hidden_original_backdrops": hidden_original_backdrops,
        "hidden_original_overview_guides": hidden_original_overview_guides,
        "highlight_objects": highlight_objects,
    }


def main() -> None:
    args = parse_args()
    args.source_blend = as_path(args.source_blend)
    args.source_report = as_path(args.source_report)
    args.output_dir = as_path(args.output_dir)
    if args.smoke_test:
        args.render_profile = "smoke"
        args.duration_seconds = min(args.duration_seconds, 8.0)
        args.fps = min(args.fps, 6)
        args.resolution_x = min(args.resolution_x, 640)
        args.resolution_y = min(args.resolution_y, 360)

    output_blend = as_path(args.output_blend) if args.output_blend else args.output_dir / "flythrough_animation.blend"
    default_mp4 = {
        "smoke": "flythrough_animation_smoke.mp4",
        "review": "flythrough_animation_review.mp4",
        "final": "flythrough_animation_1080p.mp4",
    }[args.render_profile]
    output_mp4 = as_path(args.output_mp4) if args.output_mp4 else args.output_dir / default_mp4
    report_path = args.output_dir / "flythrough_animation_report.json"
    default_frames_dir = {"smoke": "frames_smoke", "review": "frames_review", "final": "frames"}[args.render_profile]
    frames_dir = as_path(args.frames_dir) if args.frames_dir else args.output_dir / default_frames_dir

    if not args.source_blend.exists():
        raise FileNotFoundError(f"Canonical blend not found: {args.source_blend}")
    if not args.source_report.exists():
        raise FileNotFoundError(f"Canonical report not found: {args.source_report}")

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
