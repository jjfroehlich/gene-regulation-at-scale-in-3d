#!/usr/bin/env python3
"""Build Molecular Nodes DNA/RNA comparison V2.

This stays experimental. The canonical scene remains PyMOL surfaces plus our
own procedural nucleic-acid pipeline unless Molecular Nodes proves scale-correct
and path-controllable without empirical correction.
"""

from __future__ import annotations

import importlib
import json
import math
import os
import sys
from pathlib import Path

import addon_utils
import bpy
from mathutils import Vector


ROOT = Path(os.environ.get("GENE_SCENE_ROOT", Path(__file__).resolve().parents[3])).resolve()
sys.path.insert(0, str(ROOT / "scripts"))
import build_gene_expression_scene as base  # noqa: E402


EXPERIMENT_DIR = ROOT / "experiments" / "molecular_nodes"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
RCSB_DIR = ROOT / "assets" / "rcsb"
GENERATED_RNA_CIF = (
    ROOT
    / "experiments"
    / "procedural_nucleic_acids"
    / "assets"
    / "detail_route_pymol_proxy"
    / "cif"
    / "RNA_DETAIL_PROXY.cif"
)
BLEND_PATH = OUTPUT_DIR / "molecular_nodes_dna_rna_comparison_v2.blend"
PREVIEW_PATH = OUTPUT_DIR / "preview_molecular_nodes_dna_rna_comparison_v2.png"
REPORT_PATH = OUTPUT_DIR / "molecular_nodes_dna_rna_comparison_v2_report.json"
ANGSTROM_TO_MM = 0.04
DNA_BP_TO_MM = 0.136
RNA_NT_TO_MM = 0.12
BACKGROUND_COLOR = (0.985, 0.985, 0.965, 1.0)


def clean_scene() -> None:
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete()
    for collection in list(bpy.data.collections):
        bpy.data.collections.remove(collection)


def configure_scene() -> None:
    scene = bpy.context.scene
    scene.unit_settings.system = "METRIC"
    scene.unit_settings.length_unit = "MILLIMETERS"
    scene.unit_settings.scale_length = 0.001
    scene.render.engine = "BLENDER_WORKBENCH"
    scene.render.resolution_x = 2200
    scene.render.resolution_y = 1250
    scene.display.shading.light = "STUDIO"
    scene.display.shading.color_type = "MATERIAL"
    scene.display.shading.background_type = "VIEWPORT"
    scene.display.shading.background_color = BACKGROUND_COLOR[:3]
    scene.view_settings.view_transform = "Standard"
    world = scene.world or bpy.data.worlds.new("World")
    scene.world = world
    world.color = BACKGROUND_COLOR[:3]


def material(name: str, color: tuple[float, float, float, float]) -> bpy.types.Material:
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = color
    return mat


def add_text(name: str, body: str, location: tuple[float, float, float], size: float, mat: bpy.types.Material) -> None:
    curve = bpy.data.curves.new(name, "FONT")
    curve.body = body
    curve.align_x = "CENTER"
    curve.align_y = "CENTER"
    curve.size = size
    curve.materials.append(mat)
    obj = bpy.data.objects.new(name, curve)
    obj.location = location
    bpy.context.scene.collection.objects.link(obj)


def bounds(objects: list[bpy.types.Object]) -> tuple[float, float, float]:
    coords = []
    for obj in objects:
        if obj.type != "MESH":
            continue
        for corner in obj.bound_box:
            coords.append(obj.matrix_world @ Vector(corner))
    if not coords:
        return (0.0, 0.0, 0.0)
    return (
        max(c.x for c in coords) - min(c.x for c in coords),
        max(c.y for c in coords) - min(c.y for c in coords),
        max(c.z for c in coords) - min(c.z for c in coords),
    )


def atom_count(cif_path: Path) -> int:
    atoms = base.parse_atom_site(cif_path)
    return len(atoms)


def expected_bbox_mm(cif_path: Path) -> tuple[float, float, float]:
    atoms = base.parse_atom_site(cif_path)
    points = base.residue_points(atoms)
    vectors = [Vector(point["pos_A"]) for point in points]
    if not vectors:
        return (0.0, 0.0, 0.0)
    return (
        (max(v.x for v in vectors) - min(v.x for v in vectors)) * ANGSTROM_TO_MM,
        (max(v.y for v in vectors) - min(v.y for v in vectors)) * ANGSTROM_TO_MM,
        (max(v.z for v in vectors) - min(v.z for v in vectors)) * ANGSTROM_TO_MM,
    )


def interface_summary(group_name: str) -> dict:
    group = bpy.data.node_groups.get(group_name)
    if group is None:
        return {"status": "missing"}
    return {
        "status": "available",
        "inputs": [
            {"name": item.name, "socket_type": item.socket_type}
            for item in group.interface.items_tree
            if getattr(item, "item_type", "") == "SOCKET" and item.in_out == "INPUT"
        ],
        "outputs": [
            {"name": item.name, "socket_type": item.socket_type}
            for item in group.interface.items_tree
            if getattr(item, "item_type", "") == "SOCKET" and item.in_out == "OUTPUT"
        ],
    }


def set_input_default(node: bpy.types.Node, name: str, value) -> bool:
    socket = node.inputs.get(name)
    if socket is None:
        return False
    try:
        socket.default_value = value
        return True
    except Exception:  # noqa: BLE001 - report failure in caller.
        return False


def add_geometry_socket(tree: bpy.types.GeometryNodeTree, name: str, in_out: str) -> None:
    tree.interface.new_socket(name=name, in_out=in_out, socket_type="NodeSocketGeometry")


def create_input_curve(name: str, bp_count: int, location: tuple[float, float, float]) -> bpy.types.Object:
    curve = bpy.data.curves.new(name, "CURVE")
    curve.dimensions = "3D"
    curve.resolution_u = 12
    spline = curve.splines.new("POLY")
    samples = 80
    spline.points.add(samples - 1)
    length_mm = bp_count * DNA_BP_TO_MM
    for i, point in enumerate(spline.points):
        t = i / (samples - 1)
        x = (t - 0.5) * length_mm
        y = math.sin(t * math.tau * 1.25) * 0.35
        point.co = (x, y, 0.0, 1.0)
    obj = bpy.data.objects.new(name, curve)
    obj.location = location
    bpy.context.scene.collection.objects.link(obj)
    return obj


def evaluated_mesh_report(obj: bpy.types.Object) -> dict:
    depsgraph = bpy.context.evaluated_depsgraph_get()
    bpy.context.view_layer.update()
    evaluated = obj.evaluated_get(depsgraph)
    mesh = bpy.data.meshes.new_from_object(evaluated, depsgraph=depsgraph)
    coords = [obj.matrix_world @ vertex.co for vertex in mesh.vertices]
    if coords:
        bbox = (
            max(v.x for v in coords) - min(v.x for v in coords),
            max(v.y for v in coords) - min(v.y for v in coords),
            max(v.z for v in coords) - min(v.z for v in coords),
        )
        center = (
            (max(v.x for v in coords) + min(v.x for v in coords)) * 0.5,
            (max(v.y for v in coords) + min(v.y for v in coords)) * 0.5,
            (max(v.z for v in coords) + min(v.z for v in coords)) * 0.5,
        )
    else:
        bbox = (0.0, 0.0, 0.0)
        center = (0.0, 0.0, 0.0)
    report = {
        "evaluated_vertices": len(mesh.vertices),
        "evaluated_faces": len(mesh.polygons),
        "bbox_mm": bbox,
        "center_mm": center,
    }
    bpy.data.meshes.remove(mesh)
    return report


def build_mn_dna_curve_trial(text_mat: bpy.types.Material) -> dict:
    from bl_ext.blender_org.molecularnodes.nodes import nodes

    appended = []
    for name in ["MN_dna_double_helix", "MN_dna_bases", "MN_dna_style_surface"]:
        if bpy.data.node_groups.get(name) is None:
            nodes.append(name)
        appended.append(name)

    input_bp = 42
    target_x = -8.0
    curve_obj = create_input_curve("MN V2 input curve 42 bp", input_bp, (target_x, 1.1, 0.0))
    surface_mat = material("mn_dna_nodes_surface_taupe", (0.58, 0.45, 0.31, 1.0))
    tree = bpy.data.node_groups.new("MN V2 DNA curve surface wrapper", "GeometryNodeTree")
    add_geometry_socket(tree, "Geometry", "INPUT")
    add_geometry_socket(tree, "Geometry", "OUTPUT")
    group_input = tree.nodes.new("NodeGroupInput")
    group_output = tree.nodes.new("NodeGroupOutput")
    bases_node = tree.nodes.new("GeometryNodeGroup")
    bases_node.node_tree = bpy.data.node_groups["MN_dna_bases"]
    helix_node = tree.nodes.new("GeometryNodeGroup")
    helix_node.node_tree = bpy.data.node_groups["MN_dna_double_helix"]
    surface_node = tree.nodes.new("GeometryNodeGroup")
    surface_node.node_tree = bpy.data.node_groups["MN_dna_style_surface"]
    realize_node = tree.nodes.new("GeometryNodeRealizeInstances")
    group_input.location = (-780, 120)
    bases_node.location = (-780, -160)
    helix_node.location = (-470, 40)
    surface_node.location = (-160, 40)
    realize_node.location = (90, 40)
    group_output.location = (340, 40)

    input_collection = bpy.data.collections.get("prim_DNA")
    configured_inputs = {}
    if input_collection is not None:
        if input_collection.name not in {collection.name for collection in bpy.context.scene.collection.children}:
            bpy.context.scene.collection.children.link(input_collection)
        configured_inputs["prim_DNA_collection"] = set_input_default(bases_node, "Collection", input_collection)
    configured_inputs["base_colors"] = {
        "dA": set_input_default(bases_node, "dA", (0.82, 0.38, 0.16, 1.0)),
        "dC": set_input_default(bases_node, "dC", (0.74, 0.58, 0.28, 1.0)),
        "dG": set_input_default(bases_node, "dG", (0.54, 0.36, 0.21, 1.0)),
        "dT": set_input_default(bases_node, "dT", (0.90, 0.58, 0.32, 1.0)),
        "Backbone": set_input_default(bases_node, "Backbone", (0.55, 0.36, 0.23, 1.0)),
    }
    configured_inputs["style_surface"] = {
        "Resolution": set_input_default(surface_node, "Resolution", 3),
        "Radius": set_input_default(surface_node, "Radius", 1.0),
        "Probe Size": set_input_default(surface_node, "Probe Size", 0.4),
        "Subdivision Surface": set_input_default(surface_node, "Subdivision Surface", 0),
        "Color by CA": set_input_default(surface_node, "Color by CA", False),
        "Interpolate Color": set_input_default(surface_node, "Interpolate Color", 1),
        "Shade Smooth": set_input_default(surface_node, "Shade Smooth", True),
        "Material": set_input_default(surface_node, "Material", surface_mat),
    }

    tree.links.new(group_input.outputs["Geometry"], helix_node.inputs["Curve"])
    tree.links.new(bases_node.outputs["Bases"], helix_node.inputs["Bases"])
    tree.links.new(helix_node.outputs["Base Instances"], surface_node.inputs["Bases"])
    tree.links.new(surface_node.outputs["Bases"], realize_node.inputs["Geometry"])
    tree.links.new(realize_node.outputs["Geometry"], group_output.inputs["Geometry"])
    modifier = curve_obj.modifiers.new("MN DNA curve surface", "NODES")
    modifier.node_group = tree
    curve_obj["style_source"] = "Molecular Nodes DNA nodes"
    curve_obj["input_bp_count"] = input_bp
    curve_obj["input_curve_length_mm"] = input_bp * DNA_BP_TO_MM

    mesh_report = evaluated_mesh_report(curve_obj)
    add_text("label_mn_dna_nodes", "MN DNA nodes from curve", (target_x, 6.0, 0.0), 0.34, text_mat)
    return {
        "route": "MN_dna_double_helix + MN_dna_bases + MN_dna_style_surface",
        "appended_node_groups": appended,
        "input_curve_object": curve_obj.name,
        "input_bp_count": input_bp,
        "expected_axis_length_mm": input_bp * DNA_BP_TO_MM,
        "configured_inputs": configured_inputs,
        "created_prim_DNA_collection": input_collection is not None,
        "realize_instances_after_surface": True,
        "evaluated_geometry": mesh_report,
        "path_controllable": mesh_report["evaluated_vertices"] > 0,
        "scale_assessment": (
            "curve_length_is_scale_correct_input; Molecular Nodes output scale still requires visual/bbox review"
            if mesh_report["evaluated_vertices"] > 0
            else "blocked_no_evaluated_geometry_from_documented_node_chain"
        ),
    }


def set_style_surface_inputs(objects: list[bpy.types.Object], params: dict) -> dict:
    applied = {}
    for obj in objects:
        for modifier in obj.modifiers:
            tree = getattr(modifier, "node_group", None)
            if tree is None:
                continue
            for node in tree.nodes:
                node_tree = getattr(node, "node_tree", None)
                if node_tree is None or not node_tree.name.startswith("Style Surface"):
                    continue
                for name, value in params.items():
                    applied[f"{obj.name}:{name}"] = set_input_default(node, name, value)
    return applied


def style_surface_socket_values(objects: list[bpy.types.Object]) -> dict:
    values = {}
    for obj in objects:
        for modifier in obj.modifiers:
            tree = getattr(modifier, "node_group", None)
            if tree is None:
                continue
            for node in tree.nodes:
                node_tree = getattr(node, "node_tree", None)
                if node_tree is None or not node_tree.name.startswith("Style Surface"):
                    continue
                values[obj.name] = {}
                for socket in node.inputs:
                    if hasattr(socket, "default_value"):
                        try:
                            values[obj.name][socket.name] = socket.default_value
                        except Exception:  # noqa: BLE001
                            values[obj.name][socket.name] = "<unreadable>"
    return values


def import_mn_style_surface(
    label: str,
    cif_path: Path,
    target_location: tuple[float, float, float],
    params: dict | None,
    text_mat: bpy.types.Material,
    label_text: str,
) -> dict:
    before = {obj.name for obj in bpy.data.objects}
    result = bpy.ops.mn.import_local(
        filepath=str(cif_path.resolve()),
        style="surface",
        node_setup=True,
        centre=True,
        centre_type="mass",
        remove_solvent=True,
        assembly=False,
    )
    created = [obj for obj in bpy.data.objects if obj.name not in before]
    raw_bbox = bounds(created)
    expected_bbox = expected_bbox_mm(cif_path)
    raw_max = max(raw_bbox) if raw_bbox else 0.0
    expected_max = max(expected_bbox) if expected_bbox else 0.0
    scale_factor = expected_max / raw_max if raw_max > 0 and expected_max > 0 else 1.0
    applied_params = {}
    if params:
        applied_params = set_style_surface_inputs(created, params)
    for obj in created:
        obj.scale = (scale_factor, scale_factor, scale_factor)
        obj.location.x += target_location[0]
        obj.location.y += target_location[1]
        obj.location.z += target_location[2]
        obj["style_source"] = "Molecular Nodes Style Surface"
        obj["molecular_nodes_empirical_bbox_scale_factor"] = scale_factor
        obj["comparison_label"] = label
    bpy.context.view_layer.update()
    add_text(f"label_{label}", label_text, (target_location[0], 6.0, 0.0), 0.31, text_mat)
    return {
        "label": label,
        "source_cif": str(cif_path.relative_to(ROOT)).replace("\\", "/"),
        "operator_result": list(result),
        "object_count": len(created),
        "created_objects": [obj.name for obj in created],
        "atom_count": atom_count(cif_path),
        "raw_molecular_nodes_bbox": raw_bbox,
        "expected_bbox_mm_from_cif": expected_bbox,
        "applied_empirical_bbox_scale_factor": scale_factor,
        "bbox_after_scale_mm": bounds(created),
        "style_surface_parameters_requested": params or "default",
        "style_surface_parameters_applied": applied_params,
        "style_surface_socket_values": style_surface_socket_values(created),
        "scale_assessment": "requires_empirical_bbox_scale_factor_for_scene_mm",
    }


def probe_rna_specific_nodes() -> dict:
    module = importlib.import_module("bl_ext.blender_org.molecularnodes")
    nodes_yml = Path(module.__file__).resolve().parent / "assets" / "nodes.yml"
    mentions = []
    if nodes_yml.exists():
        for line in nodes_yml.read_text(encoding="utf-8", errors="ignore").splitlines():
            stripped = line.strip()
            if stripped.startswith("name:") and "rna" in stripped.lower():
                mentions.append(stripped.split(":", 1)[1].strip())
    return {
        "nodes_yml": str(nodes_yml),
        "rna_named_node_groups": sorted(set(mentions)),
        "rna_specific_procedural_node_callable": False,
        "assessment": "No RNA-specific procedural path node was found; RNA test uses Style Surface on generated and real RNA coordinates.",
    }


def add_camera() -> None:
    light_data = bpy.data.lights.new("softbox", "AREA")
    light_data.energy = 900.0
    light_data.size = 60.0
    light = bpy.data.objects.new("softbox", light_data)
    light.location = (0.0, -5.0, 65.0)
    bpy.context.scene.collection.objects.link(light)

    camera_data = bpy.data.cameras.new("Camera")
    camera = bpy.data.objects.new("Camera", camera_data)
    camera.location = (1.0, 0.0, 72.0)
    camera.rotation_euler = (0.0, 0.0, 0.0)
    camera_data.type = "ORTHO"
    camera_data.ortho_scale = 18.0
    bpy.context.scene.collection.objects.link(camera)
    bpy.context.scene.camera = camera


def disable_molecular_nodes_save_handlers(report: dict) -> list:
    removed = []
    for handler in list(bpy.app.handlers.save_post):
        module = getattr(handler, "__module__", "")
        if "molecularnodes" in module.lower():
            bpy.app.handlers.save_post.remove(handler)
            removed.append(handler)
    if removed:
        report["save_notes"] = {
            "molecular_nodes_save_handlers_temporarily_removed": [
                f"{getattr(handler, '__module__', '')}.{getattr(handler, '__name__', repr(handler))}"
                for handler in removed
            ],
            "reason": "Keep this comparison blend save independent from Molecular Nodes session hooks.",
        }
    return removed


def restore_molecular_nodes_save_handlers(handlers: list) -> None:
    for handler in handlers:
        if handler not in bpy.app.handlers.save_post:
            bpy.app.handlers.save_post.append(handler)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    clean_scene()
    configure_scene()
    text_mat = material("text_grey", (0.22, 0.22, 0.22, 1.0))
    report = {
        "title": "Molecular Nodes DNA/RNA comparison V2",
        "canonical_scene_changed": False,
        "units": {
            "angstrom_to_mm": ANGSTROM_TO_MM,
            "dna_bp_to_mm": DNA_BP_TO_MM,
            "rna_nt_to_mm": RNA_NT_TO_MM,
        },
        "sources": {
            "molecular_nodes_dna_docs": "https://bradyajohnston.github.io/MolecularNodes/nodes/DNA.html",
            "molecular_nodes_style_docs": "https://bradyajohnston.github.io/MolecularNodes/nodes/style.html",
        },
        "node_group_interfaces": {},
        "entries": {},
        "errors": [],
        "outputs": {
            "blend": str(BLEND_PATH.relative_to(ROOT)).replace("\\", "/"),
            "preview": str(PREVIEW_PATH.relative_to(ROOT)).replace("\\", "/"),
            "report": str(REPORT_PATH.relative_to(ROOT)).replace("\\", "/"),
        },
    }
    try:
        enabled = addon_utils.check("bl_ext.blender_org.molecularnodes")[1]
        if not enabled:
            addon_utils.enable("bl_ext.blender_org.molecularnodes", default_set=False, persistent=False)
        report["addon_enabled"] = addon_utils.check("bl_ext.blender_org.molecularnodes")[1]
        from bl_ext.blender_org.molecularnodes.nodes import nodes

        for name in ["MN_dna_double_helix", "MN_dna_bases", "MN_dna_style_surface", "Style Surface"]:
            if bpy.data.node_groups.get(name) is None:
                nodes.append(name)
            report["node_group_interfaces"][name] = interface_summary(name)
        report["entries"]["mn_dna_curve_surface"] = build_mn_dna_curve_trial(text_mat)
        report["entries"]["mn_style_surface_1BNA_default"] = import_mn_style_surface(
            "mn_1BNA_default",
            RCSB_DIR / "1BNA.cif",
            (-3.5, 1.1, 0.0),
            None,
            text_mat,
            "MN Style Surface 1BNA default",
        )
        report["entries"]["mn_style_surface_1BNA_tuned"] = import_mn_style_surface(
            "mn_1BNA_tuned",
            RCSB_DIR / "1BNA.cif",
            (2.2, 1.1, 0.0),
            {"Quality": 4, "Scale Radius": 1.2, "Probe Size": 0.5, "Relaxation Steps": 8, "Color Blur": 1, "Shade Smooth": True},
            text_mat,
            "MN Style Surface 1BNA tuned",
        )
        report["entries"]["mn_style_surface_1EHZ_tRNA"] = import_mn_style_surface(
            "mn_1EHZ_tRNA",
            RCSB_DIR / "1EHZ.cif",
            (8.0, 1.1, 0.0),
            {"Quality": 3, "Scale Radius": 1.15, "Probe Size": 0.5, "Relaxation Steps": 8, "Color Blur": 1, "Shade Smooth": True},
            text_mat,
            "MN Style Surface 1EHZ tRNA",
        )
        report["entries"]["mn_style_surface_generated_rna_proxy"] = import_mn_style_surface(
            "mn_generated_rna",
            GENERATED_RNA_CIF,
            (13.5, 1.1, 0.0),
            {"Quality": 3, "Scale Radius": 1.05, "Probe Size": 0.4, "Relaxation Steps": 6, "Color Blur": 1, "Shade Smooth": True},
            text_mat,
            "MN Style Surface generated RNA",
        )
        report["entries"]["rna_specific_node_probe"] = probe_rna_specific_nodes()
        report["entries"]["heavy_rna_calibrator_note"] = {
            "pdb_id": "9IOB",
            "local_atom_count": atom_count(RCSB_DIR / "9IOB.cif") if (RCSB_DIR / "9IOB.cif").exists() else None,
            "status": "not_imported_in_molecular_nodes_v2",
            "reason": "9IOB has ~82k atoms locally; this pass uses 1EHZ as the small real-RNA Style Surface calibrator and the procedural V2 pass already includes a PyMOL 9IOB calibrator.",
        }
        add_text("label_mn_summary", "Molecular Nodes V2: DNA nodes and Style Surface trials", (1.0, -6.7, 0.0), 0.35, text_mat)
        report["acceptance_checks"] = {
            "scale_constants_ok": ANGSTROM_TO_MM == 0.04 and DNA_BP_TO_MM == 0.136 and RNA_NT_TO_MM == 0.12,
            "documented_dna_node_groups_available": all(
                report["node_group_interfaces"].get(name, {}).get("status") == "available"
                for name in ["MN_dna_double_helix", "MN_dna_bases", "MN_dna_style_surface"]
            ),
            "style_surface_node_available": report["node_group_interfaces"].get("Style Surface", {}).get("status") == "available",
            "dna_nodes_path_controllable": report["entries"]["mn_dna_curve_surface"].get("path_controllable", False),
            "rna_specific_procedural_node_callable": report["entries"]["rna_specific_node_probe"].get("rna_specific_procedural_node_callable", False),
            "canonical_suitability": "experimental_only_pending_scale_and_path_control_review",
        }
    except Exception as exc:  # noqa: BLE001 - trial should still produce a file/report.
        report["errors"].append(repr(exc))
        add_text("label_mn_error", f"Molecular Nodes V2 failed: {exc!r}", (0.0, 0.0, 0.0), 0.45, text_mat)
    add_camera()
    removed_handlers = disable_molecular_nodes_save_handlers(report)
    try:
        REPORT_PATH.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
        bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_PATH))
    finally:
        restore_molecular_nodes_save_handlers(removed_handlers)
    bpy.context.scene.render.filepath = str(PREVIEW_PATH)
    bpy.ops.render.render(write_still=True)
    print(f"Wrote {BLEND_PATH}")
    print(f"Wrote {PREVIEW_PATH}")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
