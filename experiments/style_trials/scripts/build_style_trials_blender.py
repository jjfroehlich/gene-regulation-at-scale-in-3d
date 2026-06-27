#!/usr/bin/env python3
"""Build a Blender contact sheet comparing molecular representation styles."""

from __future__ import annotations

import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

import bpy
from mathutils import Vector


ROOT = Path(os.environ.get("GENE_SCENE_ROOT", Path(__file__).resolve().parents[3])).resolve()
sys.path.insert(0, str(ROOT / "scripts"))
import build_gene_expression_scene as base  # noqa: E402


EXPERIMENT_DIR = ROOT / "experiments" / "style_trials"
EXPORT_DIR = EXPERIMENT_DIR / "assets" / "pymol_exports" / "style_trials"
RCSB_DIR = ROOT / "assets" / "rcsb"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
BLEND_PATH = OUTPUT_DIR / "style_trials.blend"
PREVIEW_PATH = OUTPUT_DIR / "preview_style_trials.png"
REPORT_PATH = OUTPUT_DIR / "style_trials_report.json"
ANGSTROM_TO_MM = 0.04


STRUCTURES = {
    "1NKP": {
        "title": "TF-DNA complex (1NKP)",
        "x": -20.0,
        "bead_radius_mm": 0.10,
        "bead_max": 700,
        "protein_color": "tf_surface",
        "nucleic_color": "dna_surface",
    },
    "1KX5": {
        "title": "Nucleosome (1KX5)",
        "x": 24.0,
        "bead_radius_mm": 0.11,
        "bead_max": 1400,
        "protein_color": "histone_surface",
        "nucleic_color": "dna_surface",
    },
}

ROWS = {
    "surface": {"title": "PyMOL surface", "y": 22.0},
    "spheres": {"title": "Atom spheres", "y": 0.0},
    "beads": {"title": "Residue beads", "y": -22.0},
}

COLORS = {
    "background": (0.98, 0.98, 0.96, 1.0),
    "text": (0.12, 0.12, 0.12, 1.0),
    "muted_text": (0.42, 0.42, 0.40, 1.0),
    "dna_surface": (0.94, 0.49, 0.18, 1.0),
    "tf_surface": (0.17, 0.70, 0.72, 1.0),
    "histone_surface": (0.47, 0.70, 0.30, 1.0),
    "sphere_protein": (0.46, 0.55, 0.86, 1.0),
    "sphere_nucleic": (0.88, 0.57, 0.25, 1.0),
}


def clean_scene() -> None:
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete()
    for collection in list(bpy.data.collections):
        bpy.data.collections.remove(collection)
    for mesh in list(bpy.data.meshes):
        if mesh.users == 0:
            bpy.data.meshes.remove(mesh)


def configure_scene() -> None:
    scene = bpy.context.scene
    scene.unit_settings.system = "METRIC"
    scene.unit_settings.length_unit = "MILLIMETERS"
    scene.unit_settings.scale_length = 0.001
    scene.render.engine = "BLENDER_WORKBENCH"
    scene.display.shading.light = "STUDIO"
    scene.display.shading.color_type = "MATERIAL"
    scene.display.shading.background_type = "VIEWPORT"
    scene.display.shading.background_color = COLORS["background"][:3]
    scene.render.resolution_x = 2200
    scene.render.resolution_y = 1600
    scene.view_settings.view_transform = "Standard"
    scene.view_settings.look = "None"
    scene.view_settings.exposure = 0.0
    scene.view_settings.gamma = 1.0
    world = scene.world or bpy.data.worlds.new("World")
    scene.world = world
    world.color = COLORS["background"][:3]
    world.use_nodes = True
    background = world.node_tree.nodes.get("Background")
    if background:
        background.inputs["Color"].default_value = COLORS["background"]
        background.inputs["Strength"].default_value = 1.0


def make_collection(name: str) -> bpy.types.Collection:
    collection = bpy.data.collections.new(name)
    bpy.context.scene.collection.children.link(collection)
    return collection


def make_materials() -> dict[str, bpy.types.Material]:
    materials = {}
    for name, color in COLORS.items():
        mat = bpy.data.materials.new(name)
        mat.diffuse_color = color
        mat.use_nodes = True
        bsdf = mat.node_tree.nodes.get("Principled BSDF")
        if bsdf:
            bsdf.inputs["Base Color"].default_value = color
            bsdf.inputs["Roughness"].default_value = 0.88
            bsdf.inputs["Metallic"].default_value = 0.0
            if "Alpha" in bsdf.inputs:
                bsdf.inputs["Alpha"].default_value = color[3]
        materials[name] = mat
    return materials


def link_to_collection(obj: bpy.types.Object, collection: bpy.types.Collection) -> None:
    collection.objects.link(obj)
    for user_collection in list(obj.users_collection):
        if user_collection != collection:
            user_collection.objects.unlink(obj)


def create_label(
    text: str,
    location: tuple[float, float, float],
    size: float,
    material: bpy.types.Material,
    collection: bpy.types.Collection,
    align: str = "CENTER",
) -> bpy.types.Object:
    curve = bpy.data.curves.new(f"label_{text}", "FONT")
    curve.body = text
    curve.align_x = align
    curve.align_y = "CENTER"
    curve.size = size
    curve.materials.append(material)
    obj = bpy.data.objects.new(curve.name, curve)
    obj.location = location
    link_to_collection(obj, collection)
    return obj


def parse_pymol_obj(path: Path, scale: float) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]]]:
    vertices = []
    faces = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            if raw_line.startswith("v "):
                _, x, y, z = raw_line.split()[:4]
                vertices.append((float(x) * scale, float(y) * scale, float(z) * scale))
            elif raw_line.startswith("f "):
                face = []
                for token in raw_line.split()[1:]:
                    face.append(int(token.split("/")[0]) - 1)
                if len(face) >= 3:
                    faces.append(tuple(face))
    return vertices, faces


def import_obj_mesh(
    path: Path,
    name: str,
    material: bpy.types.Material,
    collection: bpy.types.Collection,
) -> bpy.types.Object:
    vertices, faces = parse_pymol_obj(path, ANGSTROM_TO_MM)
    mesh = bpy.data.meshes.new(f"{name}_mesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    mesh.materials.append(material)
    obj = bpy.data.objects.new(name, mesh)
    obj["source_obj"] = str(path.relative_to(ROOT)).replace("\\", "/")
    obj["angstrom_to_mm"] = ANGSTROM_TO_MM
    link_to_collection(obj, collection)
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)
    bpy.ops.object.shade_smooth()
    obj.select_set(False)
    return obj


def bounds(objects: list[bpy.types.Object]) -> tuple[Vector, Vector]:
    coords = []
    for obj in objects:
        for vertex in obj.data.vertices:
            coords.append(obj.matrix_world @ vertex.co)
    min_v = Vector((min(c.x for c in coords), min(c.y for c in coords), min(c.z for c in coords)))
    max_v = Vector((max(c.x for c in coords), max(c.y for c in coords), max(c.z for c in coords)))
    return min_v, max_v


def move_group_to(objects: list[bpy.types.Object], target: Vector) -> None:
    min_v, max_v = bounds(objects)
    center = (min_v + max_v) * 0.5
    delta = target - center
    for obj in objects:
        obj.location += delta


def import_pymol_style(
    pdb_id: str,
    style: str,
    target: Vector,
    materials: dict[str, bpy.types.Material],
    collection: bpy.types.Collection,
) -> dict:
    structure = STRUCTURES[pdb_id]
    style_dir = EXPORT_DIR / pdb_id / style
    protein_path = style_dir / f"{pdb_id}_{style}_protein.obj"
    nucleic_path = style_dir / f"{pdb_id}_{style}_nucleic.obj"
    material_prefix = "sphere" if style == "spheres" else None
    protein_material = materials["sphere_protein"] if material_prefix else materials[structure["protein_color"]]
    nucleic_material = materials["sphere_nucleic"] if material_prefix else materials[structure["nucleic_color"]]
    objects = [
        import_obj_mesh(protein_path, f"{pdb_id}_{style}_protein", protein_material, collection),
        import_obj_mesh(nucleic_path, f"{pdb_id}_{style}_nucleic", nucleic_material, collection),
    ]
    move_group_to(objects, target)
    min_v, max_v = bounds(objects)
    return {
        "pdb_id": pdb_id,
        "style": style,
        "objects": [obj.name for obj in objects],
        "bbox_mm": list(max_v - min_v),
        "source_files": [str(protein_path.relative_to(ROOT)), str(nucleic_path.relative_to(ROOT))],
    }


def build_bead_baseline(
    pdb_id: str,
    target: Vector,
    materials: dict[str, bpy.types.Material],
    collection: bpy.types.Collection,
) -> dict:
    structure = STRUCTURES[pdb_id]
    atoms = base.parse_atom_site(RCSB_DIR / f"{pdb_id}.cif")
    points = base.residue_points(atoms)
    sampled = base.sample_points(points, structure["bead_max"])
    vectors = [Vector(point["pos_A"]) for point in points]
    center_A = Vector((
        (min(v.x for v in vectors) + max(v.x for v in vectors)) * 0.5,
        (min(v.y for v in vectors) + max(v.y for v in vectors)) * 0.5,
        (min(v.z for v in vectors) + max(v.z for v in vectors)) * 0.5,
    ))
    by_kind = defaultdict(list)
    for point in sampled:
        local = (Vector(point["pos_A"]) - center_A) * ANGSTROM_TO_MM
        by_kind[point["kind"]].append(target + local)
    objects = []
    for kind, centers in by_kind.items():
        material = materials[structure["nucleic_color"] if kind == "nucleic" else structure["protein_color"]]
        obj = base.create_bead_mesh(
            f"{pdb_id}_residue_beads_{kind}",
            centers,
            structure["bead_radius_mm"],
            material,
            collection,
        )
        if obj:
            objects.append(obj)
            obj["pdb_id"] = pdb_id
            obj["style"] = "residue_beads"
            obj["angstrom_to_mm"] = ANGSTROM_TO_MM
    min_v, max_v = bounds(objects)
    return {
        "pdb_id": pdb_id,
        "style": "residue_beads",
        "objects": [obj.name for obj in objects],
        "sampled_beads": len(sampled),
        "bbox_mm": list(max_v - min_v),
    }


def create_atom_sphere_mesh(
    name: str,
    centers_and_radii: list[tuple[Vector, float]],
    material: bpy.types.Material,
    collection: bpy.types.Collection,
) -> bpy.types.Object | None:
    if not centers_and_radii:
        return None
    base_verts, base_faces = base.icosahedron()
    verts = []
    faces = []
    for center, radius in centers_and_radii:
        offset = len(verts)
        verts.extend([tuple(center + v * radius) for v in base_verts])
        faces.extend([(a + offset, b + offset, c + offset) for a, b, c in base_faces])
    mesh = bpy.data.meshes.new(f"{name}_mesh")
    mesh.from_pydata(verts, [], faces)
    mesh.update()
    mesh.materials.append(material)
    obj = bpy.data.objects.new(name, mesh)
    link_to_collection(obj, collection)
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)
    bpy.ops.object.shade_smooth()
    obj.select_set(False)
    return obj


def atom_radius_mm(element: str) -> float:
    radius_A = {
        "C": 1.70,
        "N": 1.55,
        "O": 1.52,
        "P": 1.80,
        "S": 1.80,
    }.get(element.upper(), 1.60)
    return radius_A * ANGSTROM_TO_MM


def build_atom_spheres(
    pdb_id: str,
    target: Vector,
    materials: dict[str, bpy.types.Material],
    collection: bpy.types.Collection,
) -> dict:
    atoms = [atom for atom in base.parse_atom_site(RCSB_DIR / f"{pdb_id}.cif") if atom["element"] != "H"]
    vectors = [Vector((atom["x"], atom["y"], atom["z"])) for atom in atoms]
    center_A = Vector((
        (min(v.x for v in vectors) + max(v.x for v in vectors)) * 0.5,
        (min(v.y for v in vectors) + max(v.y for v in vectors)) * 0.5,
        (min(v.z for v in vectors) + max(v.z for v in vectors)) * 0.5,
    ))
    by_kind = defaultdict(list)
    for atom in atoms:
        kind = "nucleic" if atom["comp"] in base.NUCLEIC_COMPS else "protein"
        center = target + (Vector((atom["x"], atom["y"], atom["z"])) - center_A) * ANGSTROM_TO_MM
        by_kind[kind].append((center, atom_radius_mm(atom["element"])))
    objects = []
    for kind, centers_and_radii in by_kind.items():
        material = materials["sphere_nucleic" if kind == "nucleic" else "sphere_protein"]
        obj = create_atom_sphere_mesh(f"{pdb_id}_atom_spheres_{kind}", centers_and_radii, material, collection)
        if obj:
            objects.append(obj)
            obj["pdb_id"] = pdb_id
            obj["style"] = "atom_spheres"
            obj["angstrom_to_mm"] = ANGSTROM_TO_MM
    min_v, max_v = bounds(objects)
    return {
        "pdb_id": pdb_id,
        "style": "atom_spheres",
        "objects": [obj.name for obj in objects],
        "atom_count": len(atoms),
        "bbox_mm": list(max_v - min_v),
    }


def add_camera_and_light() -> None:
    camera_data = bpy.data.cameras.new("Camera")
    camera = bpy.data.objects.new("Camera", camera_data)
    bpy.context.scene.collection.objects.link(camera)
    camera.location = (0.0, 0.0, 110.0)
    camera.rotation_euler = (0.0, 0.0, 0.0)
    camera_data.type = "ORTHO"
    camera_data.ortho_scale = 84.0
    bpy.context.scene.camera = camera

    light_data = bpy.data.lights.new("softbox", "AREA")
    light = bpy.data.objects.new("softbox", light_data)
    bpy.context.scene.collection.objects.link(light)
    light.location = (0.0, 0.0, 80.0)
    light_data.energy = 350.0
    light_data.size = 80.0


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    clean_scene()
    configure_scene()
    materials = make_materials()
    collection = make_collection("Style trials")
    label_collection = make_collection("Labels")
    report = {"angstrom_to_mm": ANGSTROM_TO_MM, "entries": []}

    for pdb_id, meta in STRUCTURES.items():
        create_label(meta["title"], (meta["x"], 32.5, 0.2), 1.6, materials["text"], label_collection)
    for row in ROWS.values():
        create_label(row["title"], (-25.5, row["y"], 0.2), 1.2, materials["muted_text"], label_collection, align="RIGHT")

    for pdb_id, meta in STRUCTURES.items():
        surface_target = Vector((meta["x"], ROWS["surface"]["y"], 0.0))
        report["entries"].append(import_pymol_style(pdb_id, "surface", surface_target, materials, collection))
        atom_target = Vector((meta["x"], ROWS["spheres"]["y"], 0.0))
        report["entries"].append(build_atom_spheres(pdb_id, atom_target, materials, collection))
        bead_target = Vector((meta["x"], ROWS["beads"]["y"], 0.0))
        report["entries"].append(build_bead_baseline(pdb_id, bead_target, materials, collection))

    create_label(
        "All imported PyMOL OBJ geometry is scaled at 1 A = 0.04 mm",
        (0.0, -34.5, 0.2),
        1.35,
        materials["muted_text"],
        label_collection,
    )
    add_camera_and_light()
    bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_PATH))
    bpy.context.scene.render.filepath = str(PREVIEW_PATH)
    bpy.ops.render.render(write_still=True)
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
