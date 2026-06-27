#!/usr/bin/env python3
"""Compare direct Blender DNA/RNA proxy meshes against the current PyMOL proxy style."""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import bpy
from mathutils import Vector


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import blender_nucleic_meshes as direct_meshes  # noqa: E402
import build_gene_expression_scene as base  # noqa: E402
import build_gene_expression_surface_scene as surface_scene  # noqa: E402
import procedural_nucleic_geometry as geom  # noqa: E402


EXPERIMENT_DIR = ROOT / "experiments" / "direct_blender_nucleic_acids"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs"
BLEND_PATH = OUTPUT_DIR / "direct_blender_nucleic_acid_proxy_comparison.blend"
PREVIEW_PATH = OUTPUT_DIR / "preview_direct_blender_nucleic_acid_proxy_comparison.png"
DNA_PREVIEW_PATH = OUTPUT_DIR / "preview_direct_blender_nucleic_acid_dna_closeup.png"
RNA_PREVIEW_PATH = OUTPUT_DIR / "preview_direct_blender_nucleic_acid_rna_closeup.png"
REPORT_PATH = OUTPUT_DIR / "direct_blender_nucleic_acid_proxy_comparison_report.json"


def experiment_manifest() -> dict:
    manifest = copy.deepcopy(base.load_manifest())
    manifest["title"] = f"{manifest['title']} - direct Blender nucleic-acid proxy experiment"
    dna = manifest["procedural_nucleic_acids"]["dna"]
    dna.update(
        {
            "style": "direct_blender_surface_proxy_polished",
            "strand_sides": 12,
            "strand_radius_mm": 0.12,
            "base_pair_every_bp": 1,
            "base_pair_radius_mm": 0.082,
            "base_pair_sides": 10,
            "base_stack_radius_mm": 0.172,
            "base_stack_sides": 12,
            "surface_bump_amplitude": 0.050,
            "direct_voxel_size_mm": 0.052,
            "direct_smooth_factor": 0.035,
            "direct_smooth_iterations": 1,
            "direct_unified_mesh": True,
        }
    )
    mrna = manifest["procedural_nucleic_acids"]["mrna"]
    mrna.update(
        {
            "style": "direct_blender_surface_proxy_polished",
            "tube_sides": 12,
            "tube_radius_mm": 0.165,
            "surface_bump_amplitude": 0.10,
            "base_lobe_every_nt": 1,
            "base_lobe_radius_mm": 0.095,
            "base_lobe_sides": 8,
            "base_offset_mm": 0.17,
            "nucleotide_detail_every_nt": 1,
            "phosphate_radius_mm": 0.09,
            "sugar_radius_mm": 0.085,
            "rna_base_ellipsoid_radii_mm": [0.135, 0.078, 0.055],
            "base_connector_radius_mm": 0.048,
            "nucleotide_detail_rings": 5,
            "nucleotide_detail_sides": 8,
            "direct_voxel_size_mm": 0.05,
            "direct_smooth_factor": 0.028,
            "direct_smooth_iterations": 1,
            "direct_unified_mesh": True,
        }
    )
    return manifest


def make_collection(name: str) -> bpy.types.Collection:
    collection = bpy.data.collections.new(name)
    bpy.context.scene.collection.children.link(collection)
    return collection


def shift_collection(collection: bpy.types.Collection, delta: Vector) -> None:
    for obj in collection.objects:
        obj.location += delta


def add_label(text: str, location: tuple[float, float, float], materials: dict, collection: bpy.types.Collection) -> None:
    base.create_text(f"label_{text.replace(' ', '_')}", text, location, 3.0, materials["black"], collection)


def add_cameras() -> dict[str, str]:
    specs = {
        "full": ("Camera_direct_blender_full", (0.0, 0.0, 285.0), 455.0),
        "dna": ("Camera_direct_blender_dna_closeup", (-125.0, -45.0, 85.0), 34.0),
        "rna": ("Camera_direct_blender_rna_closeup", (-48.0, 37.0, 95.0), 48.0),
    }
    for _key, (name, location, ortho_scale) in specs.items():
        camera_data = bpy.data.cameras.new(name)
        camera = bpy.data.objects.new(name, camera_data)
        camera.location = location
        camera.rotation_euler = (0.0, 0.0, 0.0)
        camera_data.type = "ORTHO"
        camera_data.ortho_scale = ortho_scale
        bpy.context.scene.collection.objects.link(camera)
    bpy.context.scene.camera = bpy.data.objects[specs["full"][0]]
    return {key: value[0] for key, value in specs.items()}


def build_pymol_comparison(manifest: dict, path_reports: dict, collection: bpy.types.Collection, materials: dict, offset: Vector) -> dict:
    dna_center = surface_scene.path_center(path_reports["direct_dna_path"]) + offset
    mrna_center = surface_scene.path_center(path_reports["direct_mrna_path"]) + offset
    compact_center = surface_scene.path_center(path_reports["direct_compact_mrna_path"]) + offset
    return {
        "dna": surface_scene.build_polished_nucleic_proxy(
            "DNA_PROXY",
            "DNA PyMOL proxy comparison polished",
            surface_scene.DNA_PROXY_COMPONENTS,
            materials,
            collection,
            "dna",
            dna_center,
        ),
        "mrna": surface_scene.build_polished_nucleic_proxy(
            "MRNA_PROXY",
            "mRNA PyMOL proxy comparison polished",
            surface_scene.MRNA_PROXY_COMPONENTS,
            materials,
            collection,
            "rna",
            mrna_center,
        ),
        "compact_mrna": surface_scene.build_polished_nucleic_proxy(
            "MRNA_COMPACT_PROXY",
            "compact mRNA PyMOL proxy comparison polished",
            surface_scene.MRNA_PROXY_COMPONENTS,
            materials,
            collection,
            "rna",
            compact_center,
        ),
    }


def validate(report: dict) -> None:
    dna = report["direct_blender"]["dna"]
    mrna = report["direct_blender"]["mrna"]
    compact = report["direct_blender"]["compact_mrna"]
    if report["pipeline"]["pymol_export_invoked"]:
        raise RuntimeError("Direct Blender experiment unexpectedly invoked PyMOL export")
    if dna["handedness"] != "right_handed":
        raise RuntimeError("Direct Blender DNA is not reported as right-handed")
    if abs(float(dna["bp_spacing_mm"]) - 0.136) > 1e-9:
        raise RuntimeError("DNA bp spacing is not 0.136 mm/bp")
    if abs(float(dna["estimated_envelope_diameter_mm"]) - 0.8) > 0.08:
        raise RuntimeError("DNA envelope diameter is outside tolerance")
    for label, item in (("elongated", mrna), ("compact", compact)):
        if abs(float(item["total_measured_mm"]) - 222.24) > 0.05:
            raise RuntimeError(f"{label} mRNA length is not 222.24 mm")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = experiment_manifest()

    surface_scene.clean_scene()
    surface_scene.configure_scene()
    direct_collection = make_collection("Direct Blender DNA_RNA")
    pymol_collection = make_collection("PyMOL proxy comparison")
    label_collection = make_collection("Labels")
    materials = base.make_materials()
    surface_scene.soften_materials(materials)

    direct_collections = {"DNA": direct_collection, "mRNA": direct_collection}
    direct_offset = Vector((-125.0, 0.0, 0.0))
    pymol_offset = Vector((125.0, 0.0, 0.0))

    dna_path, dna_report = direct_meshes.build_dna_meshes(manifest, direct_collections, materials)
    mrna_path, mrna_report = direct_meshes.build_mrna_meshes(manifest, direct_collections, materials)
    compact_path, compact_report = direct_meshes.build_compact_mrna_meshes(manifest, direct_collections, materials)
    shift_collection(direct_collection, direct_offset)

    pymol_reports = build_pymol_comparison(
        manifest,
        {
            "direct_dna_path": dna_path,
            "direct_mrna_path": mrna_path,
            "direct_compact_mrna_path": compact_path,
        },
        pymol_collection,
        materials,
        pymol_offset,
    )

    add_label("Direct Blender proxy", (-125.0, 84.0, 0.4), materials, label_collection)
    add_label("Current PyMOL proxy", (125.0, 84.0, 0.4), materials, label_collection)
    camera_names = add_cameras()

    report = {
        "title": manifest["title"],
        "kind": "direct_blender_nucleic_acid_experiment",
        "units": manifest["units"],
        "outputs": {
            "blend": str(BLEND_PATH),
            "preview": str(PREVIEW_PATH),
            "dna_closeup": str(DNA_PREVIEW_PATH),
            "rna_closeup": str(RNA_PREVIEW_PATH),
            "report": str(REPORT_PATH),
        },
        "pipeline": {
            "direct_blender_route": "procedural pseudoatom-derived centers to Blender mesh primitives, voxel remesh, weighted normals",
            "pymol_export_invoked": False,
            "pymol_proxy_comparison_uses_cached_reduced_assets": True,
        },
        "direct_blender": {
            "panel_offset_mm": list(direct_offset),
            "dna": dna_report,
            "mrna": mrna_report,
            "compact_mrna": compact_report,
        },
        "pymol_proxy_comparison": {
            "panel_offset_mm": list(pymol_offset),
            **pymol_reports,
        },
        "visual_comparison_notes": {
            "intended_decision": "Judge whether direct Blender polish is close enough to the PyMOL surface proxy for later editable workflows.",
            "expected_tradeoff": "Direct Blender geometry is more editable and does not need PyMOL, but may still show more component seams than PyMOL proxy surfaces.",
        },
    }
    validate(report)
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    bpy.ops.wm.save_as_mainfile(filepath=str(BLEND_PATH))
    for key, path in (("full", PREVIEW_PATH), ("dna", DNA_PREVIEW_PATH), ("rna", RNA_PREVIEW_PATH)):
        bpy.context.scene.camera = bpy.data.objects[camera_names[key]]
        bpy.context.scene.render.filepath = str(path)
        bpy.ops.render.render(write_still=True)
    print(f"Wrote {BLEND_PATH}")
    print(f"Wrote {PREVIEW_PATH}")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
