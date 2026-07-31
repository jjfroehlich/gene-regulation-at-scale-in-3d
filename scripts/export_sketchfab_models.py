#!/usr/bin/env python3
"""Export upload-ready canonical models for Sketchfab.

Run with Blender in background mode after building the canonical scene.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import bpy


ROOT = Path(__file__).resolve().parents[1]
EXPORT_DIR = ROOT / "outputs" / "sketchfab"


def output_paths() -> tuple[Path, Path, Path]:
    stem = "gene_expression_canonical_sketchfab"
    return EXPORT_DIR / f"{stem}.glb", EXPORT_DIR / f"{stem}.fbx", EXPORT_DIR / f"{stem}_export_report.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def select_molecular_meshes() -> list[bpy.types.Object]:
    bpy.ops.object.select_all(action="DESELECT")
    selected = []
    excluded_collections = {"Labels", "Scale bars", "Detail context"}
    for obj in bpy.context.scene.objects:
        collection_names = {collection.name for collection in obj.users_collection}
        include = (
            obj.type == "MESH"
            and not obj.hide_render
            and not obj.get("canonical_beauty_object")
            and not obj.get("canonical_detail_title")
            and collection_names.isdisjoint(excluded_collections)
        )
        obj.select_set(include)
        if include:
            selected.append(obj)
    if not selected:
        raise RuntimeError("No molecular meshes selected for Sketchfab export")
    bpy.context.view_layer.objects.active = selected[0]
    return selected


def main() -> None:
    glb_path, fbx_path, report_path = output_paths()
    EXPORT_DIR.mkdir(parents=True, exist_ok=True)
    selected = select_molecular_meshes()
    bpy.ops.export_scene.gltf(
        filepath=str(glb_path),
        export_format="GLB",
        use_selection=True,
        export_apply=True,
        export_cameras=False,
        export_lights=False,
        export_materials="EXPORT",
    )
    bpy.ops.export_scene.fbx(
        filepath=str(fbx_path),
        use_selection=True,
        object_types={"MESH"},
        use_mesh_modifiers=True,
        add_leaf_bones=False,
        bake_anim=False,
        apply_unit_scale=True,
        apply_scale_options="FBX_SCALE_UNITS",
        axis_forward="-Z",
        axis_up="Y",
    )
    report = {
        "source_blend": bpy.data.filepath,
        "selection_policy": "visible molecular meshes only; excludes cameras, lights, backdrop, labels, scale bars, and detail-only context",
        "scale_policy": "one common scene scale preserved across all molecular objects",
        "selected_object_count": len(selected),
        "selected_objects": [obj.name for obj in selected],
        "exports": {},
    }
    for name, path in (("glb_preferred", glb_path), ("fbx_fallback", fbx_path)):
        report["exports"][name] = {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
