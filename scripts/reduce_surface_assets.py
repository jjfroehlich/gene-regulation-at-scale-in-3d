#!/usr/bin/env python3
"""Reduce PyMOL OBJ surface assets for responsive Blender scene assembly."""

from __future__ import annotations

import json
import os
from pathlib import Path

import bmesh
import bpy


ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "assets" / "pymol_exports" / "surface_assets"
REDUCED_DIR = ROOT / "assets" / "pymol_exports" / "surface_assets_reduced"
REPORT_PATH = REDUCED_DIR / "surface_assets_reduced_manifest.json"

RIBOSOME_VALID_COMPONENTS = {
    "ribosomal_protein",
    "rRNA",
    "mRNA",
    "tRNA_AV",
    "tRNA_AW",
    "tRNA_AY",
}

RIBOSOME_TARGETS = {
    "ribosomal_protein": 80_000,
    "protein": 80_000,
    "rRNA": 120_000,
    "nucleic": 120_000,
    "mRNA": 8_000,
}

RIBOSOME_PDBS = {"1J5E", "1JJ2", "4V5D"}

GENERIC_TARGETS = {
    "protein": 80_000,
    "nucleic": 30_000,
}

PROCEDURAL_TARGETS = {
    "DNA_PROXY": {
        "strand_A": 55_000,
        "strand_B": 55_000,
        "base_pairs": 90_000,
    },
    "MRNA_PROXY": {
        "utr5": 15_000,
        "coding": 90_000,
        "utr3": 60_000,
    },
    "MRNA_COMPACT_PROXY": {
        "utr5": 20_000,
        "coding": 105_000,
        "utr3": 70_000,
    },
}


def procedural_target_key(pdb_id: str) -> str | None:
    for key in PROCEDURAL_TARGETS:
        if pdb_id == key or pdb_id.startswith(f"{key}_"):
            return key
    return None


def parse_component(path: Path) -> tuple[str, str] | None:
    pdb_id = path.parent.name.upper()
    prefix = f"{pdb_id}_surface_"
    if not path.name.startswith(prefix) or path.suffix.lower() != ".obj":
        return None
    if procedural_target_key(pdb_id) and os.environ.get("SURFACE_REDUCTION_INCLUDE_PROCEDURAL", "") != "1":
        return None
    return pdb_id, path.stem[len(prefix):]


def iter_raw_assets() -> list[tuple[str, str, Path]]:
    requested = selected_ids()
    assets = []
    for path in sorted(RAW_DIR.glob("*/*_surface_*.obj")):
        parsed = parse_component(path)
        if parsed is None or path.stat().st_size == 0:
            continue
        pdb_id, component = parsed
        if requested is not None and pdb_id not in requested:
            continue
        if pdb_id == "4V5D" and component not in RIBOSOME_VALID_COMPONENTS:
            continue
        assets.append((pdb_id, component, path))
    return assets


def selected_ids() -> set[str] | None:
    raw = os.environ.get("SURFACE_REDUCTION_IDS", "").strip()
    if not raw:
        return None
    return {item.strip().upper() for item in raw.split(",") if item.strip()}


def target_faces(pdb_id: str, component: str, raw_faces: int) -> int:
    procedural_key = procedural_target_key(pdb_id)
    if procedural_key:
        return min(raw_faces, PROCEDURAL_TARGETS[procedural_key].get(component, raw_faces))
    if pdb_id in RIBOSOME_PDBS:
        if component.startswith("tRNA_"):
            return min(raw_faces, 25_000)
        return min(raw_faces, RIBOSOME_TARGETS.get(component, raw_faces))
    return min(raw_faces, GENERIC_TARGETS.get(component, raw_faces))


def parse_obj(path: Path) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]]]:
    vertices = []
    faces = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            if raw_line.startswith("v "):
                _, x, y, z = raw_line.split()[:4]
                vertices.append((float(x), float(y), float(z)))
            elif raw_line.startswith("f "):
                face = []
                for token in raw_line.split()[1:]:
                    index = token.split("/")[0]
                    if index:
                        face.append(int(index) - 1)
                if len(face) >= 3:
                    faces.append(tuple(face))
    return vertices, faces


def write_obj(path: Path, mesh: bpy.types.Mesh) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    used_indices = sorted({index for polygon in mesh.polygons for index in polygon.vertices})
    remap = {old: new for new, old in enumerate(used_indices, start=1)}
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# Reduced PyMOL surface OBJ\n")
        for index in used_indices:
            vertex = mesh.vertices[index]
            co = vertex.co
            handle.write(f"v {co.x:.6f} {co.y:.6f} {co.z:.6f}\n")
        for polygon in mesh.polygons:
            indices = " ".join(str(remap[index]) for index in polygon.vertices)
            handle.write(f"f {indices}\n")


def create_mesh(name: str, vertices: list[tuple[float, float, float]], faces: list[tuple[int, ...]]) -> bpy.types.Object:
    mesh = bpy.data.meshes.new(f"{name}_mesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.scene.collection.objects.link(obj)
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)
    return obj


def weld_duplicate_vertices(mesh: bpy.types.Mesh, distance: float = 0.0001) -> int:
    """PyMOL OBJ files often duplicate vertices per face; weld them before decimation."""
    before = len(mesh.vertices)
    bm = bmesh.new()
    bm.from_mesh(mesh)
    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=distance)
    bm.to_mesh(mesh)
    bm.free()
    mesh.update()
    return before - len(mesh.vertices)


def reduce_one(pdb_id: str, component: str, raw_path: Path) -> dict:
    print(f"Reducing {pdb_id} {component}: {raw_path}")
    vertices, faces = parse_obj(raw_path)
    raw_face_count = len(faces)
    raw_vertex_count = len(vertices)
    target = target_faces(pdb_id, component, raw_face_count)
    out_path = REDUCED_DIR / pdb_id / raw_path.name

    obj = create_mesh(f"{pdb_id}_{component}_reduce", vertices, faces)
    welded_vertices = weld_duplicate_vertices(obj.data)
    ratio = 1.0
    if raw_face_count > target and raw_face_count > 0:
        ratio = max(0.001, target / raw_face_count)
        modifier = obj.modifiers.new("component_face_target", "DECIMATE")
        modifier.decimate_type = "COLLAPSE"
        modifier.ratio = ratio
        bpy.ops.object.modifier_apply(modifier=modifier.name)
    for polygon in obj.data.polygons:
        polygon.use_smooth = True
    write_obj(out_path, obj.data)

    reduced_faces = len(obj.data.polygons)
    reduced_vertices = len({index for polygon in obj.data.polygons for index in polygon.vertices})
    bpy.data.objects.remove(obj, do_unlink=True)
    return {
        "pdb_id": pdb_id,
        "component": component,
        "raw_obj": str(raw_path.relative_to(ROOT)).replace("\\", "/"),
        "reduced_obj": str(out_path.relative_to(ROOT)).replace("\\", "/"),
        "raw_vertices": raw_vertex_count,
        "raw_faces": raw_face_count,
        "welded_vertices": welded_vertices,
        "pre_decimate_vertices": raw_vertex_count - welded_vertices,
        "target_faces": target,
        "decimate_ratio": ratio,
        "reduced_vertices": reduced_vertices,
        "reduced_faces": reduced_faces,
        "target_met": reduced_faces <= max(target, int(target * 1.03)),
        "raw_bytes": raw_path.stat().st_size,
        "reduced_bytes": out_path.stat().st_size,
    }


def main() -> None:
    REDUCED_DIR.mkdir(parents=True, exist_ok=True)
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete()

    requested = selected_ids()
    existing_entries = []
    existing_skipped = []
    if requested is not None and REPORT_PATH.exists():
        existing = json.loads(REPORT_PATH.read_text(encoding="utf-8"))
        existing_entries = [
            entry
            for entry in existing.get("entries", [])
            if entry.get("pdb_id") not in requested
            and (ROOT / entry.get("reduced_obj", "")).is_file()
        ]
        existing_skipped = [
            entry
            for entry in existing.get("skipped", [])
            if entry.get("pdb_id") not in requested and (ROOT / entry.get("raw_obj", "")).is_file()
        ]

    entries = []
    skipped = []
    for pdb_id, component, path in iter_raw_assets():
        try:
            entries.append(reduce_one(pdb_id, component, path))
        except Exception as exc:  # noqa: BLE001 - report and continue reducing independent assets.
            skipped.append(
                {
                    "pdb_id": pdb_id,
                    "component": component,
                    "raw_obj": str(path.relative_to(ROOT)).replace("\\", "/"),
                    "error": repr(exc),
                }
            )

    report = {
        "title": "Reduced PyMOL surface assets",
        "entries": existing_entries + entries,
        "skipped": existing_skipped + skipped,
        "targets": {
            "ribosome": {
                "1J5E_1JJ2_protein": 80_000,
                "1J5E_1JJ2_nucleic": 120_000,
                "ribosomal_protein": 80_000,
                "rRNA": 120_000,
                "tRNA_each": 25_000,
                "mRNA": 8_000,
            },
            "generic": GENERIC_TARGETS,
            "procedural_nucleic_acids": "ignored by canonical workflow unless SURFACE_REDUCTION_INCLUDE_PROCEDURAL=1",
        },
    }
    REPORT_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
