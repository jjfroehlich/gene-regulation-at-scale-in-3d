#!/usr/bin/env python3
"""Export first-pass PyMOL geometry/style trials for Blender import.

Run inside PyMOL:
    pymol -cq -d "run C:/path/to/experiments/style_trials/scripts/export_pymol_style_trials.py"
"""

from __future__ import annotations

import json
import os
from pathlib import Path

from pymol import cmd


ROOT = Path(os.environ.get("GENE_SCENE_ROOT", Path(__file__).resolve().parents[3])).resolve()
RCSB_DIR = ROOT / "assets" / "rcsb"
OUTPUT_DIR = ROOT / "experiments" / "style_trials" / "assets" / "pymol_exports" / "style_trials"
ANGSTROM_TO_MM = 0.04


PDBS = {
    "1NKP": {
        "label": "TF-DNA complex",
        "cif": RCSB_DIR / "1NKP.cif",
        "protein_color": "cyan",
        "nucleic_color": "orange",
    },
    "1KX5": {
        "label": "Nucleosome",
        "cif": RCSB_DIR / "1KX5.cif",
        "protein_color": "forest",
        "nucleic_color": "orange",
    },
}

STYLES = {
    "surface": {
        "label": "PyMOL molecular surface",
        "representation": "surface",
    },
    "spheres": {
        "label": "PyMOL low-quality atom spheres",
        "representation": "spheres",
    },
}

COMPONENTS = {
    "protein": "polymer.protein",
    "nucleic": "polymer.nucleic",
}


def selected(mapping: dict, env_name: str) -> dict:
    raw = os.environ.get(env_name, "").strip()
    if not raw:
        return mapping
    wanted = {item.strip() for item in raw.split(",") if item.strip()}
    return {key: value for key, value in mapping.items() if key in wanted}


def safe_path(path: Path) -> str:
    return path.as_posix()


def configure_scene() -> None:
    cmd.bg_color("white")
    cmd.set("orthoscopic", 1)
    cmd.set("ray_opaque_background", 0)
    cmd.set("specular", 0)
    cmd.set("ambient", 0.72)
    cmd.set("direct", 0.32)
    cmd.set("shininess", 0)
    cmd.set("two_sided_lighting", 1)
    cmd.set("depth_cue", 0)


def apply_style(style: str, object_name: str) -> None:
    cmd.hide("everything", object_name)
    if style == "surface":
        cmd.set("surface_quality", 0)
        cmd.set("surface_smooth_edges", 1)
        cmd.set("solvent_radius", 1.4)
        cmd.show("surface", f"{object_name} and polymer")
    elif style == "spheres":
        cmd.set("sphere_quality", 0)
        cmd.set("sphere_scale", 0.55)
        cmd.show("spheres", f"{object_name} and polymer")
    else:
        raise ValueError(style)


def show_style_for_selection(style: str, selection: str) -> None:
    if style == "surface":
        cmd.show("surface", selection)
    elif style == "spheres":
        cmd.show("spheres", selection)
    else:
        raise ValueError(style)


def export_component(object_name: str, style: str, component: str, output_path: Path) -> int:
    selection = f"{object_name}_{component}_sel"
    cmd.select(selection, f"{object_name} and ({COMPONENTS[component]})")
    count = cmd.count_atoms(selection)
    if count:
        cmd.hide("everything", object_name)
        show_style_for_selection(style, selection)
        cmd.save(safe_path(output_path), selection, state=1)
    cmd.delete(selection)
    return count


def render_preview(object_name: str, output_path: Path) -> None:
    cmd.orient(object_name)
    cmd.zoom(object_name, 7)
    # PyMOL ray rendering can take a long time for nucleosomes and is not needed
    # for this pass; Blender will render the style contact sheet after import.
    cmd.png(safe_path(output_path), width=1000, height=700, dpi=120, ray=0)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pdbs = selected(PDBS, "STYLE_TRIAL_PDBS")
    styles = selected(STYLES, "STYLE_TRIAL_STYLES")
    pdb_manifest = {
        key: {inner_key: str(inner_value) for inner_key, inner_value in value.items()}
        for key, value in pdbs.items()
    }
    manifest = {
        "title": "PyMOL style trial exports",
        "angstrom_to_mm": ANGSTROM_TO_MM,
        "pdbs": pdb_manifest,
        "styles": styles,
        "exports": [],
    }

    for pdb_id, pdb_meta in pdbs.items():
        cif_path = pdb_meta["cif"]
        if not cif_path.exists():
            raise FileNotFoundError(cif_path)
        for style_id, style_meta in styles.items():
            object_name = f"{pdb_id}_{style_id}"
            cmd.reinitialize()
            configure_scene()
            cmd.load(safe_path(cif_path), object_name)
            cmd.remove(f"{object_name} and solvent")
            cmd.remove(f"{object_name} and hydro")
            apply_style(style_id, object_name)

            cmd.color(pdb_meta["protein_color"], f"{object_name} and polymer.protein")
            cmd.color(pdb_meta["nucleic_color"], f"{object_name} and polymer.nucleic")

            out_dir = OUTPUT_DIR / pdb_id / style_id
            out_dir.mkdir(parents=True, exist_ok=True)
            png_path = out_dir / f"{pdb_id}_{style_id}_preview.png"

            for component in COMPONENTS:
                obj_path = out_dir / f"{pdb_id}_{style_id}_{component}.obj"
                atom_count = export_component(object_name, style_id, component, obj_path)
                manifest["exports"].append(
                    {
                        "pdb_id": pdb_id,
                        "structure_label": pdb_meta["label"],
                        "style": style_id,
                        "style_label": style_meta["label"],
                        "component": component,
                        "atom_count": atom_count,
                        "obj": str(obj_path.relative_to(ROOT)).replace("\\", "/"),
                        "preview": str(png_path.relative_to(ROOT)).replace("\\", "/"),
                        "exists": obj_path.exists(),
                        "bytes": obj_path.stat().st_size if obj_path.exists() else 0,
                    }
                )

    manifest_path = OUTPUT_DIR / "style_trial_exports_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote {manifest_path}")
    cmd.quit()


main()
