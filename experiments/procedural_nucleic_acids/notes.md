# Procedural Nucleic Acid Notes

The canonical scene uses direct Blender procedural DNA/RNA meshes. PDB-derived structures still go through PyMOL first.

Current implementation:

- DNA is a B-form double helix around a flexible centerline.
- The DNA centerline can be edited in `scripts/procedural_nucleic_geometry.py`.
- mRNA is a full-length segmented ssRNA path with segment-aware coloring.
- The mRNA path is the binding/placement reference for RBPs and the ribosome.
- Visible procedural DNA/RNA meshes are built by `scripts/blender_nucleic_meshes.py`, not by PyMOL.

Current comparison outputs:

- `outputs/nucleic_acid_calibrator_analysis.json`
- `outputs/nucleic_acid_generation_comparison.blend`
- `outputs/preview_nucleic_acid_generation_comparison.png`
- `outputs/nucleic_acid_generation_comparison_report.json`

Molecular Nodes remains useful for importing and styling real structures. Its oxDNA/oxRNA path is useful as a control because it is path-like and scale-close, but it is coarse-grained rather than surface-style, so it is not used by the canonical workflow.

PyMOL proxy polish experiment:

- Rebuild with `run_pymol_proxy_polish_comparison.ps1`.
- Output blend: `outputs/pymol_proxy_polish_comparison.blend`.
- The comparison tests repaired PyMOL proxy surfaces, per-component voxel remesh, unified voxel remesh, and unified voxel remesh with nearest-component material transfer.
- Voxel remesh is the strongest Blender-side repair for small holes and ragged PyMOL mesh seams. Per-component remesh preserves colors but does not fuse inter-component transitions; unified remesh fuses the surface, and material transfer recovers approximate component colors afterward.
- DNA material transfer uses a strand-priority cleanup so base-pair color does not bleed across the outside of the fused helix. RNA keeps nearest-component transfer because the segment colors remain clean.

DNA/RNA optimization V3:

- Rebuild with `run_dna_rna_optimization_v3.ps1`.
- Output blend: `outputs/dna_rna_optimization_v3.blend`.
- The experiment keeps the canonical scene unchanged and compares current canonical DNA/RNA with optimized right-handed B-DNA, an irregular elongated full-length mRNA, a compact schematic full-length mRNA, and real PDB calibrator surfaces from `1BNA` and `9IOB`.
- V3 keeps the PyMOL surface proxy route: procedural pseudoatoms in Angstrom coordinates, PyMOL surface OBJ export, Blender reduction, unified voxel polish, and analytic component material reassignment.
- The DNA generator records the handedness metric used for `1BNA`; both mRNA versions preserve the 84/1125/643 nt segment lengths and 1852 nt total contour length.
