# Workflow Notes

## Main Build

The current V5 scene is built by `scripts/run_canonical_v5_workflow.ps1`.

Stages:

1. Write `config/scene_manifest_v5.json` from the base manifest plus the reader-order ACTB promoter-plus-gene layout.
2. Fetch mmCIF files from RCSB into `assets/rcsb/`.
3. Export component-level PyMOL surfaces for PDB-derived scene structures into `assets/pymol_exports/surface_assets/`.
4. Reduce PDB-derived OBJ surfaces into `assets/pymol_exports/surface_assets_reduced/`.
5. Build direct procedural DNA/RNA meshes in Blender and render the scene into `outputs/canonical/`.

Current V5 outputs:

- `outputs/canonical/gene_expression_surface_style_v5.blend`
- `outputs/canonical/preview_gene_expression_surface_style_v5.png`
- `outputs/canonical/gene_expression_surface_scene_v5_report.json`

Detail renders are also written under `outputs/canonical/` for the overview, polymerase/RNA start, nucleosome loop, mRNA spiral, ribosome/translation view, and compact mRNA reference.

The base manifest lives at `config/scene_manifest.json`. The V5 manifest is `config/scene_manifest_v5.json`; it is generated and should be treated as the source for V5 output reports.

## Legacy Builds

The previous V4 build can still be rebuilt with `scripts/run_canonical_v4_workflow.ps1`. Its outputs use the `gene_expression_surface_style_v4.*` and `gene_expression_surface_scene_v4_report.json` names.

The preserved V3 build can still be rebuilt by `scripts/run_canonical_workflow.ps1`. Its outputs use the `gene_expression_surface_style_v3.*` and `gene_expression_surface_scene_v3_report.json` names.

## Nucleic Acids

The main scene uses PyMOL surfaces for structures that came from PDB. Procedural DNA/RNA are generated directly in Blender:

- DNA: B-form helix along the V5 reader-order ACTB promoter-plus-gene serpentine path, `3,954 bp`, `0.34 nm/bp`, `10.5 bp/turn`, with a `2.0 nm` envelope target.
- mRNA: flexible full-length actin mRNA path, `1,852 nt * 0.30 nm/nt`, split into `84/1125/643 nt` UTR/coding/UTR segments.
- Surface proxy: strand/backbone/base-lobe meshes are built directly in millimeter coordinates by `scripts/blender_nucleic_meshes.py`.
- PDB-derived assets with DNA/RNA are path-aligned so transcription factors and RBPs visibly contact the procedural nucleic acid path.

The V5 build fails if a PDB-derived structure would fall back to beads, if scale validation fails, or if a strict-contact protein is not attached to the expected DNA/RNA path. Reports mark imported PDB components as `pymol_surface_reduced` and procedural DNA/RNA as `direct_blender_surface_proxy_polished`.

## Nucleic-Acid Comparisons

The current comparison workflow is under `experiments/procedural_nucleic_acids/`:

- `scripts/analyze_nucleic_calibrators.py` analyzes all configured DNA/RNA PDB calibrators.
- `scripts/export_pymol_nucleic_calibrator_surfaces.py` exports protein-stripped PyMOL nucleic-acid surfaces for representative calibrators.
- `experiments/procedural_nucleic_acids/scripts/build_nucleic_acid_generation_comparison.py` builds the visual comparison of custom direct meshes, PyMOL calibrators, Molecular Nodes Style Surface imports, and Molecular Nodes oxDNA/oxRNA.

Molecular Nodes remains experimental. In the current report, oxDNA/oxRNA are scale-close and path-controllable, but coarse-grained rather than surface-style. Molecular Nodes Style Surface is useful for real structures, not for the path-controlled full DNA/mRNA assets used in the main scene.
