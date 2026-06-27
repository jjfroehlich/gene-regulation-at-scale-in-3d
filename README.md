# Scale-Accurate Gene Expression Scene

This repository now has a canonical V4 workflow for the main 3D scene, preserved V3 outputs as an intermediate backup, and separate folders for experiments.

## Canonical Workflow

Run the complete V4 build from a PowerShell prompt:

```powershell
.\scripts\run_canonical_v4_workflow.ps1
```

For a quick rebuild that reuses already downloaded/exported/reduced assets:

```powershell
.\scripts\run_canonical_v4_workflow.ps1 -SkipFetch -SkipPyMolExport -SkipReduction
```

Canonical V4 outputs are written to:

- `outputs/canonical/gene_expression_surface_style_v4.blend`
- `outputs/canonical/preview_gene_expression_surface_style_v4.png`
- `outputs/canonical/gene_expression_surface_scene_v4_report.json`

V4 detail previews are also written beside the main preview:

- `outputs/canonical/preview_gene_expression_surface_style_v4_full_overview.png`
- `outputs/canonical/preview_gene_expression_surface_style_v4_polymerase_rna_start.png`
- `outputs/canonical/preview_gene_expression_surface_style_v4_nucleosome_loop.png`
- `outputs/canonical/preview_gene_expression_surface_style_v4_mrna_spiral.png`
- `outputs/canonical/preview_gene_expression_surface_style_v4_ribosome_top_translation.png`
- `outputs/canonical/preview_gene_expression_surface_style_v4_compact_mrna.png`

Preserved V3 backup outputs remain at:

- `outputs/canonical/gene_expression_surface_style_v3.blend`
- `outputs/canonical/preview_gene_expression_surface_style_v3.png`
- `outputs/canonical/gene_expression_surface_scene_v3_report.json`

Detail previews are also written beside the main preview:

- `outputs/canonical/preview_gene_expression_surface_style_v3_dna_transcription_detail.png`
- `outputs/canonical/preview_gene_expression_surface_style_v3_cas9_binding_detail.png`
- `outputs/canonical/preview_gene_expression_surface_style_v3_nucleosome_binding_detail.png`
- `outputs/canonical/preview_gene_expression_surface_style_v3_translation_detail.png`
- `outputs/canonical/preview_gene_expression_surface_style_v3_ribosome_detail.png`

The canonical V4 build uses:

- `config/scene_manifest_v4.json` for the versioned V4 scene manifest, derived from `config/scene_manifest.json`.
- `scripts/fetch_rcsb_assets.py` to fetch mmCIF files.
- `scripts/export_pymol_surface_assets.py` to export PyMOL molecular surfaces for PDB-derived structures.
- `scripts/reduce_surface_assets.py` to weld duplicate PyMOL OBJ vertices and decimate PDB-derived surfaces for Blender.
- `scripts/blender_nucleic_meshes.py` to build scale-correct direct Blender DNA/RNA meshes with nucleotide-level detail.
- `scripts/build_gene_expression_surface_scene_v4.py` to arrange the final V4 Blender scene.

## Scene Style

The main scene uses PyMOL molecular surfaces for imported PDB-derived structures. DNA and mRNA are generated directly in Blender as scale-correct polished meshes, avoiding the procedural DNA/RNA PyMOL export path in the main scene builds. V4 uses the arrangement V2 full-gene layout: 3,454 bp ACTB canonical gene span plus a 500 bp upstream promoter, with mRNA rising from the RNA polymerase II/TSS point. A separate compact mRNP-like full-length actin mRNA remains as a secondary reference. DNA- and RNA-bound proteins are registered using their co-crystallized nucleic-acid geometry with strict guide-to-path contact validation.

The DNA/RNA-only arrangement comparison is in `experiments/arrangement_variants`. It builds four full-scale direct Blender DNA/RNA panels in one `.blend` file so the base and upward RNA spiral can be chosen before protein attachment is revisited. Canonical DNA length notes are recorded in `docs/canonical_dna_scale_notes.md`.

Current canonical PDB assets include RNA polymerase II elongation complex `2E2I`, ribosome subunits `1J5E` and `1JJ2`, tRNA `4TNA`, nucleosome `1AOI`, Cas9 `4UN3`, Argonaute `1U04`, and p53 tetramer bound to DNA `3TS8`.

## Experiments

Experiments are isolated under `experiments/`:

- `experiments/style_trials/`: earlier surface/spheres/beads comparisons.
- `experiments/molecular_nodes/`: earlier Molecular Nodes import/style trial.
- `experiments/archive_v1/`: older first-pass scene outputs.
- `experiments/procedural_nucleic_acids/`: current custom-vs-PyMOL-calibrator-vs-Molecular-Nodes DNA/RNA comparison.

The current nucleic-acid comparison can be refreshed with:

```powershell
python .\scripts\analyze_nucleic_calibrators.py
$env:GENE_SCENE_ROOT=(Get-Location).Path; $env:PYMOL_CALIBRATOR_PDBS='1BNA,9IOB'; $env:PYMOL_CALIBRATOR_REPRESENTATIVE_CHAINS='1'; & "$env:LOCALAPPDATA\Schrodinger\PyMOL2\Scripts\pymol.exe" -cq -d "run scripts/export_pymol_nucleic_calibrator_surfaces.py"
& "C:\Program Files\Blender Foundation\Blender 5.1\blender.exe" --background --python .\experiments\procedural_nucleic_acids\scripts\build_nucleic_acid_generation_comparison.py
```

Generated canonical raw/reduced PyMOL surface assets remain in `assets/pymol_exports/` because they are shared by the main workflow.
