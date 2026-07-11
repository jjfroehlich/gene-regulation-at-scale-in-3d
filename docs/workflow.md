# Canonical V5 Workflow

## Build

Run the current canonical scene workflow:

```powershell
.\scripts\run_canonical_v5_workflow.ps1
```

Reuse downloaded and reduced molecular assets for a faster rebuild:

```powershell
.\scripts\run_canonical_v5_workflow.ps1 -SkipFetch -SkipPyMolExport -SkipReduction
```

The workflow normalizes `config/scene_manifest_v5.json`, fetches required RCSB structures when requested, exports and reduces PyMOL molecular surfaces, builds direct DNA/RNA meshes, validates molecular contacts, and renders the canonical overview plus five focused views.

## Scale and biological interpretation

- Shared scene scale: `1 nm = 0.4 mm`.
- ACTB promoter-plus-gene DNA: `3,954 bp`, comprising a `3,454 bp` gene span plus `500 bp` upstream promoter; contour length `537.744 mm`.
- ACTB mRNA: `1,852 nt`, comprising `84 nt` 5′ UTR, `1,125 nt` coding sequence, and `643 nt` 3′ UTR; contour length `222.24 mm`.
- Protein contour conversion: `0.36 nm` per amino acid.

The DNA uses a B-form helix along the reader-order promoter/gene path. The canonical elongated RNA keeps the deterministic conical centerline and uses the smooth `surface_twisted_groove` envelope without radial nucleotide lobes. The compact reference is the RNA-only `compact_rosette`, with 38 deterministic schematic stems and a 58% paired-fraction target. Both preserve the nucleotide allocation and contour length.

## Outputs

Canonical outputs are written under `outputs/canonical/`:

- `gene_expression_surface_style_v5.blend`
- `gene_expression_surface_scene_v5_report.json`
- `preview_gene_expression_surface_style_v5.png`
- `preview_gene_expression_surface_style_v5_full_overview.png`
- `preview_gene_expression_surface_style_v5_p53_dna.png`
- `preview_gene_expression_surface_style_v5_nucleosome_loop.png`
- `preview_gene_expression_surface_style_v5_polymerase_rna_start.png`
- `preview_gene_expression_surface_style_v5_ribosome_trna.png`
- `preview_gene_expression_surface_style_v5_actin_product.png`

The build fails on scale errors, missing surface assets, non-schematic RNA metadata, wrong RNA allocation or stem counts, visible source curves, or failed strict DNA/RNA contact checks.

## Retained experiments

- `experiments/arrangement_variants/` compares alternate DNA/RNA layouts.
- `experiments/procedural_nucleic_acids/` retains the nucleic-acid method comparisons and calibrators.
- `experiments/v5_flythrough_animation/` builds the educational camera flight and README GIF.
