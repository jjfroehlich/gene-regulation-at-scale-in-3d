# Canonical Workflow

## Build

Run the current canonical scene workflow:

```powershell
.\scripts\run_canonical_workflow.ps1
```

Reuse downloaded and reduced molecular assets for a faster rebuild:

```powershell
.\scripts\run_canonical_workflow.ps1 -SkipFetch -SkipPyMolExport -SkipReduction
```

The workflow validates `config/scene_manifest.json`, fetches required RCSB structures when requested, exports and reduces PyMOL molecular surfaces, builds direct DNA/RNA meshes, validates molecular contacts, and renders the canonical overview plus six focused views.

## Scale and biological interpretation

- Shared scene scale: `1 nm = 0.4 mm`.
- ACTB promoter-plus-gene DNA: `3,954 bp`, comprising a `3,454 bp` gene span plus `500 bp` upstream promoter; contour length `537.744 mm`.
- ACTB mRNA: `1,852 nt`, comprising `84 nt` 5′ UTR, `1,125 nt` coding sequence, and `643 nt` 3′ UTR; contour length `222.24 mm`.
- Protein contour conversion: `0.36 nm` per amino acid.

The DNA uses a B-form helix along the reader-order promoter/gene path. RNA polymerase II is attached at DNA fraction `1.0`, with the full mRNA's nascent 3′ end at the same gene-end coordinate.

The elongated RNA retains the deterministic conical centerline and smooth `surface_twisted_groove` envelope without radial nucleotide lobes. The compact reference remains the RNA-only `compact_rosette`, with 38 deterministic schematic stems and a 58% paired-fraction target.

## Outputs

Canonical outputs are written under `outputs/canonical/`:

- `gene_expression_surface_style.blend`
- `gene_expression_surface_scene_report.json`
- `preview_gene_expression_surface_style.png`
- `preview_gene_expression_surface_style_full_overview.png`
- `preview_gene_expression_surface_style_p53_dna.png`
- `preview_gene_expression_surface_style_nucleosome_loop.png`
- `preview_gene_expression_surface_style_polymerase_gene_end.png`
- `preview_gene_expression_surface_style_ribosome_trna.png`
- `preview_gene_expression_surface_style_actin_product.png`
- `preview_gene_expression_surface_style_cas9_dna.png`

The build also fails if Pol II is not at the final DNA coordinate, the RNA origin diverges from the gene endpoint, the lower-angle overview does not maintain the required projected DNA/mRNA separation, or any existing scale/contact/render invariant fails.

## Retained experiments

- `experiments/arrangement_variants/` compares alternate DNA/RNA layouts.
- `experiments/procedural_nucleic_acids/` retains the nucleic-acid method comparisons and calibrators.
- `experiments/flythrough_animation/` builds the current 66-second educational camera flight and README GIF.
