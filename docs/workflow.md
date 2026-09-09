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

The workflow validates `config/scene_manifest.json`, fetches required RCSB structures when requested, exports and reduces PyMOL molecular surfaces, builds direct DNA/RNA meshes, validates molecular contacts, and renders the canonical overview, six shared-scale focused views, and six individually framed close-ups. A complete build refreshes the documentation JPEGs from that same generation.

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

Each focused view also has an additional `_closeup.png` output. Shared-scale panels retain their 42 mm orthographic width; close-ups reserve space for titles and structural identifiers and display a measured 1–2–5-series nanometer scale bar. The overview uses right-side DNA, mRNA, and actin category labels with vertical extent lines, and desaturated comparison scales in the lower-left corner. All scale bars use a thicker grey stroke. The overview uses concise molecule names; structural identifiers remain in detail captions.

Final stills use AgX with Medium High Contrast, 256 samples, denoising, and an adaptive threshold of 0.015. The renderer selects an available GPU and falls back to CPU. Set the optional environment variable `CANONICAL_PREVIEW=1` for 64 samples and a 0.035 threshold.

The build also fails if Pol II is not at the final DNA coordinate, the RNA origin diverges from the gene endpoint, the lower-angle overview does not maintain the required projected DNA/mRNA separation, or any existing scale/contact/render invariant fails.

## Animation

Run `./animation/run_flythrough_animation.ps1` to build or resume the approved 66-second dark flythrough. This is an active project workflow, using the canonical molecular scene and its saved original camera route. Outputs are in `outputs/animation/`; the README GIF is generated from the delivered film at 1x speed. See [animation/README.md](../animation/README.md) for switches and validation.

## Rebuild files

The active animation retains its required original source in `outputs/animation/rollback/before_dark/` and its current resumable frames in `outputs/animation/staging/`.

## Git inputs and generated files

Track the manifest, source scripts, animation workflow and its `animation/assets/` source library, tests, documentation images, RCSB mmCIF inputs, and reduced molecular OBJ assets with their manifest. The canonical DNA/RNA geometry is built directly in Blender; the retired pseudoatom export pipeline and its proxy assets are not required.

Generated Blender scenes, MP4s, frame caches, build reports, and validation output stay under the ignored `outputs/` directory. Raw PyMOL exports are regenerable and ignored. A fresh animation run automatically creates its canonical scene and original camera-route baseline from the tracked inputs.

Use `scripts/validate_presentation.py` in background Blender for saved-scene checks; see [the animation validation commands](../animation/README.md). The tests use Python's standard-library unittest runner.
