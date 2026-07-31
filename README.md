# Gene Regulation at Scale in 3D

A three-dimensional sculpture of transcriptional and post-transcriptional control of gene expression. The scene follows a gene DNA path, mRNA spiral, to the final protein product for actin (ACTB). It includes RNA- and DNA-bound regulatory proteins, polymerase ii, Cas9, and the ribosomal translation machinery. And it is scale-accurate!

![Canonical overview render](docs/images/overview.jpg)

## Render Details

|  |  |
| --- | --- |
| ![Full p53 tetramer bound to DNA](docs/images/p53-dna.jpg) | ![Nucleosome core with wrapped DNA](docs/images/nucleosome-loop.jpg) |

|  |  |
| --- | --- |
| ![RNA polymerase II at the gene end with the nascent RNA 3′ end](docs/images/transcription-end.jpg) | ![Ribosome with tRNA](docs/images/translation.jpg) |

|  |  |
| --- | --- |
| ![ACTB protein product](docs/images/actin.jpg) | ![Cas9 with guide and target DNA](docs/images/cas9-dna.jpg) |

## What Is Being Built Here
- Scale-accurate 3D scene. 
- Full ACTB promoter-plus-gene DNA path: `3,454 bp` canonical ACTB span plus `500 bp` upstream promoter, for `3,954 bp` total.
- Full-length actin mRNA: `1,852 nt`, split into `5' UTR`, coding sequence, and `3' UTR` segments.
- RNA polymerase II sits at the completed gene endpoint, where the full transcript's nascent `3′` end remains attached; the downstream RNA-to-protein branch is kept as one coherent arrangement.
- Blender-generated DNA/RNA proxies at one shared physical scale.
- Reduced surfaces for PyMOL PDB-derived proteins and nucleoprotein complexes.
- Current PDB-derived scene assets include RNA polymerase II elongation complex `2E2I`, ribosome subunits `1J5E` and `1JJ2`, tRNA `4TNA`, nucleosome `1AOI`, Cas9 `4UN3`, transcription factors `6ML2`, `6KKS`, `3TS8`, and `3G73`, RNA-binding proteins `1U04`, `1CVJ`, `4ED5`, `3Q0Q`, and `1ZDH`, mCherry/RFP tag `2H5Q`, and actin protein `1J6Z`.
- PDB-derived molecular surfaces with generated DNA/RNA meshes at shared scale: `1 nm = 0.4 mm`.
- Note: this is of course not a realistic situation from a crowded cell. Also, between transcription and before translation the mRNA would need to be exported from the nucleus. 

## Flythrough

<p align="center">
  <img src="docs/images/flythrough-preview.gif" alt="Educational flythrough animation preview">
</p>

## 2017 Collage

This project builds on a collage I made in 2017.

![2017 gene expression scale poster](docs/images/poster-2017-gene-expression-scale.png)


> Transcriptional- and posttranscriptional control of gene expression in scale. All structures from David Goodsell & RCSB PDB. Main parts from https://mm.rcsb.org. The mRNA and the idea of the scale of actin mRNA to actin protein from https://book.bionumbers.org/which-is-bigger-mrna-or-the-protein-it-codes-for/ and I took other molecules from "molecules of the month" https://pdb101.rcsb.org/motm/181 https://pdb101.rcsb.org/motm/98 https://pdb101.rcsb.org/motm/112 https://pdb101.rcsb.org/motm/31.

See [docs/references.md](docs/references.md) for source links, PDB IDs, and attribution details.


## Build from Scripts

Run the complete canonical build from a PowerShell prompt:

```powershell
.\scripts\run_canonical_workflow.ps1
```

For a quick rebuild that reuses already downloaded/exported/reduced assets:

```powershell
.\scripts\run_canonical_workflow.ps1 -SkipFetch -SkipPyMolExport -SkipReduction
```

Canonical outputs are written to:

- `outputs/canonical/gene_expression_surface_style.blend`
- `outputs/canonical/preview_gene_expression_surface_style.png`
- `outputs/canonical/gene_expression_surface_scene_report.json`

Detail previews are written beside the main preview:

- `outputs/canonical/preview_gene_expression_surface_style_full_overview.png`
- `outputs/canonical/preview_gene_expression_surface_style_p53_dna.png`
- `outputs/canonical/preview_gene_expression_surface_style_polymerase_gene_end.png`
- `outputs/canonical/preview_gene_expression_surface_style_nucleosome_loop.png`
- `outputs/canonical/preview_gene_expression_surface_style_ribosome_trna.png`
- `outputs/canonical/preview_gene_expression_surface_style_actin_product.png`
- `outputs/canonical/preview_gene_expression_surface_style_cas9_dna.png`

The canonical build uses:

- `config/scene_manifest.json` as the single resolved scene manifest.
- `scripts/fetch_rcsb_assets.py` to fetch mmCIF files.
- `scripts/export_pymol_surface_assets.py` to export PyMOL molecular surfaces.
- `scripts/reduce_surface_assets.py` to weld and decimate OBJ surfaces for Blender.
- `scripts/blender_nucleic_meshes.py` to build scale-correct direct Blender DNA/RNA meshes.
- `scripts/build_gene_expression_surface_scene.py` to build the gene-end arrangement and render the final scene.

More details in [docs/workflow.md](docs/workflow.md).

## Sketchfab Exports

Upload-ready molecular-scene exports are written to `outputs/sketchfab/`:

- `gene_expression_canonical_sketchfab.glb` (preferred)
- `gene_expression_canonical_sketchfab.fbx` (fallback)

Regenerate them from the canonical blend with:

```powershell
.\scripts\run_sketchfab_export.ps1
```

## Experiments

Experiments to try out variations or new features are isolated under `experiments/`:

- `experiments/arrangement_variants/`: DNA/RNA-only layout comparisons.
- `experiments/procedural_nucleic_acids/`: custom-vs-PyMOL-calibrator-vs-Molecular-Nodes DNA/RNA comparisons.
- `experiments/rna_structure_variants/`: scale-accurate elongated and compact RNA folding candidates with explicit stems and base pairing.
- `experiments/flythrough_animation/`: current 66-second educational camera flight and README GIF preview.
