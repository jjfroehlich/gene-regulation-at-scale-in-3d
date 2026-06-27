# Arrangement Variants

This experiment compares DNA/RNA-only arrangements before protein placement is revisited.

Run:

```powershell
.\experiments\arrangement_variants\run_arrangement_variants_workflow.ps1
```

Outputs are written to `experiments/arrangement_variants/outputs/`.

The generated Blender file contains eight full-scale DNA/RNA arrangement panels:

- compact irregular 3/4 DNA base with upward RNA spiral
- compact DNA spiral pedestal using 3,454 bp ACTB plus 500 bp upstream promoter DNA
- nested irregular DNA coil with wigglier RNA
- broader organic DNA base with tighter vertical RNA spiral
- full-gene serpentine back-and-forth DNA base with rounded turns
- full-gene irregular ribbon DNA base with broader rounded folds
- full-gene compact loop DNA base
- full-gene varied serpentine base with unequal lane lengths and rounded turns

The source path curves are included as render-hidden editable curves for rebaking a chosen variant.
