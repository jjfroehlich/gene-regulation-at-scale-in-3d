# RNA Surface and Compact mRNP Variants

This is one unified experiment with two comparison areas. The first compares three ssRNA surface treatments; the second compares six compact-mRNA models. Every model represents the same `1,852 nt` ACTB mRNA contour (`222.24 mm` at scene scale). It does not modify canonical V5.

## Scene-RNA surface study

The three scene-RNA candidates share the exact same centerline and contain no explicit stem loops or base-pair rungs. Holding the shape fixed isolates the surface treatment:

- `surface_smooth_tube`: clean rounded reference surface.
- `surface_soft_molecular`: shallow, low-frequency molecular irregularity without protruding lobes.
- `surface_twisted_groove`: a subtle rotating oval groove that gives the polymer directionality without a barbed silhouette.

## Compact mRNP architecture study

The compact area includes the three original RNA-only baselines and three protein-containing architectural analogies. All retain schematic paired RNA stems; none is a sequence-resolved ACTB fold or literal protein composition:

- `compact_rosette`: RNA-only globular baseline.
- `compact_crescent`: RNA-only crescent baseline.
- `compact_multi_domain`: densely paired RNA-only stress test.

- `rnp_ejc_clamped`: several localized protein clamps distributed around a compact RNA rosette, inspired by RNA engagement in the exon-junction complex.
- `rnp_srp_scaffold`: an RNA crescent hugged by elongated protein saddles, inspired by the combined RNA-protein recognition surface of SRP.
- `rnp_telomerase_bilobal`: two asymmetric protein-rich lobes bridged and wrapped by RNA, inspired by telomerase RNP architecture.

Of these six, `rnp_ejc_clamped` is the closest overall visual analogy to a cellular mRNP: a flexible compact RNA with multiple localized protein contacts. The EJC name describes the structural inspiration only; it does not assert that mature ACTB is permanently coated by EJCs. The RNA-only candidates omit essential protein mass, while the SRP- and telomerase-inspired candidates are closer to specialized, stable RNP machines than to a typical messenger RNP.

Every candidate preserves `84 + 1125 + 643 = 1852 nt` and the exact centerline contour length. Compact base-pair rungs are subsampled only for visual legibility.

## Run

```powershell
.\experiments\rna_structure_variants\run_rna_structure_variants.ps1
```

Outputs are written to `experiments/rna_structure_variants/outputs/`, including two comparison images, nine detail renders, one unified Blender scene, and a machine-readable report.

## Structural inspiration

- [Exon-junction complex, PDB 2J0Q](https://www.rcsb.org/structure/2j0q): a local multi-protein clamp around a short RNA segment.
- [ALYREF-EJC mRNP complex, PDB 7ZNJ](https://www.rcsb.org/structure/7ZNJ): multivalent recognition of an EJC-decorated mature mRNP.
- [Signal-recognition particle core, PDB 1DUL](https://www.rcsb.org/structure/1DUL): RNA and protein create a combined recognition surface.
- [SRP:SR complex, PDB 2XXA](https://www.rcsb.org/structure/2xxa): protein assemblies occupy distinct positions along an RNA scaffold.
- [Human telomerase holoenzyme, PDB 7BG9](https://www.rcsb.org/structure/7bg9): RNA bridges a catalytic core and a second RNP lobe.
- [Human telomerase catalytic core, PDB 7TRF](https://www.rcsb.org/structure/7TRF): an asymmetric bilobal RNP architecture.

The deposited structures are used only to identify transferable visual motifs—clamps, saddles, bridges, and lobes. The experiment does not imply that ACTB mRNA contains these exact complexes or stoichiometries.
