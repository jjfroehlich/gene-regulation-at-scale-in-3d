# Storage Cleanup 2026-06-19

Canonical V5 is the current scene. Canonical V4 and the retained V3 outputs remain as backups, while V2 and old experiment-heavy artifacts are no longer needed for the active build path.

The useful arrangement work from `experiments/arrangement_v2` was ported into Canonical V4/V5. V5 further adds reader-order DNA, hidden source curves, tightened protein/RNA/DNA contact validation, and compact mRNP placement near the translation product.

Direct Blender DNA/RNA generation superseded the older procedural/PyMOL proxy mesh experiments. The old experiment reports and previews document the conclusions; the large OBJ and Blend artifacts can be regenerated if those experiments need to be revisited.

Removed categories:

- Raw PyMOL surface exports from `assets/pymol_exports/surface_assets/`.
- PyMOL scratch exports from `assets/pymol_exports/scratch/`.
- Reduced surface directories and mmCIF files not referenced by `config/scene_manifest_v5.json`.
- Heavy old experiment assets and `.blend`/`.blend1` outputs while retaining scripts, JSON reports, and previews.
- Canonical `.blend1` backups and V2 canonical blend files.

Regeneration notes:

- RCSB mmCIF files can be refetched with `scripts/fetch_rcsb_assets.py`.
- Raw PyMOL surfaces can be regenerated with `scripts/export_pymol_surface_assets.py`.
- Reduced surfaces can be regenerated with `scripts/reduce_surface_assets.py`.
- The current V5 scene can be rebuilt with:

```powershell
.\scripts\run_canonical_v5_workflow.ps1 -SkipFetch -SkipPyMolExport -SkipReduction
```

The cleanup preserved the V5-used PDB IDs:

```text
1AOI 1CVJ 1J5E 1J6Z 1JJ2 1U04 1ZDH 2E2I 2H5Q 3G73 3Q0Q 3TS8 4ED5 4TNA 4UN3 6KKS 6ML2
```
