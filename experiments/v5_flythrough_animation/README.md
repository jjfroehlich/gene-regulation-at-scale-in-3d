# V5 Flythrough Animation

This experiment adds a continuous 60-second camera flight on top of the canonical V5 `.blend`.

The current cut uses one perspective camera on a continuous Bezier path, a separately animated look-at target, keyed focal length and depth of field, brief slowdowns without intermediate stops, emission pulses on the actual molecular surfaces, soft key/rim lighting, subtle volumetric atmosphere, compositor glow/color grading, individual molecule-attached labels, and canonical-style overview callouts.

## Run

Build the animation setup, render the MP4, and refresh the README GIF preview:

```powershell
.\experiments\v5_flythrough_animation\run_v5_flythrough_animation.ps1
```

Run a fast low-resolution smoke test:

```powershell
.\experiments\v5_flythrough_animation\run_v5_flythrough_animation.ps1 -SmokeTest
```

Render the full 60-second review cut at 960×540 and 12 fps before committing to the final 1080p pass:

```powershell
.\experiments\v5_flythrough_animation\run_v5_flythrough_animation.ps1 -ReviewRender
```

Build the animation `.blend` and report without rendering video:

```powershell
.\experiments\v5_flythrough_animation\run_v5_flythrough_animation.ps1 -SkipVideoRender
```

The default export is 60 seconds at 24 fps. `-DurationSeconds` remains available to scale the complete storyboard proportionally. Render-frame directories are cleared before a fresh render and removed after successful encoding; add `-KeepFrames` when the PNG sequence is needed for diagnostics or re-encoding.

## Outputs

Generated files are written under `experiments/v5_flythrough_animation/outputs/`, which is ignored by Git:

- `v5_flythrough_animation.blend`
- `v5_flythrough_animation_review.mp4`
- `v5_flythrough_animation_report.json`
- `review_contact_sheet.png`

The tracked README preview is:

- `docs/images/v5-flythrough-preview.gif`

Full renders sample representative frames for the README GIF so the tracked preview stays compact. Smoke-test renders use the richer short preview directly.

## Storyboard

- `0-5s`: moving full-scene overview with scale bars, DNA/mRNA brackets, and a direct pointer to actin.
- `5-8s`: descend to the promoter.
- `8-21s`: follow DNA past R2R3 MYB, Cas9, ZBTB24, p53, FOXM1-DBD, RNA Pol II, and the nucleosome.
- `21-24s`: follow the transcription and mRNA-origin context.
- `24-41s`: follow mRNA past PUM2, PABP, MS2, mCherry, Argonaute, and HuR.
- `41-45s`: approach the ribosome.
- `45-49s`: move slowly across the small subunit, large subunit, and tRNA without stopping.
- `49-52s`: pass compact mRNA.
- `52-57s`: arrive at actin.
- `57-60s`: final actin hold.

## Visual Notes

- Wide shots use 20-24 mm lenses and higher f-stop so the whole scale model remains legible.
- Close-up focal lengths are eased per component, with wider traversal views used when multiple nearby structures must remain visible.
- A single Bezier path with auto-clamped animation handles provides uninterrupted camera motion. Path progress remains strictly increasing until the final actin hold.
- Molecule emphasis uses animation-local copies of the source objects' materials, so the actual protein or tRNA surface emits briefly without duplicate geometry or transform drift.
- Earlier target-ring callouts were removed because the turquoise wireframe look distracted from the molecular surfaces.
- The experiment copy hides the original static labels and adds clean, PDB-ID-free labels and leaders for all 17 named structural assets. The overview retains DNA and mRNA brackets while the protein label points directly to actin. There is no bottom caption bar and no millimeter dimensions in animation text.
- The scene uses a warm area/sun key, cool fill, rim area/point lights, low-strength camera eye light, ambient occlusion, subtle volumetric fog, restrained bloom, and compositor fog glow for readable molecular surfaces.
- The compositor applies a restrained grade with cool shadows, warm highlights, and gentle saturation so labels and pulses pop without washing out the DNA/RNA scale paths.
- Emission pulses preserve each source surface's canonical color; DNA and RNA path emphasis remains blue and orange, and actin remains coral-red.

## Useful Blender References

- [Blender Manual: Cameras](https://docs.blender.org/manual/en/latest/render/cameras.html)
- [Blender Manual: Lights](https://docs.blender.org/manual/en/latest/render/lights/index.html)
- [Blender Manual: Color Management](https://docs.blender.org/manual/en/latest/render/color_management.html)
- [Blender Studio Training](https://studio.blender.org/training/)

The repeatable export path is the headless Blender script invoked by the runner.
