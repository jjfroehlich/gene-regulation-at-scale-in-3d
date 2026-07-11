# V5 Flythrough Animation

This experiment adds a short camera flight on top of the existing canonical V5 `.blend`. It does not change the canonical V5 scene builder.

The current cut uses a perspective camera, keyed focal length, keyed depth of field, longer camera holds on important molecules, animated molecule color pulses, localized beat lights, soft key/rim lighting, subtle volumetric atmosphere, compositor glow/color grading, and camera captions so the wide scale shots and molecule close-ups read more like a guided scientific animation.

## Run

Build the animation setup, render the MP4, and refresh the README GIF preview:

```powershell
.\experiments\v5_flythrough_animation\run_v5_flythrough_animation.ps1
```

Run a fast low-resolution smoke test:

```powershell
.\experiments\v5_flythrough_animation\run_v5_flythrough_animation.ps1 -SmokeTest
```

Build the animation `.blend` and report without rendering video:

```powershell
.\experiments\v5_flythrough_animation\run_v5_flythrough_animation.ps1 -SkipVideoRender
```

The default export is 60 seconds at 24 fps. Use `-DurationSeconds 75` for a slower version.

## Outputs

Generated files are written under `experiments/v5_flythrough_animation/outputs/`, which is ignored by Git:

- `v5_flythrough_animation.blend`
- `v5_flythrough_animation_1080p.mp4`
- `v5_flythrough_animation_smoke.mp4`
- `v5_flythrough_animation_report.json`

The tracked README preview is:

- `docs/images/v5-flythrough-preview.gif`

Full renders sample representative frames for the README GIF so the tracked preview stays compact. Smoke-test renders use the richer short preview directly.

## Storyboard

- `0-10s`: high top-down scale view with a longer hold on DNA, mRNA, ribosome, compact mRNP-like RNA, and actin.
- `10-18s`: zoom into ACTB DNA and the 3,954 bp / 537.7 mm callout.
- `18-25s`: close-up hold on the p53 tetramer bound to DNA, with green surface pulse and localized light.
- `25-35s`: RNA polymerase II hold and nearby nucleosome context, with cyan/violet surface pulses.
- `35-48s`: zoom out to the full 1,852 nt / 222.2 mm actin mRNA, with an orange full-path pulse and hold.
- `48-54s`: RNA-binding proteins, ribosome, and tRNA translation region, with warm ribosome/tRNA pulse lighting.
- `54-60s`: actin protein endpoint and compact structured-RNA reference, with a longer red actin pulse hold.

## Visual Notes

- Wide shots use 20-24 mm lenses and higher f-stop so the whole scale model remains legible.
- TF, RNA Pol II, ribosome, and actin close-ups use 80-90 mm lenses and lower f-stop for stronger depth of field.
- Important positions use duplicate camera keyframes so the camera eases into the view, pauses briefly, then continues.
- Molecule emphasis uses animation-only overlay meshes and pulse lights; the canonical V5 materials are not modified.
- Earlier target-ring callouts were removed because the turquoise wireframe look distracted from the molecular surfaces.
- The experiment copy hides the original static labels and adds animation-specific world labels plus camera captions.
- The scene uses a warm area/sun key, cool fill, rim area/point lights, low-strength camera eye light, ambient occlusion, subtle volumetric fog, restrained bloom, compositor fog glow, and per-shot beat lights for readable molecular surfaces.
- The compositor applies a restrained grade with cool shadows, warm highlights, and gentle saturation so labels and pulses pop without washing out the DNA/RNA scale paths.
- Color coding is stable across the cut: DNA blue, RNA orange, transcription factor green, RNA Pol II cyan, nucleosome violet, ribosome gold, tRNA green, and actin red.

## Useful Blender References

- [Blender Manual: Cameras](https://docs.blender.org/manual/en/latest/render/cameras.html)
- [Blender Manual: Lights](https://docs.blender.org/manual/en/latest/render/lights/index.html)
- [Blender Manual: Color Management](https://docs.blender.org/manual/en/latest/render/color_management.html)
- [Blender Studio Training](https://studio.blender.org/training/)

The repeatable export path is the headless Blender script invoked by the runner.
