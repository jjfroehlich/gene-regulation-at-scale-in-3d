# Flythrough Animation

The active delivery is a 66-second dark molecular journey at 24 fps. It follows the original camera route and look target. The original eight-second ribosome approach and pass is shortened smoothly to six seconds (48-54 seconds); the two seconds are redistributed across the earlier journey. The final four-second actin hold is preserved. Molecular geometry, physical scale and biological contacts are preserved.

The presentation uses Cycles, slate protein surfaces, the native DNA segment palette, a charcoal background, and soft key/fill/rim lighting. Gentle focus shifts and light sweeps accompany existing p53, nucleosome, polymerase and ribosome passages. The blue decorative DNA centerline is hidden so native core and helix colors agree. Existing labels and RNA highlights remain; the former atmosphere and compositor glow are disabled.

## Run

```powershell
.\animation\run_flythrough_animation.ps1
```

The default run produces a 1280 × 720 review at 64 samples followed by a 1920 × 1080 final at 256 samples, both at 24 fps. Denoising is enabled, adaptive thresholds are 0.035 and 0.015, and motion blur is disabled. GPU selection is automatic. Existing Blender, Python, FFmpeg and FFprobe installations are required; no new dependencies are introduced.

Supported switches remain available:

- `-SmokeTest`: render the storyboard checkpoints without compressing the timeline further.
- `-ReviewRender`: build the complete 720p/24 fps review in staging, without final publication.
- `-SkipVideoRender`: prepare and validate the staged scene only.
- `-SkipReadmeGif`: omit GIF regeneration.
- `-KeepFrames`: retained for compatibility; staging frames are now retained by default for recovery and inspection.
- `-DurationSeconds`, `-Fps`, `-ResolutionX`, `-ResolutionY`: explicit delivery overrides; defaults are 66, 24, 1920 and 1080.

## Staging and rollback

The pipeline reuses the original camera-route baseline in `outputs/animation/rollback/before_dark/`. On a fresh clone it builds the canonical scene from tracked molecular inputs when needed, then generates an untreated 66-second baseline before applying the dark presentation. A partial baseline is rebuilt as a pair; an existing dark delivery is never used as the original baseline. No generated scene files are required. The tracked `assets/compact_mrna.blend` is a small source library containing the four compact RNA meshes used by the approved film; `assets/baseline.json` records their checksum and framing center. This preserves the film independently of the later compact RNA revision in the still scene. Keep both source files in Git.

The active files are untouched while rendering into `outputs/animation/staging/`. Complete PNG frames resume on repeated runs with matching input and code signatures. Changed signatures archive the earlier staging directory. Do not edit a staged scene during a render.

After both videos pass decoding, frame-count, resolution and duration checks, the pipeline publishes the unversioned animation `.blend`, review and final MP4s, report, contact sheet and review page. The README GIF is regenerated at the film's actual speed, with no additional acceleration. The original camera-route source remains in rollback, and the current resumable frames remain in staging.

The original inclusive endpoint is retained: the nominal 66-second, 24 fps timeline contains 1,585 frames and encodes to approximately 66.04 seconds. The geometric camera path is unchanged; animation channels are sampled at smoothly remapped source times and baked at each output frame. The 10× optical gate preserves field of view while making focus shifts visible at the scene's coordinate scale.

## Validation and outputs

`outputs/animation/flythrough_animation_report.json` records the inherited storyboard and label data together with camera-equivalence, label-timing, geometry, DNA-overlay and rendering checks. Camera world matrices and lens/sensor ratios must match the corresponding original-time samples within 1e-6. The source canonical scene and report are hash-protected.

Full outputs retain their existing names: `flythrough_animation.blend`, `flythrough_animation_review.mp4`, and `flythrough_animation_1080p.mp4`. Open `outputs/animation/index.html` for playback. The saved scene opens at the first frame with the animation camera active.

All builds and renders run in background Blender. The canonical still-image scene and unrelated open Blender session are not modified.

To check the saved scenes without rendering or saving them, run from the project root:

```powershell
& "C:/Program Files/Blender Foundation/Blender 5.1/blender.exe" --background --factory-startup --python-exit-code 1 --python scripts/validate_presentation.py
python -m unittest discover -s tests
```

The scene validator checks the saved overview, recorded close-up scale measurements against camera projection, animation geometry, every output-frame camera and label state, path progress, focus limits, pointer projection, and the final hold. Results go to `outputs/validation/presentation.json`. Use `-- --kind still` or `-- --kind animation` for one scene, and `-- --animation-dir outputs/animation/staging` to validate a staged build. Image appearance and playback remain visual review steps.
