$ErrorActionPreference = "Stop"

$Root = Resolve-Path (Join-Path $PSScriptRoot "..\..")
$env:GENE_SCENE_ROOT = $Root.Path

$BlenderCandidates = @(
  "C:\Program Files\Blender Foundation\Blender 5.1\blender.exe",
  "C:\Program Files\Blender Foundation\Blender 5.0\blender.exe",
  "C:\Program Files\Blender Foundation\Blender\blender.exe"
)

$Blender = $BlenderCandidates | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
if (-not $Blender) {
  throw "Could not find Blender. Checked: $($BlenderCandidates -join ', ')"
}

& $Blender --background --python (Join-Path $PSScriptRoot "scripts\build_molecular_nodes_nucleic_comparison_v2.py")
