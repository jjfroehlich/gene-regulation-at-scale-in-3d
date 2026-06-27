param(
    [string]$BlenderExe = "C:\Program Files\Blender Foundation\Blender 5.1\blender.exe"
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$env:GENE_SCENE_ROOT = $Root

if (-not (Test-Path -LiteralPath $BlenderExe)) {
    throw "Blender executable not found: $BlenderExe"
}

Push-Location $Root
try {
    & $BlenderExe --background --python (Join-Path $Root "experiments\procedural_nucleic_acids\scripts\build_pymol_proxy_polish_experiment.py")
    if ($LASTEXITCODE -ne 0) {
        throw "PyMOL proxy polish comparison failed with exit code $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
