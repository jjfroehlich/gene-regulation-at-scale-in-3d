param()

$ErrorActionPreference = "Stop"

$ProjectRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..")
Set-Location $ProjectRoot

$Blender = "C:\Program Files\Blender Foundation\Blender 5.1\blender.exe"
if (-not (Test-Path -LiteralPath $Blender)) {
    throw "Blender was not found at $Blender"
}

& $Blender --background --python scripts\build_direct_blender_nucleic_acid_experiment.py

Write-Host "Wrote experiments\direct_blender_nucleic_acids\outputs\direct_blender_nucleic_acid_proxy_comparison.blend"
Write-Host "Wrote experiments\direct_blender_nucleic_acids\outputs\preview_direct_blender_nucleic_acid_proxy_comparison.png"
