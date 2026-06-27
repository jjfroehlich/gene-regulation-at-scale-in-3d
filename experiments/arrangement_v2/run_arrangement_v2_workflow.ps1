param(
    [switch]$SkipPymolExport
)

$ErrorActionPreference = "Stop"

$ProjectRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..")
Set-Location $ProjectRoot

$Blender = "C:\Program Files\Blender Foundation\Blender 5.1\blender.exe"

if (-not (Test-Path -LiteralPath $Blender)) {
    throw "Blender was not found at $Blender"
}

$env:GENE_SCENE_MANIFEST = "experiments\arrangement_v2\arrangement_v2_manifest.json"

python scripts\write_arrangement_v2_manifest.py
& $Blender --background --python scripts\build_arrangement_v2_scene.py

Write-Host "Wrote experiments\arrangement_v2\outputs\gene_expression_arrangement_v2.blend"
Write-Host "Wrote experiments\arrangement_v2\outputs\preview_gene_expression_arrangement_v2.png"
