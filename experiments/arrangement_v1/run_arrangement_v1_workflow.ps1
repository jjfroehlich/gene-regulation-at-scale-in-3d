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

$env:GENE_SCENE_MANIFEST = "experiments\arrangement_v1\arrangement_v1_manifest.json"

python scripts\write_arrangement_v1_manifest.py
& $Blender --background --python scripts\build_arrangement_v1_scene.py

Write-Host "Wrote experiments\arrangement_v1\outputs\gene_expression_arrangement_v1.blend"
Write-Host "Wrote experiments\arrangement_v1\outputs\preview_gene_expression_arrangement_v1.png"
