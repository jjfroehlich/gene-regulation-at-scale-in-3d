param()

$ErrorActionPreference = "Stop"

$ProjectRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..")
Set-Location $ProjectRoot

$Blender = "C:\Program Files\Blender Foundation\Blender 5.1\blender.exe"

if (-not (Test-Path -LiteralPath $Blender)) {
    throw "Blender was not found at $Blender"
}

& $Blender --background --python scripts\build_arrangement_variants_scene.py

Write-Host "Wrote experiments\arrangement_variants\outputs\gene_expression_arrangement_variants.blend"
Write-Host "Wrote experiments\arrangement_variants\outputs\preview_gene_expression_arrangement_variants.png"
Write-Host "Wrote experiments\arrangement_variants\outputs\gene_expression_arrangement_variants_report.json"
