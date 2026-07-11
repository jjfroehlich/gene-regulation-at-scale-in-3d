param()

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$Blender = "C:\Program Files\Blender Foundation\Blender 5.1\blender.exe"

if (-not (Test-Path -LiteralPath $Blender)) {
    throw "Blender executable not found: $Blender"
}

Push-Location $Root
try {
    & $Blender --background --python "experiments\rna_structure_variants\scripts\build_rna_structure_variants.py"
    if ($LASTEXITCODE -ne 0) {
        throw "RNA structure variant build failed with exit code $LASTEXITCODE"
    }
    $Report = Join-Path $Root "experiments\rna_structure_variants\outputs\rna_structure_variants_report.json"
    if (-not (Test-Path -LiteralPath $Report)) {
        throw "Blender did not write the RNA structure variant report."
    }
}
finally {
    Pop-Location
}
