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

& $Blender --background --python scripts\export_arrangement_v1_source_paths.py

if ($SkipPymolExport) {
    & (Join-Path $PSScriptRoot "run_arrangement_v1_workflow.ps1") -SkipPymolExport
} else {
    & (Join-Path $PSScriptRoot "run_arrangement_v1_workflow.ps1")
}
