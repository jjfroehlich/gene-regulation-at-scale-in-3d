param(
    [string]$PythonExe = "python",
    [string]$PyMolExe = "$env:LOCALAPPDATA\Schrodinger\PyMOL2\Scripts\pymol.exe",
    [string]$BlenderExe = "C:\Program Files\Blender Foundation\Blender 5.1\blender.exe",
    [switch]$SkipPyMolExport
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$env:GENE_SCENE_ROOT = $Root

function Invoke-Checked {
    param(
        [string]$Label,
        [scriptblock]$Command
    )
    Write-Host ""
    Write-Host "== $Label =="
    & $Command
    if ($LASTEXITCODE -ne 0) {
        throw "$Label failed with exit code $LASTEXITCODE"
    }
}

if (-not (Test-Path -LiteralPath $BlenderExe)) {
    throw "Blender executable not found: $BlenderExe"
}

if (-not $SkipPyMolExport -and -not (Test-Path -LiteralPath $PyMolExe)) {
    throw "PyMOL executable not found: $PyMolExe"
}

Push-Location $Root
try {
    Invoke-Checked "Generate detailed DNA/RNA pseudoatom CIFs" {
        & $PythonExe "experiments\procedural_nucleic_acids\scripts\generate_detail_proxy_assets.py"
    }

    if (-not $SkipPyMolExport) {
        $exportScript = (Join-Path $Root "experiments\procedural_nucleic_acids\scripts\export_detail_pymol_surfaces.py").Replace("\", "/")
        Invoke-Checked "Export detailed DNA/RNA PyMOL surfaces" {
            & $PyMolExe -cq -d "run $exportScript"
        }
    }

    Invoke-Checked "Build detailed DNA/RNA route comparison" {
        & $BlenderExe --background --python (Join-Path $Root "experiments\procedural_nucleic_acids\scripts\build_dna_rna_detail_route_comparison.py")
    }
}
finally {
    Pop-Location
}
