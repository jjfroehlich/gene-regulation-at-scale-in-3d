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
    Invoke-Checked "Fetch required 1BNA/9IOB calibrators" {
        & $PythonExe "scripts\fetch_rcsb_assets.py" --skip-mcp --pdb-id 1BNA --pdb-id 9IOB
    }

    Invoke-Checked "Generate optimization V3 DNA/RNA CIFs" {
        & $PythonExe "experiments\procedural_nucleic_acids\scripts\optimization_v3_geometry.py"
    }

    Invoke-Checked "Generate canonical procedural comparison CIFs" {
        & $PythonExe "scripts\generate_procedural_nucleic_assets.py"
    }

    if (-not $SkipPyMolExport) {
        $canonicalExportScript = (Join-Path $Root "scripts\export_pymol_procedural_nucleic_surfaces.py").Replace("\", "/")
        $exportScript = (Join-Path $Root "experiments\procedural_nucleic_acids\scripts\export_optimization_v3_pymol_surfaces.py").Replace("\", "/")
        $calibratorExportScript = (Join-Path $Root "scripts\export_pymol_nucleic_calibrator_surfaces.py").Replace("\", "/")
        Invoke-Checked "Export canonical procedural comparison surfaces" {
            & $PyMolExe -cq -d "run $canonicalExportScript"
        }
        Invoke-Checked "Export optimization V3 PyMOL surfaces" {
            & $PyMolExe -cq -d "run $exportScript"
        }
        Invoke-Checked "Export 1BNA/9IOB calibrator surfaces" {
            & $PyMolExe -cq -d "run $calibratorExportScript"
        }

        $env:SURFACE_REDUCTION_INCLUDE_PROCEDURAL = "1"
        $env:SURFACE_REDUCTION_IDS = "DNA_PROXY,MRNA_PROXY,MRNA_COMPACT_PROXY"
        Invoke-Checked "Reduce canonical procedural comparison surfaces" {
            & $BlenderExe --background --python (Join-Path $Root "scripts\reduce_surface_assets.py")
        }
    }
    else {
        $requiredCachedAssets = @(
            "assets\pymol_exports\surface_assets_reduced\DNA_PROXY\DNA_PROXY_surface_strand_A.obj",
            "assets\pymol_exports\surface_assets_reduced\MRNA_PROXY\MRNA_PROXY_surface_coding.obj",
            "experiments\procedural_nucleic_acids\assets\optimization_v3\raw_surfaces\DNA_OPT_V3\DNA_OPT_V3_surface_strand_A.obj",
            "experiments\procedural_nucleic_acids\assets\optimization_v3\raw_surfaces\MRNA_ELONGATED_OPT_V3\MRNA_ELONGATED_OPT_V3_surface_coding.obj",
            "experiments\procedural_nucleic_acids\assets\optimization_v3\raw_surfaces\MRNA_COMPACT_OPT_V3\MRNA_COMPACT_OPT_V3_surface_coding.obj",
            "experiments\procedural_nucleic_acids\assets\pymol_calibrator_surfaces\1BNA\1BNA_nucleic_surface.obj",
            "experiments\procedural_nucleic_acids\assets\pymol_calibrator_surfaces\9IOB\9IOB_nucleic_surface.obj"
        )
        foreach ($relativePath in $requiredCachedAssets) {
            $assetPath = Join-Path $Root $relativePath
            if (-not (Test-Path -LiteralPath $assetPath)) {
                throw "-SkipPyMolExport requires cached asset: $relativePath"
            }
        }
    }

    $reportPath = Join-Path $Root "experiments\procedural_nucleic_acids\outputs\dna_rna_optimization_v3_report.json"
    if (Test-Path -LiteralPath $reportPath) {
        Remove-Item -LiteralPath $reportPath -Force
    }
    Invoke-Checked "Build optimization V3 Blender comparison" {
        & $BlenderExe --background --python (Join-Path $Root "experiments\procedural_nucleic_acids\scripts\build_dna_rna_optimization_v3.py")
    }
    if (-not (Test-Path -LiteralPath $reportPath)) {
        throw "Blender did not write the optimization V3 report; inspect the Blender traceback above."
    }
}
finally {
    Pop-Location
}
