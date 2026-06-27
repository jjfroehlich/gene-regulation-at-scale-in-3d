param(
    [string]$PythonExe = "python",
    [string]$PyMolExe = "$env:LOCALAPPDATA\Schrodinger\PyMOL2\Scripts\pymol.exe",
    [string]$BlenderExe = "C:\Program Files\Blender Foundation\Blender 5.1\blender.exe",
    [switch]$ForceDownload,
    [switch]$SkipFetch,
    [switch]$SkipPyMolExport,
    [switch]$SkipReduction
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
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
    Invoke-Checked "Write canonical V4 manifest" {
        & $PythonExe "scripts\write_canonical_v4_manifest.py"
    }

    if (-not $SkipFetch) {
        $fetchArgs = @("scripts\fetch_rcsb_assets.py", "--skip-mcp")
        if ($ForceDownload) {
            $fetchArgs += "--force"
        }
        Invoke-Checked "Fetch RCSB assets" { & $PythonExe @fetchArgs }
    }

    if (-not $SkipPyMolExport) {
        $exportScript = (Join-Path $Root "scripts\export_pymol_surface_assets.py").Replace("\", "/")
        Invoke-Checked "Export PyMOL surface assets" { & $PyMolExe -cq -d "run $exportScript" }
    }

    if (-not $SkipReduction) {
        Invoke-Checked "Reduce surface assets in Blender" {
            & $BlenderExe --background --python (Join-Path $Root "scripts\reduce_surface_assets.py")
        }
    }

    Invoke-Checked "Build canonical V4 Blender scene" {
        & $BlenderExe --background --python (Join-Path $Root "scripts\build_gene_expression_surface_scene_v4.py")
    }
}
finally {
    Pop-Location
}
