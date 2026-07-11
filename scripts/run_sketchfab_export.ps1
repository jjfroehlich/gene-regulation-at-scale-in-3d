param(
    [string]$BlenderExe = "C:\Program Files\Blender Foundation\Blender 5.1\blender.exe"
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Blend = Join-Path $Root "outputs\canonical\gene_expression_surface_style_v5.blend"
$Exporter = Join-Path $Root "scripts\export_sketchfab_models.py"
$Fbx = Join-Path $Root "outputs\sketchfab\gene_expression_canonical_v5_sketchfab.fbx"

if (-not (Test-Path -LiteralPath $Blend)) { throw "Missing canonical blend: $Blend" }
if (-not (Test-Path -LiteralPath $BlenderExe)) { throw "Missing Blender executable: $BlenderExe" }

& $BlenderExe --background $Blend --python $Exporter
if ($LASTEXITCODE -ne 0) { throw "Sketchfab model export failed with exit code $LASTEXITCODE" }

Get-Item (Join-Path $Root "outputs\sketchfab\gene_expression_canonical_v5_sketchfab.glb"), $Fbx |
    Select-Object FullName, Length, LastWriteTime
