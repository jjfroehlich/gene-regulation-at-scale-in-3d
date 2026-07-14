param(
    [string]$BlenderExe = "C:\Program Files\Blender Foundation\Blender 5.1\blender.exe",
    [ValidateSet("v5", "v6")]
    [string]$CanonicalVersion = "v6"
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Blend = Join-Path $Root "outputs\canonical\gene_expression_surface_style_${CanonicalVersion}.blend"
$Exporter = Join-Path $Root "scripts\export_sketchfab_models.py"
$Fbx = Join-Path $Root "outputs\sketchfab\gene_expression_canonical_${CanonicalVersion}_sketchfab.fbx"

if (-not (Test-Path -LiteralPath $Blend)) { throw "Missing canonical blend: $Blend" }
if (-not (Test-Path -LiteralPath $BlenderExe)) { throw "Missing Blender executable: $BlenderExe" }

& $BlenderExe --background $Blend --python $Exporter -- --canonical-version $CanonicalVersion
if ($LASTEXITCODE -ne 0) { throw "Sketchfab model export failed with exit code $LASTEXITCODE" }

Get-Item (Join-Path $Root "outputs\sketchfab\gene_expression_canonical_${CanonicalVersion}_sketchfab.glb"), $Fbx |
    Select-Object FullName, Length, LastWriteTime
