param(
    [string]$BlenderExe = "C:\Program Files\Blender Foundation\Blender 5.1\blender.exe",
    [string]$PythonExe = "python",
    [string]$FfmpegExe = "ffmpeg",
    [double]$DurationSeconds = 60,
    [int]$Fps = 24,
    [int]$ResolutionX = 1920,
    [int]$ResolutionY = 1080,
    [switch]$SmokeTest,
    [switch]$SkipVideoRender,
    [switch]$SkipReadmeGif
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$ScriptPath = Join-Path $PSScriptRoot "scripts\build_v5_flythrough_animation.py"
$GifFallbackScript = Join-Path $PSScriptRoot "scripts\build_readme_gif.py"
$OutputDir = Join-Path $PSScriptRoot "outputs"
$VideoPath = Join-Path $OutputDir ($(if ($SmokeTest) { "v5_flythrough_animation_smoke.mp4" } else { "v5_flythrough_animation_1080p.mp4" }))
$ReportPath = Join-Path $OutputDir "v5_flythrough_animation_report.json"
$ReadmeGif = Join-Path $Root "docs\images\v5-flythrough-preview.gif"

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

Push-Location $Root
try {
    New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null
    if (Test-Path -LiteralPath $ReportPath) {
        Remove-Item -LiteralPath $ReportPath -Force
    }
    if ((-not $SkipVideoRender) -and (Test-Path -LiteralPath $VideoPath)) {
        Remove-Item -LiteralPath $VideoPath -Force
    }
    $encodedVideo = $false
    $blenderArgs = @(
        "--background",
        "--python", $ScriptPath,
        "--",
        "--output-dir", $OutputDir,
        "--duration-seconds", "$DurationSeconds",
        "--fps", "$Fps",
        "--resolution-x", "$ResolutionX",
        "--resolution-y", "$ResolutionY",
        "--output-mp4", $VideoPath
    )
    if ($SmokeTest) {
        $blenderArgs += "--smoke-test"
    }
    if ($SkipVideoRender) {
        $blenderArgs += "--skip-video-render"
    }
    Invoke-Checked "Build V5 flythrough animation" {
        & $BlenderExe @blenderArgs
    }
    if (-not (Test-Path -LiteralPath $ReportPath)) {
        throw "Animation report was not written: $ReportPath"
    }
    $report = Get-Content -Raw -LiteralPath $ReportPath | ConvertFrom-Json

    if ((-not $SkipVideoRender) -and $report.rendered) {
        $framePattern = Join-Path $report.frames_dir "frame_%04d.png"
        Invoke-Checked "Encode MP4 from rendered frames" {
            & $FfmpegExe -y -framerate $report.fps -start_number 1 -i $framePattern -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -pix_fmt yuv420p -crf 19 $VideoPath
        }
        $encodedVideo = $true
    }

    if (-not $SkipReadmeGif) {
        if ($encodedVideo -and (Test-Path -LiteralPath $VideoPath)) {
            $gifFilter = if ($SmokeTest) {
                "fps=6,scale=720:-1:flags=lanczos,split[s0][s1];[s0]palettegen=max_colors=96[p];[s1][p]paletteuse=dither=bayer:bayer_scale=5"
            }
            else {
                "fps=1,setpts=N/(6*TB),scale=720:-1:flags=lanczos,split[s0][s1];[s0]palettegen=max_colors=96[p];[s1][p]paletteuse=dither=bayer:bayer_scale=5"
            }
            Invoke-Checked "Build README GIF from rendered video" {
                & $FfmpegExe -y -i $VideoPath -vf $gifFilter $ReadmeGif
            }
        }
        else {
            Invoke-Checked "Build README GIF fallback from still renders" {
                & $PythonExe $GifFallbackScript --output $ReadmeGif
            }
        }
    }
}
finally {
    Pop-Location
}
