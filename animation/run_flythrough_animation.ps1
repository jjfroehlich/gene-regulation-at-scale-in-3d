param(
    [string]$BlenderExe = "C:\Program Files\Blender Foundation\Blender 5.1\blender.exe",
    [string]$PythonExe = "python",
    [string]$FfmpegExe = "ffmpeg",
    [double]$DurationSeconds = 66,
    [int]$Fps = 24,
    [int]$ResolutionX = 1920,
    [int]$ResolutionY = 1080,
    [switch]$SmokeTest,
    [switch]$ReviewRender,
    [switch]$KeepFrames,
    [switch]$SkipVideoRender,
    [switch]$SkipReadmeGif
)
$ErrorActionPreference = "Stop"
if ($SmokeTest -and $ReviewRender) { throw "SmokeTest and ReviewRender are mutually exclusive." }
$runnerArgs = @((Join-Path $PSScriptRoot "scripts\run_dark_pipeline.py"), "--blender-exe", $BlenderExe, "--ffmpeg-exe", $FfmpegExe, "--duration-seconds", "$DurationSeconds", "--fps", "$Fps", "--resolution-x", "$ResolutionX", "--resolution-y", "$ResolutionY")
if ($SmokeTest) { $runnerArgs += "--smoke-test" }
if ($ReviewRender) { $runnerArgs += "--review-render" }
if ($KeepFrames) { $runnerArgs += "--keep-frames" }
if ($SkipVideoRender) { $runnerArgs += "--skip-video-render" }
if ($SkipReadmeGif) { $runnerArgs += "--skip-readme-gif" }
& $PythonExe @runnerArgs
if ($LASTEXITCODE -ne 0) { throw "Flythrough pipeline failed: $LASTEXITCODE" }
