param(
    [string]$RepoRoot,
    [string]$OutputDir,
    [string]$ExpectedPlatform,
    [string]$ExpectedArch,
    [switch]$KeepTemp
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
if (!$RepoRoot) {
    $RepoRoot = (Resolve-Path (Join-Path $ScriptDir "..\..")).ProviderPath
}
else {
    $RepoRoot = (Resolve-Path -LiteralPath $RepoRoot).ProviderPath
}

$runtime = [System.Runtime.InteropServices.RuntimeInformation]
$hostIsWindows = $runtime::IsOSPlatform([System.Runtime.InteropServices.OSPlatform]::Windows)
$hostIsLinux = $runtime::IsOSPlatform([System.Runtime.InteropServices.OSPlatform]::Linux)
$hostIsMacOS = $runtime::IsOSPlatform([System.Runtime.InteropServices.OSPlatform]::OSX)
if (!$hostIsWindows -and !$hostIsLinux -and !$hostIsMacOS) {
    throw "Unsupported build host: $($runtime::OSDescription)"
}

$platform = if ($hostIsWindows) { "windows" } elseif ($hostIsMacOS) { "macos" } else { "linux" }
$arch = $runtime::OSArchitecture.ToString().ToLowerInvariant()
if ($ExpectedPlatform -and $platform -ne $ExpectedPlatform.ToLowerInvariant()) {
    throw "Expected platform '$ExpectedPlatform', but the runner is '$platform'."
}
if ($ExpectedArch -and $arch -ne $ExpectedArch.ToLowerInvariant()) {
    throw "Expected architecture '$ExpectedArch', but the runner is '$arch'."
}
$sourceDir = Join-Path $RepoRoot "src\plugin"
if (!$OutputDir) {
    $OutputDir = Join-Path $RepoRoot "_dist\plugin\$platform-$arch"
}
$OutputDir = [System.IO.Path]::GetFullPath($OutputDir)

$tempRoot = Join-Path ([System.IO.Path]::GetTempPath()) (
    "rdistnn-plugin-build-" + [guid]::NewGuid().ToString("N")
)
$success = $false

try {
    New-Item -ItemType Directory -Path $tempRoot -Force | Out-Null
    foreach ($file in @(
        "Makefile",
        "Makefile.win",
        "rdistnn_core.hpp",
        "rdistnn_distance.hpp",
        "compensated_arithmetic.hpp",
        "rdistnn_plugin.cpp",
        "stplugin.c",
        "stplugin.h",
        "test_rdistnn_core.cpp",
        "test_rdistnn_blocks.cpp",
        "test_rdistnn_loo_plugin.cpp",
        "test_rdistnn_distance.cpp",
        "test_rdistnn_distance_plugin.cpp",
        "test_rdistnn_distance_rounding.cpp"
    )) {
        Copy-Item -LiteralPath (Join-Path $sourceDir $file) -Destination $tempRoot
    }

    if ($hostIsWindows) {
        $vswhere = Join-Path ${env:ProgramFiles(x86)} "Microsoft Visual Studio\Installer\vswhere.exe"
        if (!(Test-Path -LiteralPath $vswhere)) {
            throw "Visual Studio Installer's vswhere.exe was not found."
        }
        $vsRoot = (& $vswhere -latest -products "*" `
            -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
            -property installationPath | Select-Object -First 1)
        if ([string]::IsNullOrWhiteSpace($vsRoot)) {
            throw "Visual Studio C++ x64 build tools were not found."
        }
        $vcvars = Join-Path $vsRoot.Trim() "VC\Auxiliary\Build\vcvars64.bat"
        if (!(Test-Path -LiteralPath $vcvars)) {
            throw "vcvars64.bat was not found under $vsRoot."
        }

        $buildCmd = Join-Path $tempRoot "build.cmd"
        @"
@echo off
call "$vcvars"
if errorlevel 1 exit /b %errorlevel%
nmake /f Makefile.win clean
if errorlevel 1 exit /b %errorlevel%
nmake /f Makefile.win
if errorlevel 1 exit /b %errorlevel%
nmake /f Makefile.win check
exit /b %errorlevel%
"@ | Set-Content -LiteralPath $buildCmd -Encoding Ascii

        Push-Location $tempRoot
        try {
            & cmd.exe /d /c $buildCmd
            if ($LASTEXITCODE -ne 0) {
                throw "MSVC plugin build failed with exit code $LASTEXITCODE."
            }
        }
        finally {
            Pop-Location
        }
    }
    else {
        & make -C $tempRoot clean
        if ($LASTEXITCODE -ne 0) {
            throw "Plugin clean step failed with exit code $LASTEXITCODE."
        }
        & make -C $tempRoot
        if ($LASTEXITCODE -ne 0) {
            throw "Plugin build failed with exit code $LASTEXITCODE."
        }
        & make -C $tempRoot check
        if ($LASTEXITCODE -ne 0) {
            throw "Plugin core test failed with exit code $LASTEXITCODE."
        }
    }

    $builtPlugin = Join-Path $tempRoot "rdistnn.plugin"
    if (!(Test-Path -LiteralPath $builtPlugin)) {
        throw "Plugin build did not create rdistnn.plugin."
    }
    if ((Get-Item -LiteralPath $builtPlugin).Length -le 0) {
        throw "Built plugin is empty."
    }

    if ($hostIsWindows) {
        $bytes = [System.IO.File]::ReadAllBytes($builtPlugin)
        if ($bytes.Length -lt 2 -or $bytes[0] -ne 0x4d -or $bytes[1] -ne 0x5a) {
            throw "Windows plugin does not have a PE/MZ header."
        }
    }
    else {
        $fileDescription = (& file $builtPlugin 2>&1) -join "`n"
        if ($LASTEXITCODE -ne 0) {
            throw "Could not inspect the built plugin with file(1)."
        }
        if ($hostIsMacOS -and $fileDescription -notmatch "Mach-O") {
            throw "macOS plugin is not a Mach-O binary: $fileDescription"
        }
        if ($hostIsLinux -and $fileDescription -notmatch "ELF") {
            throw "Linux plugin is not an ELF binary: $fileDescription"
        }
    }

    New-Item -ItemType Directory -Path $OutputDir -Force | Out-Null
    $outputPlugin = Join-Path $OutputDir "rdistnn.plugin"
    Copy-Item -LiteralPath $builtPlugin -Destination $outputPlugin -Force
    $hash = (Get-FileHash -LiteralPath $outputPlugin -Algorithm SHA256).Hash

    Write-Host "Plugin build passed: $platform-$arch"
    Write-Host "Output: $outputPlugin"
    Write-Host "SHA256: $hash"
    if ($env:GITHUB_STEP_SUMMARY) {
        @"
## Stata plugin build

- Platform: ``$platform-$arch``
- Output: ``$outputPlugin``
- SHA-256: ``$hash``
"@ | Add-Content -LiteralPath $env:GITHUB_STEP_SUMMARY -Encoding UTF8
    }
    $success = $true
}
finally {
    if (!$KeepTemp -and (Test-Path -LiteralPath $tempRoot)) {
        Remove-Item -LiteralPath $tempRoot -Recurse -Force -ErrorAction SilentlyContinue
    }
    elseif ($KeepTemp) {
        Write-Host "Retained temporary build directory: $tempRoot"
    }
    if (!$success) {
        Write-Host "Plugin build check failed."
    }
}
