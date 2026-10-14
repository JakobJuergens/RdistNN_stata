Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function ConvertTo-StataPath {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Path
    )
    return $Path.Replace("\", "/")
}

function Resolve-StataCandidate {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Candidate
    )

    if (Test-Path -LiteralPath $Candidate -PathType Leaf) {
        return (Resolve-Path -LiteralPath $Candidate).ProviderPath
    }

    $cmd = Get-Command $Candidate -ErrorAction SilentlyContinue
    if ($null -ne $cmd -and $cmd.CommandType -eq "Application") {
        return $cmd.Source
    }

    return $null
}

function Set-RefreshedPath {
    $pieces = @(
        [Environment]::GetEnvironmentVariable("Path", "User"),
        [Environment]::GetEnvironmentVariable("Path", "Machine"),
        $env:Path
    ) | Where-Object { $_ -and $_.Trim().Length -gt 0 }

    $env:Path = ($pieces -join ";")
}

function Find-StataExecutable {
    [CmdletBinding()]
    param(
        [string]$StataExe
    )

    if ($StataExe) {
        $resolved = Resolve-StataCandidate -Candidate $StataExe
        if ($resolved) {
            return $resolved
        }
        throw "Requested Stata executable was not found: $StataExe"
    }

    if ($env:STATA_EXE) {
        $resolved = Resolve-StataCandidate -Candidate $env:STATA_EXE
        if ($resolved) {
            return $resolved
        }
        throw "STATA_EXE is set but was not found: $env:STATA_EXE"
    }

    Set-RefreshedPath

    $commandCandidates = @(
        "StataMP-64.exe",
        "StataSE-64.exe",
        "StataBE-64.exe",
        "Stata-64.exe",
        "StataMP.exe",
        "StataSE.exe",
        "StataBE.exe",
        "Stata.exe",
        "stata"
    )

    foreach ($candidate in $commandCandidates) {
        $resolved = Resolve-StataCandidate -Candidate $candidate
        if ($resolved) {
            return $resolved
        }
    }

    $programDirs = @(
        [Environment]::GetEnvironmentVariable("ProgramFiles"),
        [Environment]::GetEnvironmentVariable("ProgramFiles(x86)")
    ) | Where-Object { $_ -and (Test-Path -LiteralPath $_ -PathType Container) }

    $stataDirs = @(
        "StataNow19",
        "Stata19",
        "StataNow18",
        "Stata18",
        "StataNow17",
        "Stata17"
    )

    foreach ($programDir in $programDirs) {
        foreach ($stataDir in $stataDirs) {
            foreach ($exeName in $commandCandidates) {
                $candidate = Join-Path (Join-Path $programDir $stataDir) $exeName
                $resolved = Resolve-StataCandidate -Candidate $candidate
                if ($resolved) {
                    return $resolved
                }
            }
        }
    }

    throw "No Stata executable found. Pass -StataExe or set STATA_EXE."
}

function Invoke-StataBatchDo {
    [CmdletBinding()]
    param(
        [Parameter(Mandatory = $true)]
        [string]$StataExe,

        [Parameter(Mandatory = $true)]
        [string]$DoFile,

        [Parameter(Mandatory = $true)]
        [string]$WorkingDirectory
    )

    if (!(Test-Path -LiteralPath $DoFile -PathType Leaf)) {
        throw "Stata do-file not found: $DoFile"
    }
    if (!(Test-Path -LiteralPath $WorkingDirectory -PathType Container)) {
        throw "Working directory not found: $WorkingDirectory"
    }

    $process = Start-Process `
        -FilePath $StataExe `
        -ArgumentList @("/e", "do", ('"' + $DoFile + '"')) `
        -WorkingDirectory $WorkingDirectory `
        -Wait `
        -PassThru `
        -WindowStyle Hidden

    return $process.ExitCode
}

function Assert-StataSuccessLog {
    [CmdletBinding()]
    param(
        [Parameter(Mandatory = $true)]
        [string]$LogPath
    )

    if (!(Test-Path -LiteralPath $LogPath -PathType Leaf)) {
        throw "Expected Stata log not found: $LogPath"
    }

    $content = Get-Content -LiteralPath $LogPath -Raw

    if ($content -match "(?m)^\s+\[runner\].*aborted unexpectedly") {
        throw "Stata test runner aborted unexpectedly. See log: $LogPath"
    }
    if ($content -match "(?m)^\d+\s+of\s+\d+\s+tests\s+FAILED\.") {
        throw "Stata test failures found. See log: $LogPath"
    }

    $match = [regex]::Match($content, "All\s+(\d+)\s+tests\s+passed\.")
    if (!$match.Success) {
        throw "Stata success marker not found in log: $LogPath"
    }

    return [int]$match.Groups[1].Value
}
