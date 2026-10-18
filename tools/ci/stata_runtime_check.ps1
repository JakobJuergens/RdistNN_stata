param(
    [string]$StataExe,
    [string]$RepoRoot,
    [switch]$KeepLog
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
. (Join-Path $ScriptDir "stata_common.ps1")

if (!$RepoRoot) {
    $RepoRoot = (Resolve-Path (Join-Path $ScriptDir "..\..")).ProviderPath
}
else {
    $RepoRoot = (Resolve-Path -LiteralPath $RepoRoot).ProviderPath
}

$stata = Find-StataExecutable -StataExe $StataExe
$doFile = Join-Path $RepoRoot "tests\tests.do"
$logFile = Join-Path $RepoRoot "tests.log"

if (Test-Path -LiteralPath $logFile -PathType Leaf) {
    if ($KeepLog) {
        Copy-Item -LiteralPath $logFile -Destination (Join-Path $RepoRoot "tests.previous.log") -Force
    }
    Remove-Item -LiteralPath $logFile -Force
}

Write-Output "Using Stata executable: $stata"
Write-Output "Running Stata tests from: $RepoRoot"

$exitCode = Invoke-StataBatchDo -StataExe $stata -DoFile $doFile -WorkingDirectory $RepoRoot
Write-Output "Stata process exit code: $exitCode"
if ($exitCode -ne 0) {
    throw "Stata runtime process failed with exit code $exitCode. Log: $logFile"
}

$passed = Assert-StataSuccessLog -LogPath $logFile
Write-Output "Stata runtime check passed: All $passed tests passed."
