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
    $RepoRoot = (Resolve-Path (Join-Path $ScriptDir "..\..")).Path
}
else {
    $RepoRoot = (Resolve-Path -LiteralPath $RepoRoot).Path
}

$stata = Find-StataExecutable -StataExe $StataExe
$doFile = Join-Path $RepoRoot "tests\tests.do"
$logFile = Join-Path $RepoRoot "tests.log"

if (!$KeepLog -and (Test-Path -LiteralPath $logFile -PathType Leaf)) {
    Remove-Item -LiteralPath $logFile -Force
}

Write-Output "Using Stata executable: $stata"
Write-Output "Running Stata tests from: $RepoRoot"

$exitCode = Invoke-StataBatchDo -StataExe $stata -DoFile $doFile -WorkingDirectory $RepoRoot
Write-Output "Stata process exit code: $exitCode"

$passed = Assert-StataSuccessLog -LogPath $logFile
Write-Output "Stata runtime check passed: All $passed tests passed."
