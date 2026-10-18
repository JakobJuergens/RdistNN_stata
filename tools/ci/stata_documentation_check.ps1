param(
    [string]$StataExe,
    [string]$RepoRoot,
    [switch]$UpdateFigures
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
$buildDir = Join-Path $RepoRoot "build\documentation"
New-Item -ItemType Directory -Path $buildDir -Force | Out-Null
$doFile = Join-Path $buildDir "documentation_examples.do"
$logFile = Join-Path $RepoRoot "documentation_examples.log"
if (Test-Path -LiteralPath $logFile) { Remove-Item -LiteralPath $logFile -Force }
$repoPath = ConvertTo-StataPath -Path $RepoRoot
$template = @'
version 17
clear all
set more off
cd "__REPO__"
do "dev-setup.do"
help rdistnn
local passed 0
foreach example in getting-started cate-and-regression weights-and-tuning diagnostics-and-plots clate plqr-moments {
    capture noisily do "examples/`example'.do"
    local run_rc = _rc
    if `run_rc' {
        di as error "DOCUMENTATION_EXAMPLE_FAILED `example' rc=`run_rc'"
        exit `run_rc'
    }
    local ++passed
    di as result "DOCUMENTATION_EXAMPLE_OK `example'"
}
assert `passed'==6
if __UPDATE_FIGURES__ {
    capture mkdir "docs/figures"
    graph display docs_cate_profile
    graph export "docs/figures/cate-profile.png", width(1200) replace
    graph display docs_cate_dashboard
    graph export "docs/figures/cate-dashboard.png", width(1800) replace
    graph display docs_np_outcome
    graph export "docs/figures/np-outcome.png", width(1200) replace
    graph display docs_clate_est
    graph export "docs/figures/clate-forest.png", width(1200) replace
}
// Render installed-format SMCL without requiring a visible Results window.
translate "src/rdistnn.sthlp" "build/documentation/rdistnn-help.txt", translator(smcl2txt) replace
di as result "DOCUMENTATION_EXAMPLES_PASSED 6"
exit, clear
'@
$figureFlag = "0"
if ($UpdateFigures) { $figureFlag = "1" }
$doText = $template.Replace("__REPO__", $repoPath).Replace("__UPDATE_FIGURES__", $figureFlag)
Set-Content -LiteralPath $doFile -Value $doText -Encoding ASCII
Write-Output "Using Stata executable: $stata"
$exitCode = Invoke-StataBatchDo -StataExe $stata -DoFile $doFile -WorkingDirectory $RepoRoot
Write-Output "Stata process exit code: $exitCode"
if (!(Test-Path -LiteralPath $logFile)) { throw "Documentation log not found: $logFile" }
$content = Get-Content -LiteralPath $logFile -Raw
if ($content -notmatch '(?m)^DOCUMENTATION_EXAMPLES_PASSED 6\s*$' -or
    $content -match '(?m)^DOCUMENTATION_EXAMPLE_FAILED ') {
    throw "Documentation examples did not all pass. See $logFile"
}
Write-Output "All 6 documentation examples and overview help passed."
