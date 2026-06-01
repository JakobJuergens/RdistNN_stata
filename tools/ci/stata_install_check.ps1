param(
    [string]$StataExe,
    [string]$RepoRoot,
    [switch]$KeepTemp
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
$tempRoot = Join-Path ([System.IO.Path]::GetTempPath()) ("rdistnn-install-check-" + [guid]::NewGuid().ToString("N"))
$plusDir = Join-Path $tempRoot "plus"
$personalDir = Join-Path $tempRoot "personal"
$doFile = Join-Path $tempRoot "install_check.do"
$logFile = Join-Path $tempRoot "install_check.log"

$repoFrom = ConvertTo-StataPath -Path $RepoRoot
$repoSrcNorm = (ConvertTo-StataPath -Path (Join-Path $RepoRoot "src")).ToLowerInvariant()
$plusPath = ConvertTo-StataPath -Path $plusDir
$personalPath = ConvertTo-StataPath -Path $personalDir
$plusNorm = $plusPath.ToLowerInvariant()

$success = $false

try {
    New-Item -ItemType Directory -Path $plusDir -Force | Out-Null
    New-Item -ItemType Directory -Path $personalDir -Force | Out-Null

    $doTemplate = @'
version 17
clear all
set more off

local repo_from "__REPO_FROM__"
local repo_src_norm "__REPO_SRC_NORM__"
local plus_path "__PLUS_PATH__"
local personal_path "__PERSONAL_PATH__"
local plus_norm "__PLUS_NORM__"

sysdir set PLUS "`plus_path'"
sysdir set PERSONAL "`personal_path'"

capture findfile np_dnn.ado
if !_rc {
    local found `"`r(fn)'"'
    local found_norm = lower(subinstr(`"`found'"', `"\"', `"/"', .))
    if strpos(`"`found_norm'"', `"`repo_src_norm'"') {
        local slash = max(strrpos(`"`found'"', `"/"'), strrpos(`"`found'"', `"\"'))
        local dir = substr(`"`found'"', 1, `slash' - 1)
        di as text "Removing development adopath entry: " as result `"`dir'"'
        adopath - `"`dir'"'
    }
}

net install rdistnn, from("`repo_from'") replace

capture program drop _rdistnn_assert_from_plus
program define _rdistnn_assert_from_plus
    syntax anything(name=filename), Plusnorm(string)
    capture findfile `filename'
    if _rc {
        di as error "`filename' not found after net install."
        exit 601
    }
    local found `"`r(fn)'"'
    local found_norm = lower(subinstr(`"`found'"', `"\"', `"/"', .))
    if strpos(`"`found_norm'"', `"`plusnorm'"') != 1 {
        di as error "`filename' resolved outside temporary PLUS directory:"
        di as error `"`found'"'
        exit 459
    }
    di as text "`filename' resolved from " as result `"`found'"'
end

which weights_dnn
_rdistnn_assert_from_plus weights_dnn.ado, plusnorm("`plus_norm'")
which np_dnn
_rdistnn_assert_from_plus np_dnn.ado, plusnorm("`plus_norm'")
which cate_dnn
_rdistnn_assert_from_plus cate_dnn.ado, plusnorm("`plus_norm'")
which plqr_dnn
_rdistnn_assert_from_plus plqr_dnn.ado, plusnorm("`plus_norm'")
_rdistnn_assert_from_plus _rdistnn_mata.ado, plusnorm("`plus_norm'")
_rdistnn_assert_from_plus rdistnn.plugin, plusnorm("`plus_norm'")

findfile _rdistnn_mata.ado
quietly run "`r(fn)'"
if "$rdistnn_plugin_ok" != "1" {
    di as error "rdistnn.plugin was installed but not discoverable by Stata."
    exit 459
}

clear
set obs 160
set seed 90210
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double e = invlogit(0.25*x1 - 0.15*x2)
gen byte D = runiform() < e
gen double m0 = 1 + 0.4*x1 - 0.2*x2
gen double m1 = m0 + 0.5 + 0.1*x1
gen double y = m0 + D*(m1 - m0) + rnormal()*0.25

weights_dnn x1 x2, x0(0 0) k(3) s(12) generate(w_install)
if r(neff) >= . | r(neff) <= 1 {
    di as error "weights_dnn returned invalid r(neff)."
    exit 459
}

np_dnn y, x(x1 x2) x0(0 0) k(3) s(12)
if _b[mu] >= . | _se[mu] >= . {
    di as error "np_dnn returned non-finite results."
    exit 459
}

cate_dnn y, treat(D) x(x1 x2) x0(0 0) propensity(e) m0(m0) m1(m1) k(3) s(12)
if _b[tau] >= . | _se[tau] >= . {
    di as error "cate_dnn returned non-finite results."
    exit 459
}

gen double md1 = 0.3*x1 - 0.2*x2
gen double md2 = -0.1*x1 + 0.4*x2
gen double d1 = md1 + rnormal()*0.4
gen double d2 = md2 + rnormal()*0.4
gen double q = 0.25*x1 - 0.15*x2
gen double yq = d1 - 0.5*d2 + q + rnormal()*0.5

plqr_dnn yq, x(x1 x2) x0(0 0) d(d1 d2) q(q) md(md1 md2) beta(1 -0.5) tau(0.5) k(3) s(12)
matrix B = e(b)
matrix V = e(V)
if B[1,1] >= . | B[1,2] >= . | V[1,1] >= . | V[2,2] >= . {
    di as error "plqr_dnn returned non-finite results."
    exit 459
}

di as result "INSTALL_CHECK_PASSED"
exit
'@

    $doText = $doTemplate.
        Replace("__REPO_FROM__", $repoFrom).
        Replace("__REPO_SRC_NORM__", $repoSrcNorm).
        Replace("__PLUS_PATH__", $plusPath).
        Replace("__PERSONAL_PATH__", $personalPath).
        Replace("__PLUS_NORM__", $plusNorm)

    Set-Content -LiteralPath $doFile -Value $doText -Encoding ASCII

    Write-Output "Using Stata executable: $stata"
    Write-Output "Running fresh install check in: $tempRoot"

    $exitCode = Invoke-StataBatchDo -StataExe $stata -DoFile $doFile -WorkingDirectory $tempRoot
    Write-Output "Stata process exit code: $exitCode"

    if (!(Test-Path -LiteralPath $logFile -PathType Leaf)) {
        throw "Expected Stata install-check log not found: $logFile"
    }

    $content = Get-Content -LiteralPath $logFile -Raw
    if ($content -notmatch "(?m)^INSTALL_CHECK_PASSED\r?$") {
        throw "Fresh install check did not reach INSTALL_CHECK_PASSED. Log: $logFile"
    }

    $success = $true
    Write-Output "Fresh install check passed."
}
finally {
    if ($success -and !$KeepTemp) {
        Remove-Item -LiteralPath $tempRoot -Recurse -Force
    }
    else {
        Write-Output "Install-check temp directory retained: $tempRoot"
    }
}
