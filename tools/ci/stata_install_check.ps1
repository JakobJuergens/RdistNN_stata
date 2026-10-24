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
    $RepoRoot = (Resolve-Path (Join-Path $ScriptDir "..\..")).ProviderPath
}
else {
    $RepoRoot = (Resolve-Path -LiteralPath $RepoRoot).ProviderPath
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

_rdistnn_assert_from_plus rdistnn.sthlp, plusnorm("`plus_norm'")
help rdistnn
which weights_dnn
_rdistnn_assert_from_plus weights_dnn.ado, plusnorm("`plus_norm'")
which np_dnn
_rdistnn_assert_from_plus np_dnn.ado, plusnorm("`plus_norm'")
which np_dnn_diagnostics
_rdistnn_assert_from_plus np_dnn_diagnostics.ado, plusnorm("`plus_norm'")
which np_dnn_plot
_rdistnn_assert_from_plus np_dnn_plot.ado, plusnorm("`plus_norm'")
which cate_dnn
_rdistnn_assert_from_plus cate_dnn.ado, plusnorm("`plus_norm'")
which cate_dnn_diagnostics
_rdistnn_assert_from_plus cate_dnn_diagnostics.ado, plusnorm("`plus_norm'")
which cate_dnn_plot
_rdistnn_assert_from_plus cate_dnn_plot.ado, plusnorm("`plus_norm'")
which plqr_dnn
_rdistnn_assert_from_plus plqr_dnn.ado, plusnorm("`plus_norm'")
foreach command in clate_dnn clate_dnn_diagnostics clate_dnn_plot rdistnn_results rdistnn_candidate rdistnn_predict np_tdnn weights_tdnn np_tdnn_diagnostics np_tdnn_plot {
    which `command'
    _rdistnn_assert_from_plus `command'.ado, plusnorm("`plus_norm'")
    _rdistnn_assert_from_plus `command'.sthlp, plusnorm("`plus_norm'")
}
_rdistnn_assert_from_plus _rdistnn_clate_mata.ado, plusnorm("`plus_norm'")
_rdistnn_assert_from_plus _rdistnn_mata.ado, plusnorm("`plus_norm'")
_rdistnn_assert_from_plus _rdistnn_distance_mata.ado, plusnorm("`plus_norm'")
_rdistnn_assert_from_plus _rdistnn_plot_helpers.ado, plusnorm("`plus_norm'")
foreach helper in report report_header estimation_state sample_signature diagnostics diagnostics_save plot_dispatch model_store model_check model_mata model_fit model_plot reuse tdnn_mata rank_weights rank_loo distance_plugin native_transport accumulation covariance geometry plot_mata {
    _rdistnn_assert_from_plus _rdistnn_`helper'.ado, plusnorm("`plus_norm'")
}
_rdistnn_assert_from_plus rdistnn_license.md, plusnorm("`plus_norm'")
_rdistnn_assert_from_plus rdistnn_linking_exception.md, plusnorm("`plus_norm'")
_rdistnn_assert_from_plus rdistnn.plugin, plusnorm("`plus_norm'")

findfile _rdistnn_mata.ado
quietly run "`r(fn)'"
if "$rdistnn_plugin_ok" != "1" {
    di as error "rdistnn.plugin was installed but not discoverable by Stata."
    exit 459
}

// Run public onboarding before other plots so CATE-first autoloading is cold.
do "__REPO_FROM__/examples/getting-started.do"

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

global rdistnn_use_plugin 0
weights_dnn x1 x2, x0(0.2 -0.3) k(3) s(12) method(exact) generate(w_install_mata)
global rdistnn_use_plugin 1
weights_dnn x1 x2, x0(0.2 -0.3) k(3) s(12) method(exact) generate(w_install_plugin)
capture macro drop rdistnn_use_plugin
gen double w_install_diff = abs(w_install_mata - w_install_plugin) ///
    if !missing(w_install_mata, w_install_plugin)
quietly summarize w_install_diff, meanonly
if r(N) == 0 | r(max) > 1e-12 {
    di as error "installed plugin weights do not match Mata for signed x0()."
    exit 459
}

// Require the new distance API from the installed binary, not checkout sources.
matrix install_distance_api=(.)
plugin call rdistnn_plugin, distance_api install_distance_api
if install_distance_api[1,1] != 1 exit 459
foreach metric in l1 l2 l1_scaled l2_scaled {
    global rdistnn_distance_engine mata
    quietly weights_dnn x1 x2, x0(0.2 -0.3) k(3) s(12) distance(`metric')
    matrix install_distances=r(distances)
    matrix install_weights=r(weights)
    global rdistnn_distance_engine plugin
    quietly weights_dnn x1 x2, x0(0.2 -0.3) k(3) s(12) distance(`metric')
    if "$rdistnn_last_distance_engine" != "plugin" exit 459
    mata: assert(max(abs(st_matrix("r(distances)")-st_matrix("install_distances"))) < 1e-12)
    mata: assert(max(abs(st_matrix("r(weights)")-st_matrix("install_weights"))) < 1e-12)
}
quietly np_dnn y, xcols(x1 x2) x0(0 0) k(3) s(12) distance(l1_scaled)
tempfile native_install_fit
quietly estimates save "`native_install_fit'", replace
preserve
clear
quietly estimates use "`native_install_fit'"
matrix install_new_targets=(0.2,-0.3 \ 0.4,0.1)
quietly rdistnn_predict, x0matrix(install_new_targets) sefit
if "$rdistnn_last_distance_engine" != "plugin" | _N != 0 exit 459
matrix install_prediction=r(fit)
mata: assert(!hasmissing(st_matrix("install_prediction")))
restore
capture macro drop rdistnn_distance_engine
capture macro drop rdistnn_last_distance_engine

// Require the installed original-signal block API, including its native path.
// A rounded block mean would lose the unit survivor after deleting 2^54.
mata {
    // Explicit represented metric distances tie observations1 and2.
    install_block_distances=(0\0\1)
    install_block_loo=rdistnn_loo_from_table((1::3),(2^54\1\0),(1\0),0,"plugin",install_block_distances)
    assert(all(install_block_loo:==(1\2^54\2^53)))
    assert(st_global("rdistnn_last_jackknife_engine")=="plugin")
    install_eta=rdistnn_power2(-1074)
    assert(install_eta>0)
    install_block_loo=rdistnn_loo_from_table((1::3),(install_eta\0\0),(1\0),0,"plugin",install_block_distances)
    assert(all(install_block_loo:==(0\install_eta\0)))
}

np_dnn y, xcols(x1 x2) x0(0 0) k(3) s(12)
if _b[mu] >= . | _se[mu] >= . {
    di as error "np_dnn returned non-finite results."
    exit 459
}
np_dnn y, xcols(x1 x2) x0(0 0)
if e(k) != 1 | e(s) != floor(.4*160^(3/4)) | abs(e(lambda)-.4) > 1e-12 {
    di as error "installed np_dnn does not use the R/paper k and s defaults."
    exit 459
}

matrix X0_np_install = (0,0 \ 0.2,-0.3)
np_dnn y, xcols(x1 x2) x0matrix(X0_np_install) k(3) tuningmode(grid)
matrix B_np_batch = e(b)
matrix V_np_batch = e(V)
matrix D_np_batch = e(diagnostics)
matrix C_np_grid = e(tuning_candidates)
matrix E_np_grid = e(tuning_estimates)
if e(m) != 2 | colsof(B_np_batch) != 2 | rowsof(V_np_batch) != 2 | ///
    rowsof(D_np_batch) != 2 | colsof(D_np_batch) != 17 | ///
    rowsof(C_np_grid) != 4 | rowsof(E_np_grid) != 4 | ///
    e(tuning_default_candidate) != 2 | "`e(tuning_mode)'" != "grid" {
    di as error "installed np_dnn batch, diagnostics, or tuning contract failed."
    exit 459
}
np_dnn_diagnostics, allcandidates
if r(n_rows) != 4 | "`r(scope)'" != "candidates" {
    di as error "installed np_dnn_diagnostics candidate path failed."
    exit 459
}
foreach p in estimates diagnostics tuning weights outcome dashboard {
    local extra ""
    if "`p'" == "weights" local extra "view(kernel)"
    np_dnn_plot, type(`p') `extra' name(install_np_`p') replace nodraw
    capture graph describe install_np_`p'
    if _rc {
        di as error "installed np_dnn_plot type(`p') failed."
        exit 459
    }
}
np_dnn y, xcols(x1 x2) x0matrix(X0_np_install) s(160) nose
capture matrix list e(V)
if _rc == 0 | "`e(se_method)'" != "none" {
    di as error "installed np_dnn nose path unexpectedly posted inference."
    exit 459
}

cate_dnn y, tcol(D) xcols(x1 x2) x0(0 0) ecol(e) m0col(m0) m1col(m1) k(3) s(12)
if _b[tau] >= . | _se[tau] >= . {
    di as error "cate_dnn returned non-finite results."
    exit 459
}
local tau_aipw = _b[tau]
local se_aipw = _se[tau]

gen double eclip = min(max(e, 1e-8), 1 - 1e-8)
gen double rho_aipw = D*(y-m1)/eclip - (1-D)*(y-m0)/(1-eclip) + m1-m0
cate_dnn rho_aipw, xcols(x1 x2) x0(0 0) pseudooutcome k(3) s(12)
if abs(_b[tau] - `tau_aipw') > 1e-10 | abs(_se[tau] - `se_aipw') > 1e-10 {
    di as error "installed cate_dnn pseudooutcome path does not match AIPW."
    exit 459
}
cate_dnn rho_aipw, xcols(x1 x2) x0(0 0) pseudooutcome
if e(k) != 1 | e(s) != floor(.4*160^(3/4)) | abs(e(lambda)-.4) > 1e-12 {
    di as error "installed cate_dnn does not use the R/paper k and s defaults."
    exit 459
}

matrix X0_install = (0,0 \ 0.2,-0.3)
cate_dnn y, tcol(D) xcols(x1 x2) x0matrix(X0_install) ecol(e) ///
    m0col(m0) m1col(m1) k(3) s(12)
matrix B_cate_batch = e(b)
matrix V_cate_batch = e(V)
matrix T_cate_batch = e(targets)
matrix N_cate_batch = e(neff_targets)
matrix D_cate_batch = e(diagnostics)
if e(m) != 2 | colsof(B_cate_batch) != 2 | rowsof(V_cate_batch) != 2 | ///
    rowsof(T_cate_batch) != 2 | colsof(N_cate_batch) != 2 | ///
    rowsof(D_cate_batch) != 2 | colsof(D_cate_batch) != 25 | ///
    D_cate_batch[1,12] != 1 | D_cate_batch[1,16] != 1 {
    di as error "installed cate_dnn x0matrix() path returned invalid result shapes."
    exit 459
}
if abs(B_cate_batch[1,1] - `tau_aipw') > 1e-10 | ///
    abs(sqrt(V_cate_batch[1,1]) - `se_aipw') > 1e-10 {
    di as error "installed cate_dnn x0matrix() first target does not match x0()."
    exit 459
}
cate_dnn_diagnostics
if r(n_rows) != 2 | "`r(scope)'" != "targets" {
    di as error "installed cate_dnn_diagnostics target path failed."
    exit 459
}
foreach p in estimates diagnostics overlap weights influence {
    local extra ""
    if "`p'" == "overlap" local extra "view(mass)"
    if "`p'" == "weights" local extra "view(kernel)"
    cate_dnn_plot, type(`p') `extra' name(install_`p') replace nodraw
    capture graph describe install_`p'
    if _rc {
        di as error "installed cate_dnn_plot type(`p') failed."
        exit 459
    }
}

cate_dnn rho_aipw, xcols(x1 x2) x0matrix(X0_install) pseudooutcome ///
    tuningmode(grid)
matrix C_cate_grid = e(tuning_candidates)
matrix E_cate_grid = e(tuning_estimates)
matrix S_cate_grid = e(tuning_se)
matrix D_cate_grid = e(tuning_diagnostics)
if rowsof(C_cate_grid) != 4 | rowsof(E_cate_grid) != 4 | ///
    rowsof(S_cate_grid) != 4 | rowsof(D_cate_grid) != 4 | ///
    e(tuning_default_candidate) != 2 | "`e(tuning_mode)'" != "grid" {
    di as error "installed cate_dnn tuning-grid path returned invalid results."
    exit 459
}
local s_grid_default = e(s)
matrix B_grid_default = e(b)
matrix V_grid_default = e(V)
cate_dnn rho_aipw, xcols(x1 x2) x0matrix(X0_install) pseudooutcome ///
    s(`s_grid_default')
matrix B_grid_scalar = e(b)
matrix V_grid_scalar = e(V)
mata {
    st_numscalar("_install_grid_b", max(abs(st_matrix("B_grid_default") - st_matrix("B_grid_scalar"))))
    st_numscalar("_install_grid_v", max(abs(st_matrix("V_grid_default") - st_matrix("V_grid_scalar"))))
}
if scalar(_install_grid_b) > 1e-10 | scalar(_install_grid_v) > 1e-10 {
    di as error "installed cate_dnn grid default does not match scalar fit."
    exit 459
}

cate_dnn y, tcol(D) xcols(x1 x2) x0matrix(X0_install) ecol(e) ///
    m0col(m0) m1col(m1) k(3) tuningmode(grid)
cate_dnn_diagnostics, allcandidates
if r(n_rows) != 4 | "`r(scope)'" != "candidates" {
    di as error "installed cate_dnn_diagnostics candidate path failed."
    exit 459
}
cate_dnn_plot, type(tuning) allcandidates name(install_tuning) replace nodraw
cate_dnn_plot, type(dashboard) name(install_dashboard) replace nodraw
capture graph describe install_dashboard
if _rc {
    di as error "installed cate_dnn_plot dashboard path failed."
    exit 459
}
capture graph drop install_estimates install_diagnostics install_overlap ///
    install_weights install_influence install_tuning install_dashboard ///
    install_np_estimates install_np_diagnostics install_np_tuning ///
    install_np_weights install_np_outcome install_np_dashboard

cate_dnn rho_aipw, xcols(x1 x2) x0(0 0) pseudooutcome s(160) nose
if _b[tau] >= . | e(k) != 1 | "`e(se_method)'" != "none" {
    di as error "installed cate_dnn nose/default-k path returned invalid results."
    exit 459
}
capture matrix list e(V)
if !_rc {
    di as error "installed cate_dnn nose path unexpectedly posted e(V)."
    exit 459
}

gen double md1 = 0.3*x1 - 0.2*x2
gen double md2 = -0.1*x1 + 0.4*x2
gen double d1 = md1 + rnormal()*0.4
gen double d2 = md2 + rnormal()*0.4
gen double q = 0.25*x1 - 0.15*x2
gen double yq = d1 - 0.5*d2 + q + rnormal()*0.5

plqr_dnn yq, x(x1 x2) x0(0 0) d(d1 d2) g(q) md(md1 md2) thetaeval(1 -0.5) tau(0.5) k(3) s(12)
matrix B = e(b)
matrix V = e(V)
if B[1,1] >= . | B[1,2] >= . | V[1,1] >= . | V[2,2] >= . {
    di as error "plqr_dnn returned non-finite results."
    exit 459
}
plqr_dnn yq, x(x1 x2) x0(0 0) d(d1 d2) g(q) md(md1 md2) thetaeval(1 -0.5)
if e(k) != 1 | e(s) != floor(.4*160^(3/4)) {
    di as error "installed plqr_dnn does not use the R/paper k and s defaults."
    exit 459
}

// The CLATE example exercises only freshly installed files, including cold
// paired evaluator loading, diagnostics and instrument/influence/tuning plots.
do "__REPO_FROM__/examples/reusing-fits.do"
if e(model_version)!=4 | e(m)!=3 {
    di as error "installed saved-fit prediction example failed."
    exit 459
}
do "__REPO_FROM__/examples/clate.do"
if "`e(cmd)'" != "clate_dnn" | e(m) != 3 {
    di as error "installed CLATE example did not retain its batch results."
    exit 459
}
matrix CL_install = e(estimates)
matrix CL_cov = e(theta_vcov)
if missing(CL_install[1,1]) | missing(CL_cov[1,1]) {
    di as error "installed CLATE results are not finite."
    exit 459
}
clate_dnn_diagnostics, allcandidates
foreach p in diagnostics weights dashboard {
    clate_dnn_plot, type(`p') name(install_clate_`p') replace nodraw
}

// Check real summary/localization/plot tables from the isolated installation.
foreach family in clate cate np {
    if "`family'"=="cate" cate_dnn Y, xcols(x1 x2) x0matrix(X0) pseudooutcome
    if "`family'"=="np" np_dnn Y, xcols(x1 x2) x0matrix(X0)
    `family'_dnn, level(90)
    matrix REPORT_INSTALL = r(table)
    matrix SUMMARY_INSTALL = r(summary)
    assert rowsof(REPORT_INSTALL)==9 & colsof(REPORT_INSTALL)==3
    assert rowsof(SUMMARY_INSTALL)==3 & colsof(SUMMARY_INSTALL)==6
    `family'_dnn_diagnostics, targets(3 1)
    matrix LOCAL_INSTALL = r(localization)
    assert rowsof(LOCAL_INSTALL)==2 & LOCAL_INSTALL[1,1]==3
    `family'_dnn_plot, targets(3 1) name(install_common) replace nodraw
    matrix PLOT_INSTALL = r(plot_data)
    assert rowsof(PLOT_INSTALL)==2
}
// All eight public examples run from installed commands: quick start, saved
// reuse and CLATE ran above; these complete the article collection.
foreach example in cate-and-regression weights-and-tuning diagnostics-and-plots plqr-moments two-scale-regression {
    do "__REPO_FROM__/examples/`example'.do"
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
    if ($exitCode -ne 0) {
        throw "Stata install-check process failed with exit code $exitCode. Log: $logFile"
    }

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
