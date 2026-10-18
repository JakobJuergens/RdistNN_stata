// tests.do — master test runner for RdistNN Stata package
// Inspired by btskinner.io/code/stata-unit-testing/
//
// Run from the repository root:
//   do "tests/tests.do"
//   do "tests/tests.do" "test-np-dnn-*"
//
// The runner sources the framework once, pre-loads the Mata library, then
// executes each test-*.do file in order.  Each file updates the shared
// globals $tf_pass / $tf_fail.  A grand total is printed at the end.
// Exit code 1 is returned if any test failed (useful for CI).

version 17
set more off
args filter
if `"`filter'"' == "" local filter "*"

// Every master invocation starts fresh, including focused and repeated runs.
global tf_pass 0
global tf_fail 0

// ── adopath: give src/ priority so commands resolve to the local package ──
local pkg_root : pwd
local pkg_root = subinstr("`pkg_root'", "\", "/", .)
adopath ++ "`pkg_root'/src"

// ── load framework (defines test_assert, test_section, DGP programs) ──────
do "tests/test_framework.do"

// ── pre-load Mata library (the per-file guard prevents double-compile) ─────
quietly do "src/_rdistnn_mata.ado"

// ── ordered list of test files ────────────────────────────────────────────
local test_files ///
    test-plot-command-order ///
    test-mata-distances      ///
    test-mata-rank-r-mirror  ///
    test-mata-weights-exact  ///
    test-mata-weights-asymp  ///
    test-mata-jackknife      ///
    test-mata-jackknife-golden ///
    test-mata-jackknife-vector ///
    test-r-parity-priority   ///
    test-weights-dnn-smoke   ///
    test-weights-dnn-options ///
    test-weights-dnn-errors  ///
    test-weights-dnn-plugin  ///
    test-weights-dnn-r-invariants ///
    test-np-dnn-smoke        ///
    test-np-dnn-methods      ///
    test-np-dnn-errors       ///
    test-np-dnn-r-identities ///
    test-np-dnn-parity       ///
    test-cate-dnn-core       ///
    test-cate-dnn-nuisance   ///
    test-cate-dnn-score-inputs ///
    test-cate-dnn-batch      ///
    test-cate-dnn-tuning     ///
    test-cate-dnn-diagnostics ///
    test-cate-dnn-postestimation ///
    test-cate-dnn-naming     ///
    test-cate-dnn-errors     ///
    test-cate-dnn-r-identities ///
    test-clate-mata          ///
    test-clate-dnn-core      ///
    test-clate-dnn-errors    ///
    test-clate-dnn-tuning    ///
    test-clate-r-parity      ///
    test-clate-dnn-postestimation ///
    test-estimator-reporting ///
    test-estimator-diagnostics ///
    test-estimator-plots      ///
    test-premerge-clate       ///
    test-premerge-postestimation ///
    test-plqr-dnn-core       ///
    test-plqr-dnn-methods    ///
    test-plqr-dnn-errors     ///
    test-plqr-dnn-r-identities ///
    test-integration

// ── run each file ─────────────────────────────────────────────────────────
local selected 0
foreach f of local test_files {
    if !strmatch("`f'", "`filter'") continue
    local ++selected
    di _n as text "══ `f'.do " _dup(45) "═"
    cap noisily do "tests/`f'.do"
    local run_rc = _rc
    if `run_rc' {
        di as error "  [runner] `f'.do aborted unexpectedly (rc=`run_rc')"
        global tf_fail = $tf_fail + 1
    }
}

if `selected' == 0 {
    di as error "No test files match `filter'."
    exit 198
}

// ── grand total ───────────────────────────────────────────────────────────
local total = $tf_pass + $tf_fail
di ""
di as text _dup(55) "─"
if $tf_fail == 0 {
    di as result "All $tf_pass tests passed."
}
else {
    di as error "$tf_fail of `total' tests FAILED."
}
di as text _dup(55) "─"

if $tf_fail > 0 exit 1
