// tests.do — master test runner for RdistNN Stata package
// Inspired by btskinner.io/code/stata-unit-testing/
//
// Run from the repository root:
//   do "tests/tests.do"
//
// The runner sources the framework once, pre-loads the Mata library, then
// executes each test-*.do file in order.  Each file updates the shared
// globals $tf_pass / $tf_fail.  A grand total is printed at the end.
// Exit code 1 is returned if any test failed (useful for CI).

version 17
set more off

// ── adopath: give src/ priority so commands resolve to the local package ──
local pkg_root : pwd
adopath + "`pkg_root'/src"

// ── load framework (defines test_assert, test_section, DGP programs) ──────
do "tests/test_framework.do"

// ── pre-load Mata library (the per-file guard prevents double-compile) ─────
quietly do "src/_rdistnn_mata.ado"

// ── ordered list of test files ────────────────────────────────────────────
local test_files ///
    test-mata-distances      ///
    test-mata-weights-exact  ///
    test-mata-weights-asymp  ///
    test-mata-jackknife      ///
    test-mata-jackknife-vector ///
    test-weights-dnn-smoke   ///
    test-weights-dnn-options ///
    test-weights-dnn-errors  ///
    test-np-dnn-smoke        ///
    test-np-dnn-methods      ///
    test-np-dnn-errors       ///
    test-cate-dnn-core       ///
    test-cate-dnn-nuisance   ///
    test-cate-dnn-errors     ///
    test-plqr-dnn-core       ///
    test-plqr-dnn-methods    ///
    test-plqr-dnn-errors     ///
    test-integration

// ── run each file ─────────────────────────────────────────────────────────
foreach f of local test_files {
    di _n as text "══ `f'.do " _dup(45) "═"
    cap noisily do "tests/`f'.do"
    if _rc {
        di as error "  [runner] `f'.do aborted unexpectedly (rc=`_rc')"
        global tf_fail = $tf_fail + 1
    }
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
