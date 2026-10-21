// test-weights-dnn-plugin.do
// Tests for public weights_dnn using the native distance engine.
//
// Tags: [plugin] [parity]

quietly do "src/_rdistnn_mata.ado"
do "tests/test_framework.do"

cap program drop _rdnn_compare_plugin_weights
program define _rdnn_compare_plugin_weights
    syntax varlist(min=1) [if] [in] , ///
        x0(numlist)                    ///
        Method(string)                 ///
        k(integer)                     ///
        s(integer)                     ///
        Tol(real)                      ///
        Label(string)

    cap drop __w_mata __w_plugin __w_diff

    global rdistnn_use_plugin 0
    weights_dnn `varlist' `if' `in', ///
        x0(`x0') k(`k') s(`s') method(`method') generate(__w_mata)

    global rdistnn_use_plugin 1
    weights_dnn `varlist' `if' `in', ///
        x0(`x0') k(`k') s(`s') method(`method') generate(__w_plugin)

    quietly count if missing(__w_mata) != missing(__w_plugin)
    local miss_mismatch = r(N)

    gen double __w_diff = abs(__w_mata - __w_plugin) ///
        if !missing(__w_mata, __w_plugin)
    quietly summarize __w_diff, meanonly
    local ncomp = r(N)
    local maxdiff = r(max)

    test_assert `"[plugin parity] `label': compared observations exist"' ///
        `ncomp' > 0
    test_assert `"[plugin parity] `label': missing pattern matches"' ///
        `miss_mismatch' == 0
    test_assert `"[plugin parity] `label': max abs diff <= `tol'"' ///
        `maxdiff' <= `tol'

    cap drop __w_mata __w_plugin __w_diff
end

test_assert "[plugin] rdistnn.plugin is discoverable" "$rdistnn_plugin_ok"=="1"

test_section "weights_dnn: [plugin] signed x0 parity"
// ════════════════════════════════════════════════════════════

rdnn_dgp_np, n(140) seed(702)

local scenario1 "0 0"
local scenario2 "0.2 0.3"
local scenario3 "0.2 -0.3"
local scenario4 "-0.2 0.3"
local scenario5 "-0.2 -0.3"

foreach method in exact asymp {
    if "`method'" == "exact" local tol 1e-12
    if "`method'" == "asymp" local tol 1e-8

    forvalues j = 1/5 {
        local this_x0 "`scenario`j''"
        _rdnn_compare_plugin_weights x1 x2, ///
            x0(`this_x0') method(`method') k(5) s(15) tol(`tol') ///
            label("method=`method' x0=(`this_x0')")
    }
}

// ════════════════════════════════════════════════════════════
test_section "weights_dnn: [plugin] one-dimensional signed x0 parity"
// ════════════════════════════════════════════════════════════

rdnn_dgp_np, n(110) seed(703)

foreach method in exact asymp {
    if "`method'" == "exact" local tol 1e-12
    if "`method'" == "asymp" local tol 1e-8

    _rdnn_compare_plugin_weights x1, ///
        x0(-0.4) method(`method') k(4) s(12) tol(`tol') ///
        label("method=`method' one-dimensional x0=(-0.4)")
}

// ════════════════════════════════════════════════════════════
test_section "weights_dnn: [plugin] if/in and markout parity"
// ════════════════════════════════════════════════════════════

rdnn_dgp_np, n(150) seed(704)

_rdnn_compare_plugin_weights x1 x2 if x1 > -0.25 in 5/90, ///
    x0(0.2 -0.3) method(exact) k(4) s(13) tol(1e-12) ///
    label("restricted if/in x0=(0.2 -0.3)")

rdnn_dgp_np, n(150) seed(705)
replace x2 = . in 3/7
replace x1 = . in 11

_rdnn_compare_plugin_weights x1 x2, ///
    x0(0.2 -0.3) method(exact) k(4) s(13) tol(1e-12) ///
    label("missing covariates excluded by markout")

capture macro drop rdistnn_use_plugin
