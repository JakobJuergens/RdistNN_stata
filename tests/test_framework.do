// test_framework.do
// Minimalist unit-test harness for RdistNN Stata package.
// Inspired by btskinner.io/code/stata-unit-testing/
//
// Source this file at the top of every test-*.do file (or let tests.do
// load it once for the session).  Programs are guarded with cap program drop
// so re-sourcing is safe.
//
// Usage:
//   test_assert "descriptive message" condition_expression
//   test_section "label"
//
// Counters:
//   $tf_pass / $tf_fail accumulate across all files in a session.
//   tests.do reports the grand total.

version 17

// ── global counters (safe to initialise multiple times) ───────────────────
if "$tf_pass" == "" {
    global tf_pass 0
}
if "$tf_fail" == "" {
    global tf_fail 0
}

// ── test_assert ────────────────────────────────────────────────────────────
// Usage: test_assert "message" condition
//   condition  : any expression valid inside `assert'.  Expand r()/e() scalars
//                into locals *before* the call so they are baked into the
//                condition string (r() is cleared inside the program call).
//   Example:
//     weights_dnn x1 x2, x0(0 0) k(5) generate(w)
//     quietly sum w
//     local wsum = r(sum)
//     test_assert "[smoke] weights sum to 1" abs(`wsum' - 1) < 1e-8
cap program drop test_assert
program define test_assert
    gettoken msg cond : 0, quotes
    cap assert `cond'
    if _rc == 0 {
        di as result "    PASS " as text `"`msg'"'
        global tf_pass = $tf_pass + 1
    }
    else {
        di as error "    FAIL " as text `"`msg'"'
        global tf_fail = $tf_fail + 1
    }
end

// ── test_section ───────────────────────────────────────────────────────────
// Prints a section banner for visual grouping in the results window.
cap program drop test_section
program define test_section
    local label `"`0'"'
    di _n as text `"  -- `label' --"'
end

// ── DGP programs ──────────────────────────────────────────────────────────
// Shared synthetic datasets used across test files.

// rdnn_dgp_basic: bivariate X, outcome y, binary D, oracle nuisances m0/m1/e
cap program drop rdnn_dgp_basic
program define rdnn_dgp_basic
    syntax [, n(integer 200) seed(integer 42)]
    clear
    set obs `n'
    set seed `seed'
    gen double x1 = rnormal()
    gen double x2 = rnormal()
    gen double y  = sin(x1) + 0.5 * x2 + rnormal() * 0.3
    gen double e  = invlogit(0.5 * x1 - 0.3 * x2)
    gen byte   D  = (runiform() < e)
    gen double m0 = 0.5 * x2
    gen double m1 = sin(x1) + 0.5 * x2
end

// rdnn_dgp_np: bivariate X, continuous outcome only (no treatment)
cap program drop rdnn_dgp_np
program define rdnn_dgp_np
    syntax [, n(integer 150) seed(integer 7)]
    clear
    set obs `n'
    set seed `seed'
    gen double x1 = rnormal()
    gen double x2 = rnormal()
    gen double y  = x1 - 0.5 * x2 + rnormal() * 0.25
end

// rdnn_dgp_const: constant outcome — jackknife SE should be 0
cap program drop rdnn_dgp_const
program define rdnn_dgp_const
    syntax [, n(integer 100) seed(integer 1) yval(real 3.5)]
    clear
    set obs `n'
    set seed `seed'
    gen double x1 = rnormal()
    gen double x2 = rnormal()
    gen double y  = `yval'
end

// rdnn_dgp_plqr: bivariate X, two PLQR regressors, oracle q_tau(X), and m_D(X)
cap program drop rdnn_dgp_plqr
program define rdnn_dgp_plqr, rclass
    syntax [, n(integer 140) seed(integer 31)]
    clear
    set obs `n'
    set seed `seed'
    gen double x1  = rnormal()
    gen double x2  = rnormal()
    gen double md1 = 0.30 * x1 - 0.20 * x2
    gen double md2 = -0.10 * x1 + 0.40 * x2
    gen double d1  = md1 + rnormal() * 0.40
    gen double d2  = md2 + rnormal() * 0.40
    gen double q   = 0.25 * x1 - 0.15 * x2
    gen double y   = d1 - 0.50 * d2 + q + rnormal() * 0.50
    return scalar beta1 = 1
    return scalar beta2 = -0.5
    return scalar tau   = 0.5
end
